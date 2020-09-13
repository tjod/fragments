#MIT License
# 
#Copyright (c) 2020 TJ O'Donnell
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
import sys
import sqlite3
import pickle

#from rdkit import Chem
#from rdkit.Chem import rdMolDescriptors

def addView(cur, add):
	# Create view of extra columns' property_values
	property_names = {}
	if len(add) == 0:
		cur.execute("Select propid, propname From property_names")
		for row in cur:
			property_names[row[0]] = row[1]
	else:
		q = ["?"] * len(add)
		sql = "Select propid, propname From property_names Where propname in (%s)" % ",".join(q)
		cur.execute(sql, add)
		for row in cur:
			property_names[row[0]] = row[1]
	if len(property_names) == 0:
		cur.execute("Create Temporary View mol_properties As Select molid from molecule")
	else:
		sql = "Create Temporary View mol_properties As With rows As (%s) Select molid, %s From rows Group By molid"
		cases = []
		for (propid, propname) in property_names.items():
			cases.append("Case When propid=%d then propvalue End As p%d" % (propid, propid))
		case_clause = "Select molid, %s From property_values Join property_names Using (propid)" % ",".join(cases)
		select = []
		for (propid, propname) in property_names.items():
			select.append('Group_concat(p%d) As "%s"' % (propid, propname))
		select_clause = ",".join(select)
		cur.execute(sql % (case_clause, select_clause))
	return property_names
	
def db_predict(cur, property_name, format, add, compare):
	if add: property_names = addView(cur, add[0])
	result = []
	cur.execute("Create Temporary View If Not Exists counts As Select molid, atomid, Count(atom_index) pcount From atoms Group By molid, atomid")
	if format == "sdf":
		sql = """With atmp As (Select coefficient As intercept From atomid_coefficients Where atomid='Intercept'),
		 properties As (select molid, Group_concat(printf('> <%s>%s%s%s', propname, char(10), propvalue, char(10)), char(10) ) pvals
		   From property_values Join property_names Using (propid) Group By molid),
		 result As (Select molid, molblock, atomid, pcount, coefficient
		   From molecule Join counts Using (molid) Join atomid_coefficients Using (atomid))
		 Select printf('%s%s%s> <%s>%s%.2f%s%s$$$$',
		   molblock,pvals,char(10),?,char(10),Sum(pcount*coefficient)+intercept,char(10),char(10))
		   From atmp Join result Join properties Using (molid) Group By molid,pvals Order By molid"""
		cur.execute(sql, [property_name])
		for row in cur:
			print (row[0])
	elif format == "tsv" or format == "csv":
		sep = "\t" if format == "tsv" else ","
		if add:
			sql = """With atmp As (Select coefficient As intercept From atomid_coefficients Where atomid='Intercept'),
				property As (Select * From mol_properties),
				result As (Select molid, atomid, pcount, coefficient
				  From molecule Join counts Using (molid) Join atomid_coefficients Using (atomid))
				Select molid, printf('%.2f', Sum(pcount*coefficient)+intercept) pval, property.*
				  From atmp Join result Left Join property Using (molid) Group By molid Order By molid"""
			cur.execute(sql)
			header = ["molid", property_name]
			for d in cur.description[3:]:
				header.append(d[0])
			print (sep.join(header))
			for row in cur:
				prow = [str(p) for p in row]
				del prow[2] # dupliate molid from mol_properties.*
				print (sep.join(prow))
		else:
			sql = """With atmp As (Select coefficient As intercept From atomid_coefficients Where atomid='Intercept'),
				result As (Select molid, atomid, pcount, coefficient
				  From molecule Join counts Using (molid) Join atomid_coefficients Using (atomid))
				Select molid, printf('%.2f', Sum(pcount*coefficient)+intercept) pval
				  From atmp Join result Group By molid Order By molid"""
			cur.execute(sql)
			print (sep.join(["molid",property_name]))
			for row in cur:
				prow = [str(p) for p in row]
				print (sep.join(prow))
			
def pickle_predict(cur, pickle_file, property_name, format, add, compare):
	model = pickle.load(open(pickle_file, 'rb'))
# 	intercept = model.intercept_
# 	coefficients = model.coef_
# 	print (intercept, len(coefficients))
	cur.execute("Create Temporary View If Not Exists counts As Select molid, atomid, Count(atom_index) pcount From atoms Group By molid, atomid")
	cur.execute("Select count(molid) From molecule")
	nmols = cur.fetchone()[0]
	cur.execute("Select count(atomid) From atomid_coefficients Where atomid != 'Intercept'")
	natomid = cur.fetchone()[0]
	#sql = "Select molid, Coalesce(pcount,0) pcount From molecule Left Join counts Using (molid) Join atomid_coefficients Using (atomid) Order By molid, atomid"
	sql = """With matrix As (Select atomid, molid From atomid_coefficients Join molecule Where atomid != 'Intercept')
	 Select Coalesce(pcount,0) pcount From matrix Left Join counts Using (molid,atomid) Order By molid, atomid"""
	cur.execute(sql)
	counts = []
	for i in range(0, nmols):
		counts.append([row[0] for row in cur.fetchmany(natomid)])
	#print (len(counts), len(counts[0]))
	#print (counts)
	property_pred = model.predict(counts)
	if format == "sdf":
		sql = """With properties As (Select molid,
		   Group_concat( printf('> <%s>%s%s%s', propname, char(10), propvalue, char(10)), char(10) ) pvals
		   From property_values Join property_names Using (propid) Group By molid)
		 Select molid, molblock, pvals From molecule Join properties Using (molid) Order By molid"""
		cur.execute(sql)
		for row in cur:
			imol = int(row[0]) - 1
			print ("%s%s\n> <%s>\n%.2f\n\n$$$$" % (row[1], row[2], property_name, property_pred[imol]))
	elif format == "tsv" or format == "csv":
		sep = "\t" if format == "tsv" else ","
		print (sep.join(["molid", property_name]))
		i = 0
		for p in property_pred:
			i += 1
			print ("%d%s%.2f" % (i,sep,p))
	
def getMolMolblock(cur):
	cur.execute("Select molblock from main.molecule")
	mol = cur.fetchone()[0]
	return mol

def list_properties(cur):
	property_names = []
	cur.execute("Select propname From property_names")
	for row in cur:
		property_names.append (row[0])
	return property_names

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Create a model that fits molecular properties to atomid counts",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("database", help="input sqlite3 (db3) file from sdf2sql.py")
    parser.add_argument("property", help="name (sdf tag) of predicted property to output")
    parser.add_argument("model", help="input model sqlite3 (db3) file")
#     parser.add_argument("ndescriptors", type=int, help="number of atomid descriptors to consider")
#     parser.add_argument("level", type=int, help="maximum atomid level to consider")
    parser.add_argument("-p", "--pickle", help="input python pickle of model", default=None)
    parser.add_argument("-l", "--list", help="list properties in input db3 file, and exit", action="store_true")
    parser.add_argument("-f", "--format", help="output format", choices=["sdf", "csv", "tsv"], default="tsv")
    parser.add_argument("-a", "--add", help="output values in these sd tags; --add alone adds all", action="append", nargs="*")
    parser.add_argument("-c", "--compare", help="compare to values in this sd tag", default=None)
#     parser.add_argument("-c", "--correlated", help="keep or remove correlated atomid before fitting", choices=["remove", "keep"], default="remove")
#     parser.add_argument("-n", "--modulo_N", help="include only every nth molecule; useful for shorter tests", type=int, default=1)
#     parser.add_argument("-p", "--plot", help="output file for plot of input vs predicted values", default=None)
#     parser.add_argument("-v", "--verbosity", type=int, help="increase output verbosity", default=1)
    return parser
   
def main():
	parser = parse_args()
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	parsed = parser.parse_args()
	
	mol_db = parsed.database
	model_db = parsed.model
	pickle_file = parsed.pickle
	property_name = parsed.property
	format = parsed.format
	compare = parsed.compare
	add = parsed.add
	list = parsed.list
	con = sqlite3.connect(mol_db)
	cur = con.cursor()
	if list:
		[print(p) for p in list_properties(cur)]
		sys.exit()
	cur.execute('Attach ? As model',  [model_db])
	if pickle_file:
		pickle_predict(cur, pickle_file, property_name, format, add, compare)
	else:
		db_predict(cur, property_name, format, add, compare)

if __name__ == "__main__":
    #import cProfile
    #cProfile.run('main()', 'profile.out')
    main()
