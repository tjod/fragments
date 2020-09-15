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
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

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

def predict(cur, property_name, pickle_file):
	predicted_values = []
	cur.execute("Create Temporary View If Not Exists counts As Select molid, atomid, Count(atom_index) pcount From atoms Group By molid, atomid")
	if pickle_file:
		model = pickle.load(open(pickle_file, 'rb'))
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
		predicted_values = model.predict(counts)
		cur.execute("Create Temporary Table new_property (molid Integer, propvlaue Numeric)")
		for imol in range(0, len(predicted_values)):
			cur.execute("Insert Into new_property (molid, propvlaue) Values (?,?)", (imol+1, predicted_values[imol]))
	else:
		sql = """Create Temporary View new_property As With
				atmp As (Select coefficient As intercept From atomid_coefficients Where atomid='Intercept'),
				result As (Select molid, atomid, pcount, coefficient
				  From molecule Join counts Using (molid) Join atomid_coefficients Using (atomid))
				Select molid, Sum(pcount*coefficient)+intercept propvlaue
				  From atmp Join result Group By molid Order By molid"""
		cur.execute(sql)
		cur.execute("Select propvlaue From new_property Order By molid")
		predicted_values = cur.fetchall() # list(cur.fetchall())
	return predicted_values

def output_file(cur, property_name, fpout, format, add):
	# output requested property values
	if add: property_names = addView(cur, add[0])
	cur.execute("Select * from new_property Order By molid")
	if format == "sdf":
		if add:
			sql = "Select molblock, mol_properties.*, new_property.* From molecule Join mol_properties Using (molid) Join new_property Using (molid)"
		else:
			sql = "Select molblock, new_property.* From molecule Join new_property Using (molid)"
		cur.execute(sql)
		header = [property_name if d[0] == "propvlaue" else d[0] for d in cur.description]
		for row in cur:
			for icol in range(0, len(row)):
				p = row[icol]
				pname = header[icol]
				if pname == "molid":
					pass
				elif pname == "molblock":
					print (p, file=fpout)
				else:
					fmt = "> <%s>\n"
					fmt += "%d" if type(p) is int else "%.2f" if type(p) is float else "%s"
					fmt += "\n"
					print (fmt % (pname, p), file=fpout)
			print ("$$$$", file=fpout)
	elif format == "tsv" or format == "csv":
		sep = "\t" if format == "tsv" else ","
		if add:
			sql = "Select * from mol_properties Join new_property Using (molid) Order By molid"
		else:
			sql = "Select * From new_property Order By molid"
		cur.execute(sql)	
		header = [property_name if d[0] == "propvlaue" else d[0] for d in cur.description]
		print (sep.join(header), file=fpout)
		for row in cur:
			prow = []
			for p in row:
				fmt = "%d" if type(p) is int else "%.2f" if type(p) is float else "%s"
				prow.append(fmt % p)
			print (sep.join(prow), file=fpout)

def list_properties(cur):
	property_names = []
	cur.execute("Select propname From property_names")
	for row in cur:
		property_names.append (row[0])
	return property_names

def get_property_values(cur, compared_property_name):
	cur.execute("""Select molid, propvalue From property_values Join property_names Using (propid)
	 Where propname = ? Order By molid""", [compared_property_name])
	property_values = [float(row[1]) for row in cur]
	return property_values
	return 

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Create a model that fits molecular properties to atomid counts",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("database", help="input sqlite3 (db3) file from sdf2sql.py")
    parser.add_argument("property", help="name (sdf tag) of predicted property to output")
    parser.add_argument("model", help="input model sqlite3 (db3) file")
    parser.add_argument("out", help="output file name")
    parser.add_argument("-f", "--format", help="output format", choices=["sdf", "csv", "tsv"], default="tsv")
    parser.add_argument("-p", "--pickle", help="input python pickle of model", default=None)
    parser.add_argument("-a", "--add", help="output values in these sd tags; --add alone adds all", action="append", nargs="*")
    parser.add_argument("-c", "--compare", help="compare to values in this sd tag", default=None)
    parser.add_argument("-g", "--graph", help="output comparison graph/plot to file; requires --compare", default=None)
    parser.add_argument("-l", "--list", help="list properties in input db3 file, and exit", action="store_true")
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
	compare_tag = parsed.compare
	add = parsed.add
	list = parsed.list
	plotfile = parsed.graph
	if plotfile and not compare_tag:
		parser.print_help()
		parser.error("--graph requires --compare")
		sys.exit()
	fpout = open(parsed.out, "w")
	con = sqlite3.connect(mol_db)
	cur = con.cursor()
	if list:
		[print(p) for p in list_properties(cur)]
		sys.exit()
	cur.execute('Attach ? As model',  [model_db])
	predicted_values = predict(cur, property_name, pickle_file)
	output_file(cur, property_name, fpout, format, add)
	if compare_tag:
		property_values = get_property_values(cur, compare_tag)
		print ("%s:%s R-squared: %.3f" % (property_name, compare_tag, r2_score(property_values, predicted_values)))
		if plotfile:
			fig, ax = plt.subplots()
			ax.scatter(property_values, predicted_values, marker='.', s=16)
			plt.title("%s / %s" % (mol_db, model_db))
			plt.xlabel(compare_tag)
			plt.ylabel(property_name)
			fig.savefig(plotfile)		
	
if __name__ == "__main__":
    #import cProfile
    #cProfile.run('main()', 'profile.out')
    main()
