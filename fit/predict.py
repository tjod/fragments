import sys
import sqlite3

#from rdkit import Chem
#from rdkit.Chem import rdMolDescriptors

def predict(cur, property_name):
	result = []
	cur.execute("Create Temporary View counts As Select molid, atomid, Count(atom_index) pcount From atoms Group By molid, atomid")
	cur.execute("With atmp As (Select coefficient As intercept From atomid_coefficients Where atomid='Intercept'), properties As (select molid, Group_concat(printf('> <%s>%s%s%s', propname, char(10), propvalue, char(10)), char(10) ) pvals from property_values join property_names using (propid) Group By molid), result As (Select molid, molblock, atomid, pcount, coefficient From molecule Join counts Using (molid) Join atomid_coefficients Using (atomid)) Select printf('%s%s%s> <%s>%s%.2f%s%s$$$$', molblock,pvals,char(10),?,char(10),Sum(pcount*coefficient)+intercept,char(10),char(10)) From atmp Join result Join properties Using (molid) Group By molid,pvals Order By molid", [property_name])
	for row in cur:
		print (row[0])

def getMolMolblock(cur):
	cur.execute("Select molblock from main.molecule")
	mol = cur.fetchone()[0]
	return mol

mol_db = sys.argv[1]
model_db = sys.argv[2]
property_name = sys.argv[3]
con = sqlite3.connect(mol_db)
cur = con.cursor()
cur.execute('Attach ? As model',  [model_db])
predict(cur, property_name)
