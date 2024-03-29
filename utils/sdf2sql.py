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
#from __future__ import print_function
import sys
import re
import sqlite3
import os
#from io import StringIO
from os import path


from rdkit import Chem


BOND = [None, "-", "=", "#", ":"]
PERIODIC_TABLE = Chem.GetPeriodicTable()
def makeTables(cursor):
    # make tables
    cursor.execute("Create Table If Not Exists property_names (propid integer primary key autoincrement, propname text Unique)")
    cursor.execute("Create Table If Not Exists molecule (molid integer primary key autoincrement, molname text, smiles text, molblock text)")
    cursor.execute("Create Table If Not Exists property_values (molid integer, propid integer, propvalue text)")
    # training set and test set each with half
    #cursor.execute("Create View training_set as select * from molecule where molid%2=0")
    #cursor.execute("Create View test_set as select * from molecule where molid%2>0")
    # training set with all molecules
    #cursor.execute("Create View training_set as select * from molecule")
    # test set with one molecule just to keep R scripts happy
    #cursor.execute("Create View test_set as select * from molecule where molid<100")

def parentage(con, atomid):
    cur = con.cursor()
    result = []
    sql = """With Recursive parent_of(atomid, parentid, iteration) As (
              Select atomid, parentid, iteration From parents Where atomid = ?
             Union
              Select parents.atomid, parents.parentid, parents.iteration From parents, parent_of Where parents.atomid = parent_of.parentid
             )
             Select atomid, parentid, iteration, smarts From parent_of Left Join smarts Using (atomid) Order By iteration Desc"""
    cur.execute(sql, [atomid])
    for row in cur:
        result.append({"atomid":row[0], "parentid":row[1], "iteration":row[2], "smarts":row[3]})
    return result

def addMol(cursor, mol, molblock, nmol):
    if mol:
        try:
            molname = mol.GetProp("_Name")
        except  Exception as e:
            molname = str(e)
            mol.SetProp("_Name", molname)
            print ("Warning processing name of molecule #%d: %s" % (nmol, molname))
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        molname = "Error"
        smiles = None

    cursor.execute("Insert Into molecule (molname, smiles, molblock) Values (?,?,?)", [molname, smiles, molblock])
    return cursor.lastrowid

def addProp(cursor, imol, prop, val):
    cursor.execute("Select propid From property_names Where propname = (?)", [prop])
    row = cursor.fetchone()
    if row:
        propid = row[0]
    else:
        cursor.execute("Insert Into property_names (propname) Values (?)", [prop])
        propid = cursor.lastrowid
    cursor.execute("Insert Into property_values (molid, propid, propvalue) Values (?,?,?)", [imol, propid, val])

def processSDF(cursor, filename, split, storeMolblock):
    nmol = 0
    # capture stderr when processing mol
    #sio = sys.stderr = StringIO()
    #fp = open(filename, 'rb')
    #suppl = Chem.ForwardSDMolSupplier(fp)
    
    if filename.endswith(".gz"):
        import gzip
        inf = gzip.open(filename)
        suppl = Chem.ForwardSDMolSupplier(inf)
        if storeMolblock:
            print ("Warning: cannot store molblock from gzip file", file=sys.stderr)
    else:
        suppl = Chem.SDMolSupplier(filename)

    for mol in suppl:
        if hasattr(suppl, "GetItemText"):
            (molblock, sep, moldata) = suppl.GetItemText(nmol).partition('M  END')
            if storeMolblock:
                molstore = molblock+sep
            else:
                molstore = None
        else:
            molstore = None
        nmol += 1
        imol = addMol(cursor, mol, molstore , nmol)
        if mol:
            for p in mol.GetPropNames():
                if split:
                    for sp in mol.GetProp(p).split(","):
                        addProp(cursor, imol, p, sp)
                else:
                    addProp(cursor, imol, p, mol.GetProp(p))
        else:
            print ("Error adding molecule #%d" % nmol, file=sys.stderr)
            # molname stores the stderr when processing the mol
            #cursor.execute("Update molecule Set molname = ? Where molid = ?", [sio.getvalue(), imol])
            #sio = sys.stderr = StringIO() # reset

    return nmol

def atomSet(ofile, atom_add, molid, iteration, atomid, atom_index, atom_attr):
    atom = (molid, iteration, atomid, atom_index+atom_add, atom_attr)
    atoms = [atom]
    for line in ofile:
        got = re.match("^([0-9]*) \\[(.*)\\] *(.*)$", line)
        if got:
            atom_index = int(got.group(1))
            atom_attr = got.group(2).split(",")
            atomid = got.group(3)
            atom = (molid, iteration, atomid, atom_index+atom_add, atom_attr)
            atoms.append(atom)
        else:
            break
    return atoms

def makeSMARTS(cursor, iteration, atom, smarts_dict):
    smarts = ""
    symbol = ""
    (molid, iter, atomid, atom_index, attr) = atom
    if iteration == 0:
        atom_spec = {
            "in_ring":"" if int(attr[0]) == 1 else "!",
            "aromatic":"a" if int(attr[1]) == 1 else "A",
             "degree":int(attr[2]),
            # "valence":int(attr[3]),
            # "valence":int(attr[3]) + int(attr[4]),
            "hcount":int(attr[4]),
            "charge":int(attr[5]),
            "atnum":int(attr[6]),
            "mass":int(attr[7])
        }
        # print (atom_spec)
        symbol = Chem.PeriodicTable.GetElementSymbol(PERIODIC_TABLE, atom_spec["atnum"])
        if atom_spec["aromatic"] == "a" : symbol = symbol.lower()
        smarts = ("[#{atnum};{aromatic};{in_ring}R;D{degree};H{hcount};{charge:+d}]").format(**atom_spec)
        cursor.execute ("Insert Or Ignore Into smarts (atomid, smarts, symbol, in_ring, is_aromatic, hvy_degree, hvy_valence, hcount, formal_charge, atomic_number, mass) Values (?,?,?,?,?,?,?,?,?,?,?)", [atomid, smarts, symbol] + attr)
    elif iteration == 1:
        i=0
        iter = attr[i]
        root_atomid = attr[i+1].strip()
        smarts = smarts_dict[root_atomid]["smarts"]
        symbol = smarts_dict[root_atomid]["symbol"]
        i += 2
        while i < len(attr):
            bond_order = int(attr[i])
            bonded_atomid = attr[i+1].strip()
            #print (i, bond_order, bonded_atomid)
            smarts += "(%s%s)" % (BOND[bond_order], smarts_dict[bonded_atomid]["smarts"])
            symbol += "(%s%s)" % (BOND[bond_order], smarts_dict[bonded_atomid]["symbol"])
            i += 2
        cursor.execute ("Insert Or Ignore Into smarts (atomid, smarts, symbol) Values (?,?,?)", (atomid, smarts, symbol))
#     else:
#         atomid = None
#         smarts = None

    return (atomid, {"smarts": smarts, "symbol": symbol})

def add_parent(cursor, iteration, atom):
    atomid = atom[2].strip()
    atom_attr = atom[4]
    parentid = None if iteration == 0 else atom_attr[1].strip()
    cursor.execute("Insert Or Ignore Into parents (atomid, parentid, iteration) Values (?,?,?)", (atomid, parentid, iteration))

def processOutfile(cursor, outfile):
    if outfile.close:
        # already opened file or StringIO
        ofile = outfile
    else:
        try:
            if outfile.endswith(".gz"):
                import gzip
                ofile = gzip.open(outfile, "rt")
            else:
                ofile = open(outfile, "r")
        except:
            return 0

    molid = 0
    iteration = 0
    cursor.execute("Drop Table If Exists atoms")
    cursor.execute("Drop Table If Exists smarts")
    cursor.execute("Drop Table If Exists parents")
    cursor.execute("Drop View If Exists atom_types")
    cursor.execute("Create Table atoms (molid integer, iteration integer, atomid text, atom_index integer)")
    cursor.execute("Create Table smarts (atomid text Unique, smarts Text, symbol Text, in_ring integer, is_aromatic integer, hvy_degree integer, hvy_valence integer, hcount integer, formal_charge integer, atomic_number integer, mass integer)")
    cursor.execute("Create Table parents (atomid text Unique, parentid text, iteration integer)")
    #cursor.execute("Create Table counts (molid integer, iteration integer, propvalue text, pcount integer)")
    cursor.execute("Create View atom_types As Select atomid, iteration, /*Count(distinct iteration) Niterations,*/ Count(distinct molid) frequency From atoms Group By atomid")
    smarts_dict = {}
    for line in ofile:
        #fp = re.search("fp:.*{(.*)}", line)
        #if fp:
        #    #print (iteration, line)
        #    for f in fp.group(1).split(","):
        #        (atomid, count) = f.replace(" ","").split("=")
        #        cursor.execute ("Insert Into counts (molid,iteration,propvalue,pcount) values (?,?,?,?)" , (molid, iteration-1, atomid, int(count)))
        got = re.match("^([0-9]*) \\[(.*)\\] *(.*)$", line)
        if got:
            atom_index = int(got.group(1))
            atom_attr = got.group(2).split(",")
            atomid = got.group(3)
            if atom_index == 1:
                molid += 1
                atom_add = 0
                iteration = 0
            if atom_index == 0: atom_add = 1 # first set of indices are reported using atomid starting with 0
            for atom in atomSet(ofile, atom_add, molid, iteration, atomid, atom_index, atom_attr):
                cursor.execute ("Insert Into atoms (molid, iteration, atomid, atom_index) Values (?,?,?,?)",  atom[:4])
                (atom_id, smarts) = makeSMARTS(cursor, iteration, atom, smarts_dict)
                add_parent(cursor, iteration, atom)
                if atom_id: smarts_dict[atom_id] = smarts
            iteration += 1
    return molid

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Create sqlite (db3) file from input SD file and/or sdfCFP verbose output file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("db", help="output sqlite3 file")
    parser.add_argument("-v", "--verbosity", type=int, help="increase output verbosity", default=0)
    parser.add_argument("--sdf", help="input sd file, or - for stdin")
    parser.add_argument("-m", "--molblock", help="store molblocks in output db", action="store_true")
    parser.add_argument("--cfp", help="verbose output file from sdfCFP, or - for stdin")
    parser.add_argument("-k", "--keepH", help="keep explicit H atoms in sd file", action="store_true")
    parser.add_argument("-s", "--split", help="split property values on comma", action="store_true")

    return parser

def main():
    parser = parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    parsed = parser.parse_args()
    
    db = parsed.db
    sdf = parsed.sdf
    split = parsed.split
    keepH = parsed.keepH
    storeMolblock = parsed.molblock
    cfp = parsed.cfp

    if sdf:
        if sdf == "-" or path.exists(sdf):
            pass
        else:
            print ("Error: no file %s" % sdf, file=sys.stderr)
            sys.exit()   
    if cfp:
        if cfp == "-" or path.exists(cfp):
            pass
        else:
            print ("Error: no file %s" % cfp, file=sys.stderr)
            sys.exit()
    if sdf and cfp and sdf == "-" and cfp == "-":
        print ("sdf and cfp cannot both be stdin" % cfp, file=sys.stderr)

    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    makeTables(cursor)
    if sdf:
        Chem.WrapLogs() # to get error messages normally headed for stderr
        nprocessed = processSDF(cursor, sdf, split, storeMolblock)
        
    if cfp:
        nout = processOutfile(cursor, cfp)
        cursor.execute("Select count(molid) From molecule")
        nmol = cursor.fetchone()[0]
        if nmol == nout:
            print ("%d sdf molecules; %d in cfp file" % (nmol, nout))
        else:
            print ("Warning: %d sdf molecules, but %d in cfp file" % (nmol, nout))
        
    conn.commit()

if __name__ == "__main__":
    main()
     
