'''
Created on May 4, 2020

@author: tj o'donnell
'''
from __future__ import print_function
import openbabel as ob
import json
import re
import sqlite3
import sys

tbl = ob.OBElementTable()
bond_symbol = ["","-","=","#","~",":"]

def properties(atom):
    # return list of select atom properties
    return (atom.GetAtomicNum(), atom.GetFormalCharge(),
            #atom.ImplicitHydrogenCount() + atom.ExplicitHydrogenCount(), \
            atom.ImplicitHydrogenCount(), \
            atom.GetHvyValence(), \
            'a' if atom.IsAromatic() else 'A', atom.MemberOfRingCount())
    
# construct atom symbol, optionally using smarts with atom properties
def atom_type(atom, as_smarts):

    if as_smarts:
        atom_string = "[#%d;%+d;H%d;D%d;%s;R%d]" % properties(atom)
    else:
        atom_symbol = tbl.GetSymbol(atom.GetAtomicNum())
        atom_string = atom_symbol.lower() if atom.IsAromatic() else atom_symbol
    return atom_string

# get mol object from smiles
def MolFromSmiles(smi):
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    mol = ob.OBMol()
    obConversion.ReadString(mol, smi)
    return mol

# get neighbors of atom
def neighbors(atom):
    #return ob.OBAtomAtomIter(atom)
    #return [n for n in ob.OBAtomAtomIter(atom)]
    # will neighbors visited in sorted order make canonical smarts?
    return sorted([n for n in ob.OBAtomAtomIter(atom)], key=lambda a: hash(properties(a)))

# create tree of neighbors of atom up to depth
def branched_paths(atom, depth, root_atom, visited):
    paths = [{"atom": atom, "depth": depth}]
    visited.add(atom.GetIndex())
    if depth > 0:
        #print (json.dumps(neighbors(atom), cls=AtomEncoder))
        for nbr in neighbors(atom):
            #if nbr.IsHydrogen(): continue
            # don't branch back to root
            #if nbr == root_atom: continue
            # don't branch back to visited atoms
            if nbr.GetIndex() in visited: continue
            visited.add(nbr.GetIndex())
            nn = branched_paths(nbr, depth-1, atom, visited)[0]
            paths.append(nn)
    return (paths, visited)

# turn neighbour list (path) into smarts string
def path_to_string(atoms, with_properties, root_atom):
    s = ""
    atom0 = atoms[0]["atom"]
    atom_depth = atoms[0]["depth"]

    # lots of output choices
    if with_properties == 3:
        # all atoms
        as_smarts = True
    elif with_properties == 2:
        # tail atoms
        as_smarts = False if atom_depth > 0 else True
    elif with_properties == 1:
        # if root_atom is None, then this is Ur-root_atom, primary branch point atom
        as_smarts = False if root_atom else True
    elif with_properties == 4:
        # if root_atom is None, then this is Ur-root_atom, primary branch point atom
        as_smarts = False if (root_atom and atom_depth > 0) else True
    else:
        # no atoms
        as_smarts = False

    if root_atom:
        # has a root atom, so include bond symbol
        bond_order = 5 if root_atom.GetBond(atom0).IsAromatic() else root_atom.GetBond(atom0).GetBondOrder()
    else:
        # Ur-root_atom, primary branch point atom; so no bond symbol
        bond_order = 0
    
    s += "%s%s" % (bond_symbol[bond_order], atom_type(atom0, as_smarts) )    
    for atom in atoms[1:]:
        s += "(%s)" % path_to_string(atom, with_properties, atom0)
    return s

# allow json encoding of OBAtom
class AtomEncoder(json.JSONEncoder):
    def default(self, z):
        if isinstance(z, ob.OBAtom):
            return z.GetIndex()+1
        else:
            return super().default(z)

def show_results(results):
    # pick your output format(s)
    print ()
    print ("\n".join(results)) # for human consumption?
    print ()
    print (".".join(results)) # paste into MarvinSketch
    print ()
    print ("Values('" + "'),('".join(set(results)) + "')") # useful in SQL
# e.g. at https://www.gnova.com/misc/Site/chembl.html
# With amol As (Select rd.rdmol('c1ccccc1C(=O)NC')),
# s As (SELECT * FROM (
# 
#Values('[c;+0;D2;H1;a;R1](:c(:c))(:c(:c)(-C))'),('[C;+0;D1;H3;A;R0](-N(-C))'),('[N;+0;D2;H1;A;R0](-C(-c)(=O))(-C)'),('[c;+0;D2;H1;a;R1](:c(:c)(-C))(:c(:c))'),('[c;+0;D2;H1;a;R1](:c(:c))(:c(:c))'),('[O;+0;D1;H0;A;R0](=C(-c)(-N))'),('[C;+0;D3;H0;A;R0](-c(:c)(:c))(=O)(-N(-C))'),('[c;+0;D3;H0;a;R1](:c(:c))(:c(:c))(-C(=O)(-N))')
# 
# ) AS t (smarts)),
# r As (Select smarts, rd.cansmiles(rd.smarts_to_rdmol(smarts)) From s)
# Select smarts, cansmiles,
#  rd.count_matches(rdmol, smarts), rd.count_matches(rdmol, cansmiles),
#  replace(rd.svg(rdmol, 400, 400, rd.smarts_to_rdmol(smarts, true), false, false, false, true), 'svg:', '') Structure,
#  replace(rd.svg(rdmol, 400, 400, rd.smarts_to_rdmol(cansmiles, true), true, false, false, true), 'svg:', '') Structure2
# From r, amol

def addMol(cursor, mol, molblock, smiles, nmol):
    if mol:
        try:
            molname = mol.GetTitle()
        except  Exception as e:
            molname = str(e)
            mol.SetTitle(molname)
            print ("Warning processing name of molecule #%d: %s" % (nmol, molname))
    else:
        molname = "Error"
        smiles = None
    (mb, sep, moldata) = molblock.partition('M  END\n')
    cursor.execute("Insert Into molecule (molname, smiles, molblock) Values (?,?,?)", [molname, smiles, mb+sep])
    return cursor.lastrowid

def fragment_mol(mol, depth, with_properties):
    fragments = []
    #visits = []
    for atom in ob.OBMolAtomIter(mol):
        #print (atom.GetIndex() + 1)
        #if atom.IsHydrogen(): continue # and depth > 0: continue
        (path, visited) = branched_paths(atom, depth, atom, set())
        #print (json.dumps(path, cls=AtomEncoder, indent=2))
        if path:
            s = path_to_string(path, with_properties, None)
            fragments.append(s)
        else:
            print (mol.GetTitle(), depth, atom.GetIndex()+1)
    return fragments

def addProp(cursor, imol, prop, val):
    cursor.execute("Select propid From property_names Where propname = (?)", [prop])
    row = cursor.fetchone()
    if row:
        propid = row[0]
    else:
        cursor.execute("Insert Into property_names (propname) Values (?)", [prop])
        propid = cursor.lastrowid
    cursor.execute("Insert Into property_values (molid, propid, propvalue) Values (?,?,?)", [imol, propid, val])
    
def processSDF(cur, sdfile, max_depth, with_properties):
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("sdf", "smi")
    mol = ob.OBMol()
    smol = []
    tags = []
    nmol = 0
    with open(sdfile) as fp:
        for line in fp:
            if line.startswith('$$$$'):
                molblock = "".join(smol)
                conv.ReadString(mol, molblock)
                smi = conv.WriteString(mol)
                nmol += 1
                molid = addMol(cur, mol, molblock, smi, nmol)
                for tag in tags:
                    d = mol.GetData(tag)
                    #print (d.GetAttribute(), d.GetDataType())
                    addProp(cur, molid, d.GetAttribute(), d.GetValue())
                if not "MW" in tags: addProp(cur, molid, "MW", mol.GetMolWt())

                previous_fragments = None
                for depth in range(0, max_depth+1):
                    fragments = fragment_mol(mol, depth, with_properties)
                    for i in range(len(fragments)):
                        if previous_fragments and fragments[i] == previous_fragments[i]:
                            # atom has no neighbors past previous depth, e.g. water at depth > 0, CO at depth > 1, et. al.
                            ff = None
                        else:
                            ff = fragments[i]
                        cur.execute ("Insert Into atoms (molid, iteration, atomid, atom_index) Values (?,?,?,?)", (molid, depth, ff, i+1))
                    previous_fragments = fragments
                
                smol = []
                tags = []
            else:
                smol.append(line)
                if line.startswith("> <"):
                    m = re.match("^> <(.*)>", line)
                    tags.append(m.group(1))

def makeTables(cursor):
    # make tables
    cursor.execute("Create Table If Not Exists property_names (propid integer primary key autoincrement, propname text Unique)")
    cursor.execute("Create Table If Not Exists molecule (molid integer primary key autoincrement, molname text, smiles text, molblock text)")
    cursor.execute("Create Table If Not Exists property_values (molid integer, propid integer, propvalue text)")
    cursor.execute("Create Table If Not Exists atoms (molid integer, iteration integer, atomid text, atom_index integer)")
    #cursor.execute("Create Table If Not Exists smarts (atomid text Unique, smarts Text, symbol Text, in_ring integer, is_aromatic integer, hvy_degree integer, hvy_valence integer, hcount integer, formal_charge integer, atomic_number integer, mass integer)")
    #cursor.execute("Create Table If Not Exists parents (atomid text Unique, parentid text, iteration integer)")
    cursor.execute("Create View If Not Exists atom_types As Select atomid, min(iteration) iteration, /*Count(distinct iteration) Niterations,*/ Count(distinct molid) frequency From atoms Group By atomid")


def main():
    if len(sys.argv) > 4:
        # input sdf
        sdfile = sys.argv[1]
        # output sqlite file
        db = sys.argv[2]
        #maximum depth of neighbor search
        max_depth = int(sys.argv[3])
        # which atom to include atom_properties in output smarts string.
        # 0 for none; 1 for root, 2 for tail, 3 for all, 4 for root and tail
        with_properties = int(sys.argv[4])
        
        con = sqlite3.connect(db)
        cur = con.cursor()
        makeTables(cur)
        processSDF(cur, sdfile, max_depth, with_properties)
        con.commit()
    else:
        if len(sys.argv) > 1:
            smi = sys.argv[1]
        else:
            smi = "c1ccccc1C(=O)NC"
        mol = MolFromSmiles(smi)
        depth = 2
        with_properties = 1
        fragments = fragment_mol(mol, depth, with_properties)
        # pick your poison
        show_results(fragments)
        #show_results(set(fragments))
        show_results(["[$(%s)]" % recursive_smarts for recursive_smarts in set(fragments)])

if __name__ == "__main__":
    main()
