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
'''
Created on May 4, 2020

@author: tj o'donnell
'''
#from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
import json
import sqlite3
from graphs import mol_to_graph, compare, fragment_bfs
import networkx as nx

bond_symbol = ["","-","=","#","~",":"]

# return list of select atom properties
def properties(atom):
    nrings = 0
    if atom.IsInRing():
        sssr = Chem.GetSymmSSSR(atom.GetOwningMol())
        atom_index = atom.GetIdx()
        for ring in sssr:
            if len(ring) < 7 and atom_index in ring: nrings += 1
    return (atom.GetAtomicNum(), atom.GetFormalCharge(),
            atom.GetNumImplicitHs() + atom.GetNumExplicitHs(), \
             #atom.GetNumImplicitHs(), \
            atom.GetDegree(), \
            'a' if atom.GetIsAromatic() else 'A', \
            nrings \
            )
    
# construct atom symbol, optionally using smarts with atom properties
def atom_type(atom, as_smarts=True):
    
    if as_smarts:
        if atom.HasProp("smarts"):
            return atom.GetProp("smarts")
        atom_string = "[#%d;%+d;H%d;D%d;%s;R%d]" % properties(atom)
        atom.SetProp("smarts", atom_string)
    else:
        atom_symbol = atom.GetSymbol()
        atom_string = atom_symbol.lower() if atom.GetIsAromatic() else atom_symbol
    
    return atom_string

# get neighbors of atom
def neighbors(atom):
    #return ob.OBAtomAtomIter(atom)
    #return [n for n in ob.OBAtomAtomIter(atom)]
    return atom.GetNeighbors()
    #will neighbors visited in sorted order make canonical smarts? NO! It causes different orders each run
    #return sorted([n for n in atom.GetNeighbors()], key=lambda a: hash(properties(a)))

# create breadth-first graph of neighbors of atom up to depth
def branched_paths(atom, depth, root_atom, visited):
    paths = [{"atom": atom, "depth": depth}]
    visited.add(atom.GetIdx())
    if depth > 0:
        for nbr in neighbors(atom):
            #if nbr.IsHydrogen(): continue
            # don't branch back to root
            #if nbr == root_atom: continue
            # don't branch back to visited atoms; causes excess recursion, e.g. in rings
            if nbr.GetIdx() in visited: continue
            visited.add(nbr.GetIdx())
            nn = branched_paths(nbr, depth-1, atom, visited)[0]
            paths.append(nn)
    return (paths, visited)

# turn neighbour list (path) into smarts string
def path_to_string(atoms, with_properties, root_atom):
    s = ""
    atom0 = atoms[0]["atom"]
    atom_depth = atoms[0]["depth"]

    # lots of output choices
    if with_properties == "all":
        # all atoms
        as_smarts = True
    elif with_properties == "tail":
        # tail atoms
        as_smarts = False if atom_depth > 0 else True
    elif with_properties == "root":
        # if root_atom is None, then this is Ur-root_atom, primary branch point atom
        as_smarts = False if root_atom else True
    elif with_properties == "root+tail":
        # if root_atom is None, then this is Ur-root_atom, primary branch point atom
        as_smarts = False if (root_atom and atom_depth > 0) else True
    elif with_properties == "none":
        # no atoms
        as_smarts = False

    if root_atom:
        # has a root atom, so include bond symbol
        mol = root_atom.GetOwningMol()
        bond = mol.GetBondBetweenAtoms(root_atom.GetIdx(), atom0.GetIdx())
        bond_order = 5 if bond.GetIsAromatic() else int(bond.GetBondTypeAsDouble())
        #bond_order = 5 if root_atom.GetBond(atom0).IsAromatic() else root_atom.GetBond(atom0).GetBondOrder()
    else:
        # Ur-root_atom, primary branch point atom; so no bond symbol
        bond_order = 0
    
    s += "%s%s" % (bond_symbol[bond_order], atom_type(atom0, as_smarts) )    
    for atom in atoms[1:]:
        s += "(%s)" % path_to_string(atom, with_properties, atom0)
    return s

# allow json encoding of Mol, Atom
class MolAtomEncoder(json.JSONEncoder):
    def default(self, z):
        if isinstance(z, Chem.Mol):
            try:
                molname = z.GetProp("_Name")
            except Exception as e:
                molname = str(e)
            return molname
        if isinstance(z, Chem.Atom):
            return z.GetIdx()+1
        else:
            return super().default(z)

def show_results(results):
    # pick your output format(s)
    print ()
    print ("\n".join(["%d: %s" % (i,results[i]) for i in range(0,len(results))]))  # for human consumption?
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

def addMol(cursor, mol, molblock, nmol, verbosity):
    if mol:
        smiles = Chem.MolToSmiles(mol)
        try:
            molname = mol.GetProp("_Name")
        except  Exception as e:
            molname = str(e)
            mol.SetProp("_Name", molname)
            print ("Warning processing name of molecule #%d: %s" % (nmol, molname))
    else:
        molname = "Error"
        smiles = None
    if verbosity > 1: print ("%d: %s" % (nmol, smiles))
    cursor.execute("Insert Into molecule (molname, smiles, molblock) Values (?,?,?)", [molname, smiles, molblock])
    return cursor.lastrowid

def fragment_mol(mol, depth, with_properties, mol_graph):
    fragments = []
    graphs = []
    visits = []
    for atom in mol.GetAtoms():
        #print (atom.GetIndex() + 1)
        #if atom.IsHydrogen(): continue # and depth > 0: continue
        (path, visited) = branched_paths(atom, depth, atom, set())
        visits.append(visited)
        if mol_graph:
            subgraph = mol_graph.subgraph(visited).copy()
            subgraph.graph["atom"] = atom.GetIdx()
            graphs.append(subgraph)
        #print (json.dumps(path, cls=MolAtomEncoder))
        #print ([v+1 for v in visited])
        s = path_to_string(path, with_properties, None)
        fragments.append(s)
    return (fragments, visits, graphs)

def addProp(cursor, imol, prop, val):
    cursor.execute("Select propid From property_names Where propname = (?)", [prop])
    row = cursor.fetchone()
    if row:
        propid = row[0]
    else:
        cursor.execute("Insert Into property_names (propname) Values (?)", [prop])
        propid = cursor.lastrowid
    cursor.execute("Insert Into property_values (molid, propid, propvalue) Values (?,?,?)", [imol, propid, val])
    
def processSDF(cur, sdfile, max_depth, with_properties, limit, verbosity, store_graphs, removeH):
    nmol = 0
    suppl = Chem.SDMolSupplier(sdfile, removeHs=removeH)
    for mol in suppl:
        if limit and nmol >= limit: break
        #print (json.dumps(mol, cls=MolAtomEncoder))
        (molblock, sep, moldata) = suppl.GetItemText(nmol).partition('M  END\n')
        nmol += 1
        molid = addMol(cur, mol, molblock+sep, nmol, verbosity)
        for p in mol.GetPropNames():
            addProp(cur, molid, p, mol.GetProp(p))
        # MW provides a nice test of prediction at depth = 0
        if not mol.HasProp("MW"): addProp(cur, molid, "MW", Descriptors.MolWt(mol))

        previous_fragments = None
        mol_graph = mol_to_graph(mol, attr_fn=atom_type, smarts=True) if store_graphs else None
        #print (mol_graph.graph, mol_graph.nodes())
        for depth in range(0, max_depth+1):
            (fragments, visited, graphs) = fragment_mol(mol, depth, with_properties, mol_graph)

            if verbosity > 1:
                if store_graphs:
                    bfs_fragments = fragment_bfs(mol_graph, depth)
        #             for f in bfs_fragments:
        #                 print (f.graph, f.nodes())
                    print ("compare graphs vs nx bfs fragemented graphs")
                    print (mol_graph.graph, depth)
                    for (a,b) in compare (graphs, bfs_fragments):
                        if a.graph["atom"] != b.graph["atom"]: print (a.graph["atom"], b.graph["atom"])

            for i in range(len(fragments)):
                if previous_fragments and fragments[i] == previous_fragments[i]:
                    # atom has no neighbors past previous depth, e.g. water at depth > 0, CO at depth > 1, et. al.
                    ff = None
                    jit_graph = None
                    cansmi = None
                    nodes = None
                else:
                    ff = fragments[i]
                    cansmi = Chem.MolToSmiles(Chem.MolFromSmarts(ff))
                    nodes = str(visited[i])
                    if store_graphs:
                        nodes = str(graphs[i].nodes())
                        #cansmi = Chem.MolToSmiles(Chem.MolFromSmarts(ff))
                        jit_graph = nx.jit_data(graphs[i])
                cur.execute ("Insert Into atoms (molid, iteration, atomid, atom_index, cansmiles, nodes) Values (?,?,?,?,?,?)", (molid, depth, ff, i+1, cansmi, nodes))
                if store_graphs:
                    cur.execute ("Insert or Ignore Into graphs (iteration, atomid, jit_graph, cansmiles) Values (?,?,?,?)",
                              (depth, ff, jit_graph, cansmi))
            previous_fragments = fragments
    return nmol

def makeTables(cursor):
    # make tables
    cursor.execute("PRAGMA foreign_keys = ON")
    cursor.execute("Create Table If Not Exists property_names (propid Integer Primary Key Autoincrement, propname text Unique)")
    cursor.execute("Create Table If Not Exists molecule (molid Integer Primary Key Autoincrement, molname Text, smiles Text, molblock Text)")
    cursor.execute("Create Table If Not Exists property_values (molid Integer References molecule(molid) On Update Cascade On Delete Cascade DEFERRABLE INITIALLY DEFERRED, propid Integer References property_names(propid) On Update Cascade On Delete Cascade DEFERRABLE INITIALLY DEFERRED, propvalue Text)")
    cursor.execute("Create Table If Not Exists atoms (molid Integer References molecule(molid) On Update Cascade On Delete Cascade DEFERRABLE INITIALLY DEFERRED, iteration Integer, atomid Text, atom_index Integer, cansmiles Text, nodes Text)")
    cursor.execute("Create Table If Not Exists graphs (graphid Integer Primary Key, iteration Integer, atomid Text Unique, jit_graph Text, cansmiles Text)")
    #cursor.execute("Create Table If Not Exists smarts (atomid text Unique, smarts Text, symbol Text, in_ring integer, is_aromatic integer, hvy_degree integer, hvy_valence integer, hcount integer, formal_charge integer, atomic_number integer, mass integer)")
    #cursor.execute("Create Table If Not Exists parents (atomid text Unique, parentid text, iteration integer)")
    cursor.execute("Create View If Not Exists atom_types As Select atomid, cansmiles, Min(iteration) iteration, Count(distinct molid) frequency From atoms Group By atomid")

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="create fragments of input molecule(s)",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--verbosity", type=int, help="increase output verbosity", default=0)
    parser.add_argument("sdf", help="input sd file")
    parser.add_argument("-db", help="output sqlite3 file", default=":memory:")
    parser.add_argument("-d", "--depth", type=int, help="maximum depth to search", default=1)
    parser.add_argument("-k", "--keepH", help="keep explicit H atoms in sd file", action="store_true")
    parser.add_argument("-l", "--limit", type=int, help="limit number of input molecules processed")
    parser.add_argument("-t", "--test", help="one smiles string as a test; no db output")
    parser.add_argument("-p", "--properties", help="include smarts properties on some atoms; root, tail, root+tail, all, none", default="all")
    parser.add_argument("-g", "--graphs", help="test and store networkx graphs", action="store_true")
    return parser

def main():
    parser = parse_args()
    parsed = parser.parse_args()

    # input sdf
    sdfile = parsed.sdf
    # remove H
    removeH = not parsed.keepH
    # limit number of input molecules
    limit = parsed.limit
    # output sqlite file
    db = parsed.db
    #maximum depth of neighbor search
    max_depth = parsed.depth
    # which atom to include atom_properties in output smarts string.
    # n for none; r for root, t for tail, a for all, rt for root and tail
    with_properties = parsed.properties
    # verbose
    verbosity = parsed.verbosity
    # store graphs
    store_graphs = parsed.graphs
    # single smiles to process
    smi = parsed.test
    if smi and sdfile:
        print ("using single smiles in place of sdfile")

    if sdfile:
        con = sqlite3.connect(db)
        cur = con.cursor()
        makeTables(cur)
        nmols = processSDF(cur, sdfile, max_depth, with_properties, limit, verbosity, store_graphs, removeH)
        if verbosity > 0: print ("%d molecules processed" % nmols)
        con.commit()

    elif smi:
        # test single smiles, prodide lots of output
        print (smi)
        mol = Chem.MolFromSmiles(smi)
        mol.SetProp("_Name", "atest")
        mol_graph = mol_to_graph(mol)
        #print (nx.jit_data(mol_graph, 2))
        (fragments, visited, graphs) = fragment_mol(mol, max_depth, with_properties, mol_graph)
        
        # pick your poison
        show_results(fragments)
        if verbosity > 2:
            show_results(list(set(fragments)))
            show_results(["[$(%s)]" % recursive_smarts for recursive_smarts in set(fragments)])

        if verbosity > 1:
            # some tests
            print ("graphs")
            for g in graphs:
                print (max_depth, g.graph, g.nodes())
                print (nx.jit_data(g))
#                 if not nx.is_tree(g):
#                     print (max_depth, g.graph, g.nodes(), nx.cycle_basis(g))
            print ("nx bfs fragmented graphs")
            bfs_fragments = fragment_bfs(mol_graph, max_depth)
            for f in bfs_fragments:
                print (f.graph, f.nodes())
            print ("compare graphs vs nx bfs fragemented graphs")
            for (a,b) in compare (graphs, bfs_fragments):
                #print (a.graph["atom"], nx.is_tree(a), b.graph["atom"], nx.is_tree(b))
                print (a.graph["atom"], b.graph["atom"])

        if verbosity > 0:
            print ()
            print ("identical framents")
            for f in set(fragments):
                same = [i for i, x in enumerate(fragments) if x == f]
                if len(same) > 1: print  (same)
            if max_depth > 0:
                # doesn't handle depth 0, which has no graph edges
                print ("isomorphs")
                isomorphs = compare(graphs, graphs)
                for (a,b) in isomorphs:
                    atoma = a.graph["atom"]
                    atomb = b.graph["atom"]
                    if (atoma != atomb):
                        print (atoma, atomb,
                               a.nodes[atoma]["attr"], 
                               b.nodes[atomb]["attr"],
                               fragments[atoma], fragments[atomb])

    else:
        parser.print_help()

if __name__ == "__main__":
    #import cProfile
    #cProfile.run('main()', 'profile.out')
    main()
