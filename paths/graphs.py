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
from rdkit import Chem
import sys
import networkx as nx
from networkx.utils import pairwise

from networkx.algorithms import isomorphism

def makeGraphsFromSupplier(file, supplier):
    glist = []
    for mol in supplier:
        glist.append(mol_to_graph(mol))
    return glist

def mol_to_graph(mol, attr_fn=None, smarts=False):
    # create nx.Graph with atom attributes from optional function; bonds attr are bond orders alone
    g = nx.Graph(name=mol.GetProp("_Name"))
#     print (json.dumps(nx.to_dict_of_dicts(nx.minimum_spanning_tree(g)), indent=2))
#     print (g.nodes.data())
#     print (json.dumps(nx.to_dict_of_dicts(g), indent=2))
#     print (nx.get_node_attributes(g, "attr"))
#     print (nx.info(g))
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if attr_fn and smarts:
            atom_attr = attr_fn(atom)
            g.add_node(idx, smarts=atom_attr)
            for bond in mol.GetBonds():
                g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_order=bond.GetBondTypeAsDouble())
        else:
            nrings = 0
            if atom.IsInRing():
                sssr = Chem.GetSymmSSSR(atom.GetOwningMol())
                atom_index = atom.GetIdx()
                for ring in sssr:
                    if len(ring) < 7 and atom_index in ring: nrings += 1
            atom_attr = {
              "number":    atom.GetAtomicNum() \
            , "charge":    atom.GetFormalCharge() \
            , "implicitH": atom.GetNumImplicitHs() \
            , "degree":    atom.GetDegree() \
            , "aromatic":  atom.GetIsAromatic() \
            , "nrings":    nrings \
            , "chirality": None if not atom.HasProp("_CIPCode") else atom.GetProp("_CIPCode")
             }
            g.add_node(idx, attr=atom_attr)
            for bond in mol.GetBonds():
                g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), attr={"bond_order":bond.GetBondTypeAsDouble()})
    return g
            
def atom_match(a,b):
    # helper function when comparing mol graphs for isomorphism
    if "smarts" in a:
        return a["smarts"] == b["smarts"]
    elif "attr" in a:
        return a["attr"] == b["attr"]

def bond_match(a,b):
    # helper function when comparing mol graphs for isomorphism
    if "bond_order" in a:
        return a["bond_order"] == b["bond_order"]
    elif "attr" in a:
        return a["attr"] == b["attr"]

def compare(g1,g2):
    # compare lists of graphs and return list of isomorphs
    isomorphs = []
    for g in g2:
        num_nodes = nx.number_of_nodes(g)
        num_edges = nx.number_of_edges(g)
        #print (g.graph)
        for gg in g1:
            #print (" ",gg.graph)
            #if nx.is_isomorphic(g, gg):
            if   num_nodes == nx.number_of_nodes(gg) \
             and num_edges == nx.number_of_edges(gg):
                GM = isomorphism.GraphMatcher(g, gg, node_match=atom_match, edge_match=bond_match)
                #if nx.is_isomorphic(g, gg, node_match=atom_match, edge_match=bond_match):
                if GM.is_isomorphic():
                    isomorphs.append((g,gg))
                    #print (g.graph, gg.graph, GM.mapping)
                    #print (">>>", gg.graph, GM.mapping)
    return isomorphs

def is_isomorphic(g1,g2):
    GM = isomorphism.GraphMatcher(g1, g2, node_match=atom_match, edge_match=bond_match)
    return GM.is_isomorphic()

def cycles(glist):
    for g in glist:
        print (g.graph)
#         for bridge in nx.bridges(g):
#             print ([x+1 for x in bridge])
        for cycle in nx.cycle_basis(g):
            print ([x+1 for x in cycle])
#         print ([x+1 for x in nx.minimum_spanning_tree(g)])
#         print ([x for x in nx.chain_decomposition(g)])
#         for chain in nx.chain_decomposition(g):
#             print ([x+1 for x in chain])

# def make_nest(dict_of_lists, root):
#     # make nested lists of tree; fails for graph with cycles
#     #print (dict_of_lists[root])
#     nodes = dict_of_lists[root]
#     nest = []
#     for node in nodes:
#         nest.append(node)
#         nn = make_nest(dict_of_lists, node)
#         if len(nn) > 0: nest.append(nn)
#     return nest
# 
# def nest_to_smiles(g, nest):
#     smiles = []
#     for atom in nest:
#         if type(atom) == type(list()):
#             smiles.append(nest_to_smiles(g, atom))
#         else:
#             smiles.append(g.nodes[atom]["attr"]["symbol"])
#     return smiles
        
def fragment_bfs(g, depth):
    # make breadth first tree and nest of atoms in graph
    fragments = []
    #mol_atom_attr = nx.get_node_attributes(g, "attr")
    for atom in range(0, g.number_of_nodes()):
        #print (atom, g.nodes[atom]["attr"])
        tree = nx.bfs_tree(g, atom, depth_limit=depth)
        subgraph = g.subgraph(tree.nodes()).copy()
        subgraph.graph["atom"] = atom
        fragments.append(subgraph)
    return fragments
             
def main():
    # find isomorphs in sd file #1 of each smiles in file #2]
    sdfile = sys.argv[1]
    supplier = Chem.SDMolSupplier(sdfile, removeHs=False)
    g1 = makeGraphsFromSupplier(sdfile, supplier)
    smifile = sys.argv[2]
    supplier = Chem.SmilesMolSupplier(smifile, titleLine=False)
    g2 = makeGraphsFromSupplier(smifile, supplier)
    
    g = g2[0]
    #print (g.graph)
    # list atoms
    # for i in g.nodes:
    #     print (g.nodes[i]["attr"])
    
    #print (nx.dominating_set(g, 0))
    #find isomorphs
    isomorphs = compare(g1,g2)
    for (a,b) in isomorphs:
        if a.graph != b.graph: print (a.graph, b.graph)
    
    # find rings
    #cycles(g2[:5])
    
    # find fragments rooted at each atom
    depth = 0
    fragments = fragment_bfs(g, depth)
    #for f in fragments:
    for atom in range(0, g.number_of_nodes()):
        f = fragments[atom]
        print (f.graph, f.nodes())
    
    # make nest (smiles-like) from mol - fails for rings :(
    #dl = nx.to_dict_of_lists(g2[1])
    #print(dl)
    #nest = make_nest(dl, 0)
    #print (nest)

if __name__ == "__main__":
    main()
