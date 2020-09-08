import sqlite3
import networkx as nx
import sys
from graphs import atom_match, bond_match

def is_isomorphic(g1,g2):
    #if nx.could_be_isomorphic(g1, g2):
    #if nx.number_of_nodes(g1) ==  nx.number_of_nodes(g2) and  nx.number_of_edges(g1) ==  nx.number_of_edges(g2):
    #if nx.faster_could_be_isomorphic(g1, g2) and nx.fast_could_be_isomorphic(g1, g2) and nx.could_be_isomorphic(g1, g2):
    return nx.is_isomorphic(g1,g2, node_match=atom_match, edge_match=bond_match)
#     else:
#         return False

db = sys.argv[1]
con = sqlite3.connect(db)
cur = con.cursor()
cur.execute ("Create Table If Not Exists isomorphs(a_graphid Integer References graph(graphid), b_graphid Integer References graph(graphid))")
cur.execute("Create View If Not Exists isosmarts As Select a.atomid||'.'||b.atomid From isomorphs Join graphs a On (a_graphid=a.graphid) Join graphs b On (b_graphid=b.graphid)")
cur2 = con.cursor()
cur.execute("Select cansmiles, Count(atomid) cc From graphs Where iteration > 0 Group By cansmiles Having cc > 1")
for row in cur:
    cansmi = row[0]
    graphs = []
    graphid = []
    cur2.execute("Select graphid, jit_graph From graphs Where atomid Is Not Null and cansmiles=?", [cansmi])
    
    for grow in cur2:
        graphid.append(grow[0])
        graphs.append(nx.jit_graph(grow[1]))
    ngraphs = len(graphs)
    
    if ngraphs > 1:
        nmax = ngraphs * (ngraphs - 1) / 2
        niso = 0
        for i in range(ngraphs):
            g1 = graphs[i]
            for j in range(0,i):
                g2 = graphs[j]
                if is_isomorphic(g1, g2):
                    cur2.execute("Insert Into isomorphs (a_graphid, b_graphid) Values (?,?)", (graphid[i], graphid[j]))
                    cur2.execute("Insert Into isomorphs (a_graphid, b_graphid) Values (?,?)", (graphid[j], graphid[i]))
                    niso += 1
                    #print (i, j)
        print ("%s %d of %d/%d" %(cansmi, niso, ngraphs, nmax))
                
con.commit()