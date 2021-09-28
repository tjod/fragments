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
#SOFTWARE.import sqlite3
import networkx as nx
import sys
import sqlite3
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
cur.execute ("Drop Table If Exists isomorphs")
cur.execute ("Create Table isomorphs(a_graphid Integer References graph(graphid), b_graphid Integer References graph(graphid))")
#cur.execute("Create View If Not Exists isosmarts As Select a.atomid||'.'||b.atomid From isomorphs Join graphs a On (a_graphid=a.graphid) Join graphs b On (b_graphid=b.graphid)")
cur.execute("Create View If Not Exists first_isomorph As With gid As (Select graphid, atomid,cansmiles From graphs), amin As (Select a_graphid, a.cansmiles, Min(a.atomid, Min(b.atomid)) atomid From isomorphs Join gid a On (a_graphid=a.graphid) Join gid b On (b_graphid=b.graphid) where a_graphid < b_graphid Group By a_graphid) Select a_graphid As graphid, cansmiles, atomid From amin")
cur.execute("Create View If Not Exists usmarts As With uall As (Select graphid,cansmiles,atomid From graphs Left Join isomorphs On graphid=b_graphid where a_graphid Is Null Union Select * from first_isomorph) Select atomid,graphid,a_graphid,frequency From uall Join atom_types Using(atomid) Left Join isomorphs On graphid=a_graphid")
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
