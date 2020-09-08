import sys
import sqlite3
import json

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

mol_db = sys.argv[1]
con = sqlite3.connect(mol_db)
cur = con.cursor()
if len(sys.argv) > 2:
    atomid = sys.argv[2]
    p = parentage(con, atomid)
    print (json.dumps(p, indent=3))
else:
    cur.execute("Select atomid, iteration, frequency From atom_types /* Where iteration > 4 And frequency > 75 */ Order By iteration Desc, frequency Desc")
    for row in cur:
        atomid = row[0]
        iteration = row[1]
        frequency = row[2]
        p = parentage(con, atomid)
        print (iteration, frequency)
        print (json.dumps(p, indent=3))
