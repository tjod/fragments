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
