import sys
import sqlite3
import math

if len(sys.argv) < 2:
    print ("usage: log_value db3 property_name")
    sys.exit()
db = sys.argv[1]
propname = sys.argv[2]

con = sqlite3.connect(db)
cur = con.cursor()
icur = con.cursor()
cur.execute("Insert Into property_names (propname) Values (?)", ["log_"+propname])
propid = cur.lastrowid
cur.execute("Select molid, propvalue From property_values Join property_names Using(propid) Where propname=?", [propname])
for row in cur:
    molid = row[0]
    val = float(row[1])
    logval = math.log(val)
    icur.execute ("Insert Into property_values (molid, propid, propvalue) Values (?,?,?)", (molid, propid, logval))

con.commit()