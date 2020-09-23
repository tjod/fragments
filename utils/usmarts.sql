With
 gid  As (Select graphid, atomid,cansmiles From graphs),
 amin As (Select a_graphid, a.cansmiles, Min(a.atomid, Min(b.atomid)) atomid From isomorphs
      	Join gid a On (a_graphid=a.graphid)
       	Join gid b On (b_graphid=b.graphid)
       	Group By a_graphid)
Select a_graphid As graphid, cansmiles, atomid From amin;
