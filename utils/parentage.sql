With Recursive parent_of(atomid, parentid, iteration) As (
 Select atomid, parentid, iteration From parents Where atomid = -707346364
  Union
 Select parents.atomid, parents.parentid, parents.iteration From parents, parent_of Where parents.atomid = parent_of.parentid
)
Select atomid, parentid, iteration, smarts From parent_of Left Join smarts Using (atomid) Order By iteration Desc;
