# Leicester

Data files for Leicester network. The network models a region ("RA") of the Leicester road network as it was in 1997. The region is essential bounded by London Road (A6), Waterloo Way (since renamed Tigers Way), Lancaster Road and Granville Road.

nodes.txt - node ID, longitude, latitute
links.txt - Link ID, from node, to node, Cost0 (free flow travel time in seconds), capacity (vehicles/15 minutes)
counts.txt - 15 minute traffic counts on a subset of network links between 19 and 23 May 1997
monitored-links.txt - ID of links equipped with vehicle counters 
priormean.txt - Origin, Destination, Prior mean of OD traffic volumes (based on dodgy earlier survey!)
Afull.txt - link-path incidence matrix for all network links (built using buildA function in library transportation)
A.txt - link-path incidence matrix for only monitored links (built using buildA function in library transportation)
O.txt - origin node corresponding to each path (column) of A (and Afull)
D.txt - destination node corresponding to each path (column) of A (and Afull)





