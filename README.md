# Leicester

Data files for Leicester network. The network models a region ("RA") of the Leicester road network as it was in 1997. The region is essentially bounded by London Road (A6), Waterloo Way (since renamed Tigers Way), Lancaster Road and Granville Road.

nodes.txt - node ID, longitude, latitute

links.txt - Link ID, from node, to node, Cost0 (free flow travel time in seconds), capacity (vehicles/15 minutes)

counts.txt - 15 minute traffic counts on a subset of network links between 4 and 8 May 1998

monitored-links.txt - ID of links equipped with vehicle counters 

priormean.txt - Origin, Destination, Prior mean of OD traffic volumes (based on dodgy earlier survey!)

Afull.txt - link-path incidence matrix for all network links (built using buildA function in library transportation)

A.txt - link-path incidence matrix for only monitored links (built using buildA function in library transportation)

O.txt - origin node corresponding to each path (column) of A (and Afull)

D.txt - destination node corresponding to each path (column) of A (and Afull)

Files with s and l extensions (e.g. Al.txt, Ds.txt) are variations on the above with smaller and large versions of the A matrix respectively (corresponding to less or more routes).





