# PyATria
# Language: Python
# Input: CSV (network)
# Output: NOA (central nodes and centrality values)
# Tested with: PluMA 1.1, Python 3.6
# Dependency: numpy==1.16.0, networkx==2.2

PluMA plugin that computes Ablatio Triadum centrality (Cickovski et al, 2015, 2017) in Python.
The plugin accepts as input a signed and weighted network in CSV format, with rows and columns
both corresponding to nodes in the network and position (i, j) holding the weight of the edge
from node i to node j.

The plugin will then produce as output a ranked set of central nodes, in NOde Attribute (NOA)
file format which can be imported into Cytoscape.  Centrality then becomes an attribute of
every node in the network, and is available for incorporation into further downstream analysis and/or
visualization.

The plugin uses a slightly modified version of PythonDS (to obtain it, run 'unzip pythonds.zip')
