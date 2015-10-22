#This file demonstrates the use of fast fruchterman reingold layout from
#fastlayout extension

import igraph as ig
import numpy as np
import fastlayout as fl
import time

# load graph data
# 
# there are two graphs with this example
# data/test.ncol is a dummy graph with 5 nodes
# data/medium_example_network.txt is a graph with ~5000 nodes
# 
g = ig.load( filename='data/medium_example_network.txt', format='ncol' )
#g = ig.load( filename='data/test.ncol', format = 'ncol' )

# set various constants 
dimension = 2
max_it = 500
temp = 1

# get number of nodes
nodes = len(g.vs)
# get edges in the right format
edgelist = np.array( [edge.tuple for edge in g.es], np.int32)

# calculate positions using fastlayout
pos = fl.layout_fr( nodes*dimension, edgelist, max_it, temp )
# convert positions to igraph Layout object
myLayout = ig.Layout( tuple(zip(pos[0:nodes], pos[nodes:2*nodes])) )

# plot to file
ig.plot(g, "output/plot_fr_omp_stride_simd.png", layout=myLayout)
