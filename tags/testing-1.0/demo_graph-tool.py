#This file demonstrates the use of fast fruchterman reingold layout from
#fastlayout extension in graph tool

import graphtool as gt
import numpy as np
import fastlayout as fl
import math
import time

# create a price network to have some graph data
g = gt.price_network(300)

# set various constants 
dimension = 2
max_it = 500
temp = 1
# get number of nodes and an edgelist
nodes = len(g.vertices())
edgelist = np.array( [edge.tuple for edge in g.edges()], np.int32)

start = time.clock()
# place nodes initially uniformly at random in 
# [-sqrt(nodes)/2, +sqrt(nodes)/2]^2
pos = math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
pos = pos.astype( np.float32) )

# call fastlayout function
fl.layout_fr_omp_simd( edgelist, pos, max_it, temp )
fastlayout = time.clock()
fastlayout = fastlayout - start

# draw the graph in graph_tool
gt.graph_draw(g, pos=pos, output="output/price_network-fr-fastlayout.pdf")

# call and plot the reference implementation of graph-tool
start = time.clock()
reference = gt.fruchterman_reingold_layout(g, n_iter=1000)
reflayout = time.clock()
reflayout = reflayout - start
gt.graph_draw(g, pos=reference, output="output/price_network-fr-reference.pdf")

print "Time for fastlayout is: ", fastlayout
print "Time for reference implementaion in gt is: ", reflayout

# Graph-Tool references
# =====================
# * vertices and edgelist
#   https://graph-tool.skewed.de/static/doc/quickstart.html#iterating-over-all-vertices-or-edges
# * calling a layout function
#   https://graph-tool.skewed.de/static/doc/draw.html#graph_tool.draw.fruchterman_reingold_layout
# * draw graphs
#   https://graph-tool.skewed.de/static/doc/draw.html#graph_tool.draw.graph_draw
