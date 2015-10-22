# This file demonstrates the use of fast fruchterman reingold layout from
# fastlayout extension in graph tool

import graph_tool.all as gt
import numpy as np
import fastlayout as fl
import time

# create a network to have some graph data
g = gt.collection.data["celegansneural"]

# set various constants 
dimension = 2
max_it = 500
temp = 1
# get number of nodes and an edgelist
nodes = g.num_vertices()
get_vids = np.vectorize(lambda v:  g.vertex_index[v], otypes=[np.int32])
edgelist = np.array( [[edge.source(), edge.target()] for edge in g.edges()])
edgelist =  get_vids(edgelist)

# call fastlayout function
start = time.clock()
pos = fl.layout_fr_omp_simd( nodes*dimension, edgelist, max_it, temp )

pos = np.reshape(pos, (nodes, dimension))
pos_prop = g.new_vertex_property('vector<double>')
for v in g.vertices():
    i = g.vertex_index[v]
    pos_prop[v] = pos[i]#[::-1]
    
fastlayout = time.clock()
fastlayout = fastlayout - start
    
# draw the graph in graph_tool
gt.graph_draw(g, pos=pos_prop, output="output/fr-fastlayout.pdf")



# call fastlayout function with new data structure
start = time.clock()
pos = fl.layout_fr_stride_omp_simd( nodes*dimension, edgelist, max_it, temp )

pos = zip(pos[0:nodes], pos[nodes:2*nodes])
pos_prop2 = g.new_vertex_property('vector<double>')
for v in g.vertices():
    i = g.vertex_index[v]
    pos_prop2[v] = pos[i]#[::-1]
    
fastlayout_strided = time.clock()
fastlayout_strided = fastlayout - start
    
# draw the graph in graph_tool
gt.graph_draw(g, pos=pos_prop2, output="output/fr-fastlayout-strided.pdf")



# call and plot the reference implementation of graph-tool
start = time.clock()
reference = gt.fruchterman_reingold_layout(g, n_iter=max_it)
reflayout = time.clock()
reflayout = reflayout - start
gt.graph_draw(g, pos=reference, output="output/fr-reference.pdf")

print "Time for fastlayout is: ", fastlayout
print "Time for reference implementaion in gt is: ", reflayout
print "Speedup is x%.2f" % (reflayout / fastlayout)

# Graph-Tool references
# =====================
# * vertices and edgelist
#   https://graph-tool.skewed.de/static/doc/quickstart.html#iterating-over-all-vertices-or-edges
# * calling a layout function
#   https://graph-tool.skewed.de/static/doc/draw.html#graph_tool.draw.fruchterman_reingold_layout
# * draw graphs
#   https://graph-tool.skewed.de/static/doc/draw.html#graph_tool.draw.graph_draw
