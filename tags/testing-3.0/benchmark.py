#This file demonstrates the use of fast fruchterman reingold layout from
#fastlayout extension

import igraph as ig
import numpy as np
import fastlayout as fl
import time

# load graph data
g = ig.load( filename='data/medium_example_network.txt', format='ncol' )

# set various constants 
dimension = 2
max_it = 500
temp = 1
# get number of nodes
nodes = len(g.vs)
# get edges in the right format
edgelist = np.array( [edge.tuple for edge in g.es], np.int32)

# call igraph first
# this is only to compare - timing as well as output
print 'igraph'
start = time.clock()
igLayout = g.layout_fruchterman_reingold()
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fruchterman_reingold) is: ", elapsed
ig.plot(g, "output/plot_igraph.png", layout=igLayout)


# call all versions with old data format
# this is only for reference
print 'layout_fr (float baseline)'
start = time.clock()
pos = fl.layout_fr(nodes*dimension, edgelist, max_it, temp )
pos = np.reshape(pos, (nodes, dimension))
floatLayout = ig.Layout( tuple(map(tuple, pos )) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr) is: ", elapsed
ig.plot(g, "output/plot_fr_float_baseline.png", layout=floatLayout)

print 'layout_fr_omp'
start = time.clock()
pos = fl.layout_fr_omp( nodes*dimension, edgelist, max_it, temp )
pos = np.reshape(pos, (nodes, dimension))
float2Layout = ig.Layout( tuple(map(tuple, pos)) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr_omp) is: ", elapsed
ig.plot(g, "output/plot_fr_omp.png", layout=float2Layout)

print 'layout_fr_simd'
start = time.clock()
pos = fl.layout_fr_simd( nodes*dimension, edgelist, max_it, temp )
pos = np.reshape(pos, (nodes, dimension))
float3Layout = ig.Layout( tuple(map(tuple, pos)) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr_simd) is: ", elapsed
ig.plot(g, "output/plot_fr_simd.png", layout=float3Layout)

print 'layout_fr_omp_simd'
start = time.clock()
pos = fl.layout_fr_omp_simd( nodes*dimension, edgelist, max_it, temp )
pos = np.reshape( pos, (nodes, dimension) )
float6Layout = ig.Layout( tuple(map(tuple, pos)) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr_omp_simd) is: ", elapsed
ig.plot(g, "output/plot_fr_omp_simd.png", layout=float6Layout)

# call new implementations with strided data layout

# NOTE
# this one works different as it uses a different layout 
# for the data, so no reshaping has to be done, but you
# have to zip first and second half of the array
print 'layout_fr_stride (strided float baseline)'
start = time.clock()
pos = fl.layout_fr_stride( nodes*dimension, edgelist, max_it, temp )
float6Layout = ig.Layout( tuple(zip(pos[0:nodes], pos[nodes:2*nodes])) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr_stride) is: ", elapsed
ig.plot(g, "output/plot_fr_stride.png", layout=float6Layout)

print 'layout_fr_stride_simd'
start = time.clock()
pos = fl.layout_fr_stride_simd( nodes*dimension, edgelist, max_it, temp )
myLayout = ig.Layout( tuple(zip(pos[0:nodes], pos[nodes:2*nodes])) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr_stride_simd) is: ", elapsed
ig.plot(g, "output/plot_fr_simd_stride.png", layout=myLayout)

print 'layout_fr_stride_omp'
start = time.clock()
pos = fl.layout_fr_stride_omp( nodes*dimension, edgelist, max_it, temp )
myLayout = ig.Layout( tuple(zip(pos[0:nodes], pos[nodes:2*nodes])) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr_stride_omp) is: ", elapsed
ig.plot(g, "output/plot_fr_omp_stride.png", layout=myLayout)

print 'layout_fr_stride_omp_simd'
start = time.clock()
pos = fl.layout_fr_stride_omp_simd( nodes*dimension, edgelist, max_it, temp )
myLayout = ig.Layout( tuple(zip(pos[0:nodes], pos[nodes:2*nodes])) )
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (layout_fr_stride_omp_simd) is: ", elapsed
ig.plot(g, "output/plot_fr_omp_stride_simd.png", layout=myLayout)
