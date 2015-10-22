#This file demonstrates the use of fast fruchterman reingold layout from
#fastlayout extension

import igraph as ig
import numpy as np
import fastlayout as fl
import math
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
nodes = len(g.vs)
edges = np.array( [edge.tuple for edge in g.es], np.int32)
# place nodes initially uniformly at random in 
# [-sqrt(nodes)/2, +sqrt(nodes)/2]^2
pos = math.sqrt(nodes) * np.random.random_sample( (nodes, dimension) ) - math.sqrt(nodes)/2
floatPos = np.copy( pos.astype(np.float32) )
float2Pos = np.copy( pos.astype(np.float32) )
float3Pos = np.copy( pos.astype(np.float32) )
float4Pos = np.copy( pos.astype(np.float32) )
float5Pos = np.copy( pos.astype(np.float32) )
float6Pos = np.copy( pos.astype(np.float32) )

#print np.array_str( floatPos )
#print np.array_str( float2Pos )

print 'igraph'
start = time.clock()
igLayout = g.layout_fruchterman_reingold()
elapsed = time.clock()
elapsed = elapsed - start
print "Time spent in (function name) is: ", elapsed
ig.plot(g, "output/plot_igraph.png", layout=igLayout)

print 'float baseline:'
fl.layout_fr( edges, floatPos, max_it, temp )
floatLayout = ig.Layout( tuple(map(tuple, floatPos )) )
ig.plot(g, "output/plot_float_baseline.png", layout=floatLayout)

#print 'layout_fr_omp'
#fl.layout_fr_omp( edges, float2Pos, max_it, temp )
#float2Layout = ig.Layout( tuple(map(tuple, float2Pos)) )
#ig.plot(g, "output/plot_fr_omp.png", layout=float2Layout)

#print 'layout_fr_simd'
#fl.layout_fr_simd( edges, float3Pos, max_it, temp )
#float3Layout = ig.Layout( tuple(map(tuple, float3Pos)) )
#ig.plot(g, "output/plot_fr_simd.png", layout=float3Layout)

#print 'layout_fr_simd_1x4'
#fl.layout_fr_simd_1x4( edges, float4Pos, max_it, temp )
#float4Layout = ig.Layout( tuple(map(tuple, float4Pos)) )
#ig.plot(g, "output/plot_fr_simd_1x4.png", layout=float4Layout)

#print 'layout_fr_2x2'
#fl.layout_fr_2x2( edges, float5Pos, max_it, temp )
#float5Layout = ig.Layout( tuple(map(tuple, float5Pos)) )
#ig.plot(g, "output/plot_fr_2x2.png", layout=float5Layout)

print 'layout_fr_omp_simd'
fl.layout_fr_omp_simd( edges, float6Pos, max_it, temp )
float6Layout = ig.Layout( tuple(map(tuple, float6Pos)) )
ig.plot(g, "output/plot_fr_omp_simd.png", layout=float6Layout)

#print np.array_str( floatPos )
#print np.array_str( float2Pos )
