.. _implementations:

Documentation of all implementations
====================================

To speedup the computation of the Fruchterman Reingold layout two strategies have been adopted: On one hand the performance on one single core was enhanced using vectorization (also known as Intel Intrinsics), on the other hand basic OpenMP pragmas help to use all available computing nodes.

There are two basic categories: Implementations using a strided data layout and implementations using an interleaved data layout. In general strided data layout is expected to be faster because of better locality.

Implementaions with strided data layout
---------------------------------------

There are four versions with this layout

* :py:func:`layout_fr`: version using OpenMP parallelism and vector instructions in addition. This is supposed to be the fastest version.
* :py:func:`layout_fr_simd_only`: vectorized version
* :py:func:`layout_fr_omp_only`: version using OpenMP to use parallelism
* :py:func:`layout_fr_baseline`: baseline implementation

Passing positions to igraph
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
  :linenos:
  
  myLayout = ig.Layout( tuple(zip(positions[0:number_of_nodes], positions[number_of_nodes:2*number_of_nodes])) )

Passing positions to graph-tool
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
  :linenos:
  
  positions = zip(positions[0:number_of_nodes], positions[number_of_nodes:2*number_of_nodes])
  pos_prop2 = g.new_vertex_property('vector<double>')
  for v in g.vertices():
      i = g.vertex_index[v]
      pos_prop2[v] = positions[i]

Implementations with interleaved data layout
--------------------------------------------

Also here, four versions are around, adopting a similar naming convention:

* :py:func:`layout_fr_interleaved`: version using OpenMP parallelism and vector instructions in addition. This is supposed to be the fastest version.
* :py:func:`layout_fr_interleaved_simd_only`: vectorized version
* :py:func:`layout_fr_interleaved_omp_only`: version using OpenMP to make use of parallelism
* :py:func:`layout_fr_interleaved_baseline`: baseline implementation

All but two versions use single precision floating point as calculation type. The functions

* :py:func:`layout_fr_double_baseline` 
* :py:func:`layout_fr_double_simd_only`

are only around for comparison.

Passing positions to igraph
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
  :linenos:
  
  positions = np.reshape(positions, (number_of_nodes, dimension))
  myLayout = ig.Layout( tuple(map(tuple, positions )) )

Passing positions to graph-tool
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
  :linenos:
  
  positions = np.reshape(positions, (number_of_nodes, dimension))
  pos_prop = g.new_vertex_property('vector<double>')
  for v in g.vertices():
      i = g.vertex_index[v]
      pos_prop[v] = positions[i]


List of all implementations
---------------------------

Strided data layout
^^^^^^^^^^^^^^^^^^^

.. function:: layout_fr(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the strided data layout. This function combines vectorization and OpenMP parallelism and is therefore supposed to be the fastest.

.. function:: layout_fr_simd_only(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the strided data layout. This function serves as test and uses vectorization only.
  
  Use this function, if your hardware has no support for parallelism through OpenMP.
  
.. function:: layout_fr_omp_only(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the strided data layout. This function serves as test and uses parallelism through OpenMP only.
  
  Use this function, if your hardware does not support vector instructions up to SSE3.
  
.. function:: layout_fr_baseline(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the strided data layout and is therefore supposed to be faster than :py:func:`layout_fr_interleaved` using the interleaved data layout because of better locality. This function serves as baseline implementation for comparison.
  
  Use this function only, if your hardware neither supports vector instructions nor parallelism through OpenMP.
  
Interleaved data layout
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: layout_fr_interleaved(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the interleaved data layout. This function combines vectorization and OpenMP parallelism and is therefore supposed to be the fastest among the implementations using the interleaved data layout.

.. function:: layout_fr_interleaved_simd_only(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the interleaved data layout. This function serves as test and uses vectorization only.
  
.. function:: layout_fr_interleaved_omp_only(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the interleaved data layout. This function serves as test and uses parallelism through OpenMP only.
  
.. function:: layout_fr_interleaved_baseline(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the interleaved data layout. This function serves as baseline implementation for comparison.

Reference implementations using double precision
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: layout_fr_double_baseline(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the interleaved data layout. This function serves only reference and baseline implementation for comparison.

.. function:: layout_fr_double_simd_only(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

  Returns positions according to Fruchterman Reingold layout algorithm. Returned positions use the interleaved data layout. This function serves as test case for vectorization.
  
