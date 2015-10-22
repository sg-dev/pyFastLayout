.. _usage:

Usage
=====

This document explains how to use Parallel FR in combination with either `the python version of igraph <http://igraph.org/python/>`_ or `graph-tool <https://graph-tool.skewed.de/>`_.

Using Parallel FR with igraph
-----------------------------

General Remarks
^^^^^^^^^^^^^^^

Apart from `igraph <http://igraph.org/python/>`_, Parallel FR depends only on the `numpy <http://www.numpy.org/>`_ python package to efficiently transfer arrays (or data in a broader sense) from python to C and back.

To use Parallel FR within your `igraph <http://igraph.org/python/>`_ project you will only have to construct an edgelist and pass it with some parameters to :py:func:`layout_fr`. The resulting positions are then feed to igraph's Layout constructor. In detail it can be broken down to six steps:

1. load the module, `igraph <http://igraph.org/python/>`_ and `numpy <http://www.numpy.org/>`_, which is needed for some array manipulation afterwards

.. code-block:: python
  :linenos:

  import fastlayout as fl
  import igraph
  import numpy

2. load some graph data to be displayed

.. code-block:: python
  :linenos:
  
  g = igraph.load( ... )
  
3. configure the layout function and extract necessary information of the `igraph <http://igraph.org/python/>`_ graph object

.. code-block:: python 
  :linenos:
  
  # define some constants
  dimension = 2
  number_of_iterations = 500
  temperature = 1
  # get number of nodes
  number_of_nodes = len(g.vs)
  # get edges in the right format
  edgelist = np.array( [edge.tuple for edge in g.es], np.int32)

4. call the layout function

.. code-block:: python
  :linenos:
  
  positions = fl.layout_fr(number_of_nodes*dimension, edgelist, number_of_iterations, temperature )

5. reshape the returned positions and generate an `igraph <http://igraph.org/python/>`_ layout object

.. code-block:: python
  :linenos:
  
  myLayout = ig.Layout( tuple(zip(pos[0:nodes], pos[nodes:2*nodes])) )

6. call `igraph <http://igraph.org/python/>`_'s plot function with the computed layout

.. code-block:: python
  :linenos:
  
  igraph.plot(g, "output.png", layout=myLayout)

Example code
^^^^^^^^^^^^
The whole example file for `igraph <http://igraph.org/python/>`_ can be found in :doc:`fastlayout/demo_igraph.py`. Here is its content

.. literalinclude:: ../fastlayout/demo_igraph.py
  :language: python
  :emphasize-lines: 6, 37
  :linenos:

Using Parallel FR with graph-tool
---------------------------------

General Remarks
^^^^^^^^^^^^^^^

Apart from `graph-tool <https://graph-tool.skewed.de/>`_, Parallel FR depends only on the `numpy <http://www.numpy.org/>`_ python package to efficiently transfer arrays (or data in a broader sense) from python to C and back.

To use Parallel FR within your `graph-tool <https://graph-tool.skewed.de/>`_ project you will only have to construct an edgelist and pass it with some parameters to :py:func:`layout_fr`. The resulting positions are then feed as vertex property to `graph-tool <https://graph-tool.skewed.de/>`_'s drawing routine. In detail it can be broken down to six steps:

1. load the necessary modules
2. create some graph data
3. set various constants, get the number of nodes and extract an edgelist
4. call the fast layout function
5. convert output positions into vertex properties
6. draw the graph using the created vertex properties

Example code
^^^^^^^^^^^^

The full example file for `graph-tool <https://graph-tool.skewed.de/>`_ is located in :doc:`fastlayout/graph_tool.py`. Its content is below:

.. literalinclude:: ../fastlayout/demo_graphtool.py
  :language: python
  :emphasize-lines: 6,24
  :linenos: