.. _install:

Installation
============

This document is about how to install ParallelFR. There are two options: either download the wheel file or build it from source.

Getting the wheel file
----------------------

...

Building from source
--------------------

To build this extension from source you will need

* `SWIG <http://www.swig.org/>`_ to make the link from C to python
* hardware and compiler supporting at least SSE2 instruction set

  For Intel chips, SSE2 was introduced with *Pentium 4* in 2001, AMD supports SSE2 since *Opteron*/*Athlon 64* released in 2003.
* a `compiler supporting OpenMP <http://openmp.org/wp/openmp-compilers/>`_, version 2.5. This is i.e. fullfilled in GCC >= 4.2

Then, there is a `setup.py` for building and installing the extension like

.. code-block:: bash
  
  $ python setup.py build
  $ python setup.py install
  
Todo
----

* whole part on wheel files
