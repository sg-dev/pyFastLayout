#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# fastlayout extension module
_fastlayout = Extension("_fastlayout",
                   ["fastlayout.i","fastlayout.c"],
                   include_dirs = [numpy_include],
                   extra_compile_args=["-O3", "-fopenmp", "-msse3"],
                   extra_link_args=["-fopenmp"],
                   )

# fastlayout setup
setup(  name        = "fastlayout",
        description = "fastlayout implements the Fruchterman-Reingold layout algorithm in a fast, parallel way",
        author      = "Roman Cattaneo",
        author_email = "romanc@ethz.ch",
        url         = "https://www.sg.ethz.ch/",
        version     = "1.0",
        ext_modules = [_fastlayout]
        )
