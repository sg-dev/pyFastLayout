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
        version     = "0.3.0",
        py_modules  = ['fastlayout'],
        ext_modules = [_fastlayout],
        
        # Choose your license
        license='MIT',
        
        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Information Analysis',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: MIT License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 2.7',
        ],
        )
