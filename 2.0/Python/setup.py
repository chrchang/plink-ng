#!/usr/bin/python

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

import numpy as np

ext_modules = [
    Extension('pgenlib',
              sources = ['pgenlib.pyx', '../pgenlib_python_support.cpp', '../pgenlib_internal.cpp', '../plink2_base.cpp'],
              language = "c++",
              # do not compile as c++11, since cython doesn't yet support
              # overload of uint32_t operator
              # extra_compile_args = ["-std=c++11", "-Wno-unused-function"],
              # extra_link_args = ["-std=c++11"],
              extra_compile_args = ["-std=c++98", "-Wno-unused-function"],
              extra_link_args = ["-std=c++98"],
              include_dirs = [np.get_include()]
              )
    ]

setup(name = 'Pgenlib',
      version = '0.7',
      description = "Wrapper for pgenlib's basic reader and writer.",
      ext_modules = cythonize(ext_modules))
