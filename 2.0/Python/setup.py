#!/usr/bin/env python3

import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

import numpy as np

ext_modules = [
    Extension('pgenlib',
              sources = ['src/pgenlib/pgenlib.pyx', 'src/plink2/pgenlib_ffi_support.cc', 'src/plink2/include/pgenlib_misc.cc', 'src/plink2/include/pgenlib_read.cc', 'src/plink2/include/pgenlib_write.cc', 'src/plink2/include/plink2_base.cc', 'src/plink2/include/plink2_bits.cc'],
              language = "c++",
              # do not compile as c++11, since cython doesn't yet support
              # overload of uint32_t operator
              # extra_compile_args = ["-std=c++11", "-Wno-unused-function"],
              # extra_link_args = ["-std=c++11"],
              extra_compile_args = ["-std=c++98", "-Wno-unused-function", "-Wno-macro-redefined"],
              extra_link_args = ["-std=c++98"],
              include_dirs = [np.get_include()]
              )
    ]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Pgenlib",
    version="0.81.3",
    author="Christopher Chang",
    author_email="chrchang@alumni.caltech.edu",
    description="Python wrapper for pgenlib's basic reader and writer.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chrchang/plink-ng",
    project_urls={
        "Bug Tracker": "https://github.com/chrchang/plink-ng/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.5",
    ext_modules=cythonize(ext_modules)
)
