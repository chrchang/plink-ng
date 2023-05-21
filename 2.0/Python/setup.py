#!/usr/bin/env python3

import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

import numpy as np

# I'd prefer not to include static zstd here, but I get build problems on macOS
# without it.
clib = ('clib',
        {'sources': ['src/plink2/libdeflate/lib/adler32.c',
                     'src/plink2/libdeflate/lib/crc32.c',
                     'src/plink2/libdeflate/lib/deflate_compress.c',
                     'src/plink2/libdeflate/lib/deflate_decompress.c',
                     'src/plink2/libdeflate/lib/gzip_compress.c',
                     'src/plink2/libdeflate/lib/gzip_decompress.c',
                     'src/plink2/libdeflate/lib/utils.c',
                     'src/plink2/libdeflate/lib/zlib_compress.c',
                     'src/plink2/libdeflate/lib/zlib_decompress.c',
                     'src/plink2/libdeflate/lib/arm/arm_cpu_features.c',
                     'src/plink2/libdeflate/lib/x86/x86_cpu_features.c',
                     'src/plink2/zstd/lib/common/debug.c',
                     'src/plink2/zstd/lib/common/entropy_common.c',
                     'src/plink2/zstd/lib/common/error_private.c',
                     'src/plink2/zstd/lib/common/fse_decompress.c',
                     'src/plink2/zstd/lib/common/pool.c',
                     'src/plink2/zstd/lib/common/threading.c',
                     'src/plink2/zstd/lib/common/xxhash.c',
                     'src/plink2/zstd/lib/common/zstd_common.c',
                     'src/plink2/zstd/lib/compress/fse_compress.c',
                     'src/plink2/zstd/lib/compress/hist.c',
                     'src/plink2/zstd/lib/compress/huf_compress.c',
                     'src/plink2/zstd/lib/compress/zstd_compress.c',
                     'src/plink2/zstd/lib/compress/zstd_compress_literals.c',
                     'src/plink2/zstd/lib/compress/zstd_compress_sequences.c',
                     'src/plink2/zstd/lib/compress/zstd_compress_superblock.c',
                     'src/plink2/zstd/lib/compress/zstd_double_fast.c',
                     'src/plink2/zstd/lib/compress/zstd_fast.c',
                     'src/plink2/zstd/lib/compress/zstd_lazy.c',
                     'src/plink2/zstd/lib/compress/zstd_ldm.c',
                     'src/plink2/zstd/lib/compress/zstd_opt.c',
                     'src/plink2/zstd/lib/compress/zstdmt_compress.c',
                     'src/plink2/zstd/lib/decompress/huf_decompress.c',
                     'src/plink2/zstd/lib/decompress/zstd_ddict.c',
                     'src/plink2/zstd/lib/decompress/zstd_decompress.c',
                     'src/plink2/zstd/lib/decompress/zstd_decompress_block.c'],
         'macros': [('STATIC_ZSTD', None), ('ZSTD_DISABLE_ASM', None), ('LIBDEFLATE_STATIC', None)],
         'include_dirs': ['src/plink2/zstd/lib', 'src/plink2/zstd/lib/common', 'src/plink2/libdeflate']
         })

ext_modules = [
    Extension('pgenlib',
              sources = ['src/pgenlib/pgenlib.pyx',
                         'src/plink2/pgenlib_ffi_support.cc',
                         'src/plink2/pvar_ffi_support.cc',
                         'src/plink2/include/pgenlib_misc.cc',
                         'src/plink2/include/pgenlib_read.cc',
                         'src/plink2/include/pgenlib_write.cc',
                         'src/plink2/include/plink2_base.cc',
                         'src/plink2/include/plink2_bgzf.cc',
                         'src/plink2/include/plink2_bits.cc',
                         'src/plink2/include/plink2_string.cc',
                         'src/plink2/include/plink2_text.cc',
                         'src/plink2/include/plink2_thread.cc',
                         'src/plink2/include/plink2_zstfile.cc'],
              language = "c++",
              # do not compile as c++11, since cython doesn't yet support
              # overload of uint32_t operator
              # extra_compile_args = ["-std=c++11", "-Wno-unused-function"],
              # extra_link_args = ["-std=c++11"],
              extra_compile_args = ["-std=c++98", "-Wno-unused-function", "-Wno-macro-redefined", "-DSTATIC_ZSTD", "-DZSTD_DISABLE_ASM", "-DLIBDEFLATE_STATIC"],
              extra_link_args = ["-std=c++98"],
              include_dirs = [np.get_include()]
              )
    ]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Pgenlib",
    version="0.90.0",
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
    libraries=[clib],
    ext_modules=cythonize(ext_modules)
)
