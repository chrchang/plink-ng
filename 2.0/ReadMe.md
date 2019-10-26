This directory contains two major libraries, as well as the PLINK 2.0
application built on top of them.  These are carefully written to be valid C99
(from gcc and clang's perspective, anyway) to simplify FFI development, while
still taking advantage of quite a few C++-specific affordances to improve
safety and occasionally performance.

The first library is plink2_text.  This is a text file reader that is designed
to replace std::getline(), fgets(), and similar ways of iterating over text
lines.  Key properties:
* Instead of copying every line to your buffer, one at a time, it just returns
  a pointer to the beginning of each line, and gives you access to a pointer to
  the end.  In exchange, the line is invalidated when you iterate to the next
  one; it's like being forced to pass the same string to std::getline(), or the
  same buffer to fgets(), on every call.  But whenever that's problematic, you
  can always just copy the line before iterating.  In the many situations where
  there's no need to copy, you have a substantially lower-cost abstraction than
  the aforementioned standard library functions.
* It automatically detects and decompresses gzipped and Zstd-compressed files.
  This works with streams.
* It automatically reads AND DECOMPRESSES ahead for you.  Decompression is even
  multithreaded by default when the file is BGZF-compressed.

The second library is pgenlib.  This supports reading and writing of PLINK 2.x
genotype files (".pgen").  A draft specification for this format is under
../pgen_spec/ ; here are some key properties:
* A PLINK 1 .bed is a valid .pgen.
* In addition, .pgen can represent multiallelic, phased, and/or dosage
  information.  As of this writing, software support for multiallelic dosages
  does not exist yet, but it does for any other combination of these attributes
  (e.g. multiallelic+phased or phased+dosage).
* .pgen CANNOT represent genotype probability triplets.  It also cannot store
  read depths, per-call quality scores, etc.  While plink2 can *filter* on the
  aforementioned BGEN/VCF fields during import, it cannot re-export or do
  anything else with them.  Use other software, such as bcftools
  (https://samtools.github.io/bcftools/bcftools.html ) or qctool2
  (https://www.well.ox.ac.uk/~gav/qctool_v2/ ) when you must retain any of
  these fields.
* .pgen is compressed, but in a custom manner (based on SNPack) that supports
  very fast compression and decompression.  It is even practical to perform
  several key computations (e.g. allele frequency) directly on the compressed
  representation, and this capability is exposed by the pgenlib library.
* Python/pgenlib.pyx is the Python wrapper (see Python/python_api.txt for
  details), and pgenlibr/ is the R wrapper.  These are somewhat incomplete as
  of this writing, but it would not take much effort to fill in key components;
  that work is scheduled for roughly the time of the beta release, but if you
  could really use a specific feature earlier, you have good odds of getting it
  by asking on plink2-dev.  (plink2-dev is also the place to ask other
  questions about any of this code.)

As for the PLINK 2.0 application:
* build_dynamic/ contains a Makefile suitable for producing Linux and OS X
  dynamic builds.  On Linux, if Intel MKL is installed using the instructions
  at e.g.
  https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo ,
  you can dynamically link to it.
* build_win/ contains a Makefile for producing static Windows builds.  This
  requires MinGW[-w64] and zlib; a prebuilt OpenBLAS package from
  sourceforge.net/projects/openblas/files/ is also strongly recommended.
* More documentation is at www.cog-genomics.org/plink/2.0/ .
