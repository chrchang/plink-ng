This is a library for reading and querying PLINK 2.x genotype files (".pgen"),
along with an alpha version of PLINK 2.0 built on top of the library.  A few
notes:

* The draft .pgen specification is under ../pgen_spec/ .  Important properties:
    1. A PLINK 1 .bed is a valid .pgen.
    2. In addition, .pgen can represent multiallelic, phased, and/or dosage
       information.  As of this writing, software support for multiallelic
       dosages does not exist yet, but it does for any other combination of
       these attributes (e.g. multiallelic+phased or phased+dosage).
    3. .pgen CANNOT represent genotype probability triplets.  It also cannot
       store read depths, per-call quality scores, etc.  While plink2 can
       *filter* on the aforementioned BGEN/VCF/BCF fields during import, it
       cannot re-export or do anything else with them.  Use other software,
       such as bcftools (https://samtools.github.io/bcftools/bcftools.html ) or
       qctool2 (https://www.well.ox.ac.uk/~gav/qctool_v2/ ) when you must
       retain any of these fields.
    4. .pgen is compressed, but in a manner that supports very fast compression
       and decompression.  It is even practical to perform several key
       computations (e.g. allele frequency) directly on the compressed
       representation, and this capability is exposed by the pgenlib library.
* The same specification document describes the .pvar (variant info) and .psam
  (sample info) formats understood by plink2.  The most important thing to know
  is that a typical VCF file is a valid .pvar, and so is a PLINK 1 .bim file.
* pgenlib_{misc,read,write}.{h,cpp} is lowest-common-denominator code designed
  to compile on every plausibly relevant platform (e.g. Windows;
  incompatibilities with ARM are also flagged now).  It has three modes:
    1. fread()-based one-variant-at-a-time,
    2. mmap()-based one-variant-at-a-time (warning: not tested recently), and
    3. fread()-based large blocks (most commonly 64k variants at a time).
* Python/pgenlib.pyx is the Python wrapper (see Python/python_api.txt for
  details), and pgenlibr/ is the R wrapper.  These are somewhat incomplete as
  of this writing, but it would not take much effort to fill in key components;
  that work is scheduled for roughly the time of the beta release, but if you
  could really use a specific feature before then you have good odds of getting
  it by asking on plink2-users.

build_dynamic/ contains a Makefile suitable for producing Linux and OS X
dynamic builds.  On Linux, if Intel MKL is installed using the instructions at
e.g.
https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo ,
you can dynamically link to it.

build_win/ contains a Makefile for producing static Windows builds.  This
requires MinGW[-w64] and zlib; a prebuilt OpenBLAS package from
sourceforge.net/projects/openblas/files/ is also strongly recommended.

More documentation is at www.cog-genomics.org/plink/2.0/ .