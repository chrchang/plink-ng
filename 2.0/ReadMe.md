This is a library for reading and querying PLINK 2.x genotype files (".pgen"),
along with an alpha version of PLINK 2.0 built on top of the library.  A few
notes:

* pgenlib_internal.h/pgenlib_internal.cpp is lowest-common-denominator code
  designed to compile on every conceivably relevant platform (e.g. Windows;
  incompatibilities with ARM are also flagged now).  It has three modes:
    1. fread()-based one-variant-at-a-time,
    2. mmap()-based one-variant-at-a-time, and
    3. fread()-based large blocks (most commonly 64k variants at a time).
* Python/pgenlib.pyx is the Python wrapper.  See Python/python_api.txt for
  details.
* Instead of 00 = hom minor, 01 = missing, 10 = het, 11 = hom major, the basic
  2-bit encoding is 00 = hom ref, 01 = het ref/alt1, 10 = hom alt1,
  11 = missing/other.
  PLINK 1 .bed files are grandfathered in as valid .pgen files; there is a tiny
  efficiency penalty associated with accessing them via this library (since,
  after the data is loaded, it is rotated to the new 2-bit encoding).  The
  PgrPlink2ToPlink1InplaceUnsafe() function can be used to rotate the new
  representation back to PLINK 1-format, to simplify porting of existing code.
* The rest of the core format is essentially an extension of SNPack ideas to
  multiallelic data.
* The variant and sample counts are now stored in the file (unless it's a
  .bed).  Also, there are two optional data tracks in the header:
    1. Alt allele counts.  Presumably, these can also be inferred from the
       variant info file, but occasionally it may be convenient to work with
       the .pgen without loading all the variant info.
    2. Untrusted reference allele flags.  If you merge .bed and VCF data, it's
       nice to distinguish reference alleles known from the VCF vs.
       reference-allele-unknown variants.
* Phased data is now supported; see UNIT_TEST_PHASED_VCF for plink2 usage, and
  the Python interface's read_phased() function.
* Dosage data is also supported.  It can be phased.

build_dynamic/ contains a Makefile suitable for producing Linux and OS X
dynamic builds.  On Linux, if Intel MKL is installed using the instructions at
e.g.
https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo ,
you can dynamically link to it.

build_win/ contains a Makefile for producing static Windows builds.  This
requires MinGW[-w64] and zlib; a prebuilt OpenBLAS package from
sourceforge.net/projects/openblas/files/ is also strongly recommended.

More documentation is at www.cog-genomics.org/plink/2.0/ .