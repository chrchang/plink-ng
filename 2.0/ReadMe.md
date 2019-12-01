The include/ subdirectory contains two (LGPL3-licensed) major libraries, while
the immediate directory contains the PLINK 2.0 application built on top of
them.  These are carefully written to be valid C99 (from gcc and clang's
perspective, anyway) to simplify FFI development, while still taking advantage
of quite a few C++-specific affordances to improve safety and occasionally
performance.  They are currently x86-specific, but there are annotations to
facilitate a possible future port to ARM.

The first library is plink2_text, which provides a pair of classes designed to
replace std::getline(), fgets(), and similar ways of iterating over text lines.
Key properties:
* Instead of copying every line to your buffer, one at a time, these classes
  just return a pointer to the beginning of each line in the underlying binary
  stream, and give you access to a pointer to the end.  In exchange, the line
  is invalidated when you iterate to the next one; it's like being forced to
  pass the same string to std::getline(), or the same buffer to fgets(), on
  every call.  But whenever that's problematic, you can always copy the line
  before iterating to the next; on all systems I've seen, that *still* exhibits
  better throughput than getline/fgets.  And in the many situations where
  there's no need to copy, you get a fundamentally lower-latency abstraction.
* They automatically detect and decompress gzipped and Zstd-compressed
  (https://facebook.github.io/zstd/ ) files, in a manner that works with pipe
  file descriptors.
* The primary TextStream class automatically reads *and decompresses* ahead for
  you.  Decompression is even multithreaded by default when the file is
  BGZF-compressed.  (And the textFILE class covers the setting where you don't
  want to launch any more threads.)
* They do not support network input as of this writing, but that would not be
  difficult to add.  The existing code uses FILE* in a very straightforward
  manner.
* As for text parsing, the ScanadvDouble() utility function in the
  plink2_string component is a very efficient string-to-double converter.
  While it does not support perfect string<->double round-trips (that's what
  C++17 std::from_chars is for; https://abseil.io/ has a working implementation
  while we wait for gcc/clang...), or long-tail features like locale-specific
  decimal separators or hex floats, it has been incredibly useful for speeding
  up the basic job of scanning standard-locale printf("%g")-formatted and
  similar output.  (Note that you lose roughly a billion times as much accuracy
  to %g's 6-digit limit as you do to imperfect string->double conversion in
  that setting.)

(Coming soon: example text-processing programs using plink2_text.)

The second library is pgenlib.  This supports reading and writing of PLINK 2.x
genotype files (".pgen").  A draft specification for this format is under
https://github.com/chrchang/plink-ng/tree/master/pgen_spec ; here are some key
properties:
* A PLINK 1 .bed is a valid .pgen.
* In addition, .pgen can represent multiallelic, phased, and/or dosage
  information.  As of this writing, software support for multiallelic dosages
  does not exist yet, but it does for the other attribute pairs
  (multiallelic+phased, phased+dosage).
* **.pgen CANNOT represent genotype probability triplets.  It also cannot store
  read depths, per-call quality scores, etc.**  While plink2 can *filter* on
  the aforementioned BGEN/VCF fields during import, it cannot re-export or do
  anything else with them.  Use other software, such as bcftools
  (https://samtools.github.io/bcftools/bcftools.html ) or qctool2
  (www.well.ox.ac.uk/~gav/qctool_v2/ ) when you must retain any of these
  fields.
* .pgen is compressed, but in a domain-specific manner that supports very fast
  compression and decompression.  It is even practical to perform several key
  computations (e.g. allele frequency) directly on the compressed
  representation, and this capability is exposed by the pgenlib library.
* Python/pgenlib.pyx is the Python wrapper (see Python/python_api.txt for
  details), and pgenlibr/ is the R wrapper.  These are somewhat incomplete as
  of this writing, but it would not take much effort to fill in key components;
  that work is scheduled for roughly the time of the beta release, but if you
  could really use a specific feature earlier, you have good odds of getting it
  by asking at https://groups.google.com/forum/#!forum/plink2-dev .
  (plink2-dev is also the place to ask other questions about any of this code.)

As for the PLINK 2.0 application:
* build_dynamic/ contains a Makefile suitable for Linux and macOS dynamic
  builds.  On Linux, if Intel MKL is installed using the instructions at e.g.
  https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo ,
  you can dynamically link to it.
* build_win/ contains a Makefile for producing static Windows builds.  This
  requires MinGW[-w64] and zlib; a prebuilt OpenBLAS package from
  https://sourceforge.net/projects/openblas/files/ is also strongly
  recommended.
* GPUs are not exploited, and there are currently no plans to write much
  GPU-specific code before PLINK 2.0's core function set is completed around
  2021.  However, a few linear-algebra-heavy workloads may benefit
  significantly from a simple replacement of Intel MKL by cuBLAS + cuSOLVER.
  This can probably be supported earlier; feel free to open a GitHub issue
  about it if it would make a big difference to you.
* The LGPL3-licensed plink2_stats component may be of independent interest.  It
  includes a function for computing the 2x2 Fisher's exact test p-value in
  approximately O(sqrt(n)) time--much faster than the O(n) algorithms employed
  by other libraries as of this writing--as well as several log-p-value
  computations (Z-score/chi-square, T-test, F-test) that remain accurate well
  beyond the limits of most other statistical library functions.  (No, you
  don't want to take a 10^{-1000000} p-value literally, but it can be useful to
  distinguish it from 10^{-325}, and both of these numbers can naturally arise
  when analyzing biobank-scale data.)
* More documentation is at www.cog-genomics.org/plink/2.0/ .
