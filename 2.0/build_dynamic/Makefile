# Linux/OS X Makefile for PLINK 2.00.
#
# Compilation options (leave blank after "=" to disable, put "= 1" to enable):
#   Do not use SSE 4.2 instructions: NO_SSE42
#   Do not link to LAPACK: NO_LAPACK
#   Use cblas_f77 instead of cblas: FORCE_CBLAS_F77
#   Use only -O2 optimization for zstd (may be necessary for gcc 4.x): ZSTD_O2
#   Statically link zlib: STATIC_ZLIB
#   Link to MKL with 64-bit indexes (dynamically): DYNAMIC_MKL
#     (this also requires MKLROOT and MKL_IOMP5_DIR to be defined, and
#     LD_LIBRARY_PATH to include the appropriate directories)
#   32-bit binary (also sets STATIC_ZLIB, ZSTD_O2, and NO_SSE42): FORCE_32BIT
#     (warning: you may need to add a zconf.h symlink to make that work)
NO_SSE42 =
NO_LAPACK =
PREFER_CBLAS_F77 =
ZSTD_O2 = 1
STATIC_ZLIB =
DYNAMIC_MKL =
MKLROOT = /home/ubuntu/intel/mkl
MKL_IOMP5_DIR = /home/ubuntu/intel/compilers_and_libraries_2017.2.174/linux/compiler/lib/intel64
FORCE_32BIT =

BASEFLAGS=-Wall -Wextra
# ***** end configuration *****

LINKFLAGS=-lm -lpthread
ZLIB=
ARCH32=

ifdef FORCE_32BIT
  # this is targeted at Scientific Linux 6.
  NO_SSE42 = 1
  STATIC_ZLIB = 1
  ZSTD_O2 = 1
  ARCH32 = -m32 -march=i686
  CXXFLAGS = -std=c++0x
else
  CXXFLAGS = -std=c++11
endif
BASEFLAGS += ${ARCH32}

CFLAGS=-O2 -std=gnu99
# zstd appears to be seriously targeted at -O3; see 26 Jul 2016 entry at
# cbloom.com/rants.html
ifdef ZSTD_O2
  ZCFLAGS=-O2 -std=gnu99
else
  ZCFLAGS=-O3 -std=gnu99
endif
# this actually needs to be named "CXXFLAGS"
CXXFLAGS += -O2

ifndef NO_SSE42
  BASEFLAGS += -msse4.2
endif

ifdef FORCE_CBLAS_F77
  BASEFLAGS += -DFORCE_CBLAS_F77
  BLASFLAGS=-llapack -lf77blas -latlas
else
  BLASFLAGS=-llapack -lblas -lcblas -latlas
endif

ifdef STATIC_ZLIB
  BASEFLAGS += -DSTATIC_ZLIB
  LINKFLAGS += -L. ../../zlib-1.2.11/libz.a
else
  LINKFLAGS += -lz
endif

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  ifdef FORCE_32BIT
    $(error 32-bit OS X builds are not supported)
  endif
  ifdef DYNAMIC_MKL
    $(error MKL is not currently supported on OS X)
  endif
  BLASFLAGS=-framework Accelerate
else
  ifdef DYNAMIC_MKL
    ifdef NO_LAPACK
      $(error DYNAMIC_MKL and NO_LAPACK conflict)
    endif
    ifdef FORCE_32BIT
      $(error DYNAMIC_MKL + FORCE_32BIT not supported)
    endif
    BASEFLAGS = -DDYNAMIC_MKL -DLAPACK_ILP64 -I${MKLROOT}/include
    BLASFLAGS = -L${MKLROOT}/lib/intel64 -L${MKL_IOMP5_DIR} -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
    LINKFLAGS += -ldl
  endif
endif

ifdef NO_LAPACK
  BASEFLAGS += -DNOLAPACK
  BLASFLAGS=
endif

ZSTD_INCLUDE = -I../zstd/lib -I../zstd/lib/common -I../zstd/zlibWrapper
ZCFLAGS += ${ZSTD_INCLUDE}

CFLAGS += ${BASEFLAGS}
ZCFLAGS += ${BASEFLAGS}
CXXFLAGS += ${BASEFLAGS}

CSRC = ../SFMT.c ../hfile.c ../bgzf.c
ZCSRC = ../zstd/zlibWrapper/zstd_zlibwrapper.c ../zstd/zlibWrapper/gzclose.c ../zstd/zlibWrapper/gzlib.c ../zstd/zlibWrapper/gzread.c ../zstd/zlibWrapper/gzwrite.c ../zstd/lib/common/entropy_common.c ../zstd/lib/common/zstd_common.c ../zstd/lib/common/error_private.c ../zstd/lib/common/xxhash.c ../zstd/lib/common/fse_decompress.c ../zstd/lib/compress/fse_compress.c ../zstd/lib/compress/huf_compress.c ../zstd/lib/compress/zstd_compress.c ../zstd/lib/decompress/huf_decompress.c ../zstd/lib/decompress/zstd_decompress.c
CPPSRC = ../pgenlib_internal.cpp ../plink2.cpp ../plink2_adjust.cpp ../plink2_common.cpp ../plink2_compress_stream.cpp ../plink2_data.cpp ../plink2_decompress.cpp ../plink2_filter.cpp ../plink2_glm.cpp ../plink2_help.cpp ../plink2_ld.cpp ../plink2_matrix.cpp ../plink2_matrix_calc.cpp ../plink2_misc.cpp ../plink2_psam.cpp ../plink2_pvar.cpp ../plink2_random.cpp ../plink2_set.cpp ../plink2_stats.cpp

OBJ = SFMT.o hfile.o bgzf.o zstd_zlibwrapper.o gzclose.o gzlib.o gzread.o gzwrite.o entropy_common.o zstd_common.o error_private.o xxhash.o fse_decompress.o fse_compress.o huf_compress.o zstd_compress.o huf_decompress.o zstd_decompress.o pgenlib_internal.o plink2.o plink2_adjust.o plink2_common.o plink2_compress_stream.o plink2_data.o plink2_decompress.o plink2_filter.o plink2_glm.o plink2_help.o plink2_ld.o plink2_matrix.o plink2_matrix_calc.o plink2_misc.o plink2_psam.o plink2_pvar.o plink2_random.o plink2_set.o plink2_stats.o

all: plink2 pgen_compress

plink2: $(CSRC) $(ZCSRC) $(CPPSRC)
	gcc $(CFLAGS) $(CSRC) -c
	gcc $(ZCFLAGS) $(ZCSRC) -c
	g++ $(CXXFLAGS) $(CPPSRC) -c
	g++ $(OBJ) $(ARCH32) -o plink2 $(BLASFLAGS) $(LINKFLAGS)

pgen_compress: ../pgenlib_internal.cpp ../pgen_compress.cpp
	g++ $(CXXFLAGS) ../pgenlib_internal.cpp ../pgen_compress.cpp -o pgen_compress

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink2
	rm -f pgen_compress
