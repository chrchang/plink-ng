# Linux/OS X Makefile for PLINK 2.00.
#
# Compilation options (leave blank after "=" to disable, put "= 1" to enable):
#   Do not link to LAPACK: NO_LAPACK
#   Use cblas_f77 instead of cblas: FORCE_CBLAS_F77
#   Use only -O2 optimization for zstd (may be necessary for gcc 4.x): ZSTD_O2
#   Statically link zlib: STATIC_ZLIB
#   Link to MKL with 64-bit indexes (dynamically): DYNAMIC_MKL
#     (this also requires MKLROOT and MKL_IOMP5_DIR to be defined, and
#     LD_LIBRARY_PATH to include the appropriate directories)
#   32-bit binary (also sets STATIC_ZLIB and ZSTD_O2):
#     FORCE_32BIT (warning: you may need to add a zconf.h symlink to make that
#     work)
NO_AVX2 = 1
NO_SSE42 =
NO_LAPACK =
PREFER_CBLAS_F77 =
ZSTD_O2 = 1
STATIC_ZLIB =
DYNAMIC_MKL =
MKLROOT = /home/ubuntu/intel/mkl
MKL_IOMP5_DIR = /home/ubuntu/intel/compilers_and_libraries_2017.2.174/linux/compiler/lib/intel64
FORCE_32BIT =

BASEFLAGS=-Wall -Wextra -Wshadow -DZSTD_MULTITHREAD
# ***** end configuration *****

LINKFLAGS=-lm -lpthread
ZLIB=
ARCH32=

ifdef FORCE_32BIT
  # this is targeted at Scientific Linux 6.
  STATIC_ZLIB = 1
  ZSTD_O2 = 1
  ARCH32 = -m32 -march=i686
  CXXFLAGS = -std=c++0x
else
  ifdef NO_AVX2
    ifndef NO_SSE42
      BASEFLAGS += -msse4.2
    endif
  else
    BASEFLAGS += -mavx2 -mbmi -mbmi2 -mlzcnt
  endif
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
    BASEFLAGS += -DDYNAMIC_MKL -DLAPACK_ILP64 -I${MKLROOT}/include
    BLASFLAGS = -L${MKLROOT}/lib/intel64 -L${MKL_IOMP5_DIR} -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
    LINKFLAGS += -ldl
  endif
endif

ifdef NO_LAPACK
  BASEFLAGS += -DNOLAPACK
  BLASFLAGS=
endif

CFLAGS += ${BASEFLAGS}
ZCFLAGS += ${BASEFLAGS}
CXXFLAGS += ${BASEFLAGS}

include ../Makefile.src
ZCFLAGS += ${ZSTD_INCLUDE2}

all: plink2 pgen_compress

plink2: $(CSRC2) $(ZCSRC2) $(CCSRC2)
	gcc $(CFLAGS) $(CSRC2) -c
	gcc $(ZCFLAGS) $(ZCSRC2) -c
	g++ $(CXXFLAGS) $(CCSRC2) -c
	g++ $(OBJ2) $(ARCH32) -o plink2 $(BLASFLAGS) $(LINKFLAGS)

pgen_compress: ../plink2_base.cc ../pgenlib_internal.cc ../pgen_compress.cc
	g++ $(CXXFLAGS) ../plink2_base.cc ../pgenlib_internal.cc ../pgen_compress.cc -o pgen_compress

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink2
	rm -f pgen_compress
