# Linux/OS X Makefile for PLINK 2.00.
#
# Compilation options (leave blank after "=" to disable, put "= 1" to enable):
#   Do not use AVX2 instructions: NO_AVX2
#   Do not use SSE4.2 instructions: NO_SSE42
#   Print clear error message if SSE42/AVX2 needed but missing: CPU_CHECK
#   Do not link to LAPACK: NO_LAPACK
#   Use cblas_f77 instead of cblas: FORCE_CBLAS_F77
#   Use only -O2 optimization for zstd (may be necessary for gcc 4.x): ZSTD_O2
#   Statically link zlib: STATIC_ZLIB
#   Statically link zstd: STATIC_ZSTD
#   Link to MKL with 64-bit indexes (dynamically): DYNAMIC_MKL
#     (this also requires MKLROOT and MKL_IOMP5_DIR to be defined, and
#     LD_LIBRARY_PATH to include the appropriate directories)
#   32-bit binary (also sets STATIC_ZLIB and ZSTD_O2):
#     FORCE_32BIT (warning: you may need to add a zconf.h symlink to make that
#     work)
#   Debug symbols: set DEBUG to -g

AVX2_CUDA = 1
NO_AVX2 =
NO_SSE42 =
CPU_CHECK = 1
NO_LAPACK =
PREFER_CBLAS_F77 =
ZSTD_O2 = 1
STATIC_ZLIB =
STATIC_ZSTD = 1
# export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64 needed in .profile
DYNAMIC_MKL = 1
MKLPREROOT = /opt/intel
CUDAROOT = /usr/local/cuda
FORCE_32BIT =
DEBUG =
CC ?= gcc
CXX ?= g++
BASEFLAGS=-DZSTD_MULTITHREAD
# ***** end configuration *****

BASEFLAGS += ${DEBUG}
MKLROOT = ${MKLPREROOT}/mkl
MKL_IOMP5_DIR = ${MKLPREROOT}/lib/intel64

include ../Makefile.src

LINKFLAGS=-lm -lpthread ${DEBUG}
ZLIB=
ARCH32=
CPUCHECK_FLAGS = ${DEBUG}

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
      ifdef CPU_CHECK
        BASEFLAGS += -DCPU_CHECK_SSE42
        CPUCHECK_FLAGS = -O2 -DCPU_CHECK_SSE42 ${CXXWARN2}
      endif
    endif
  else
    BASEFLAGS += -mavx2 -mbmi -mbmi2 -mlzcnt
    ifdef AVX2_CUDA
      BASEFLAGS += -DUSE_CUDA
      LINKFLAGS += -L${CUDAROOT}/lib64 -lcudart -lcublas
    endif
    ifdef CPU_CHECK
      BASEFLAGS += -DCPU_CHECK_AVX2
      CPUCHECK_FLAGS = -O2 -DCPU_CHECK_AVX2 ${CXXWARN2}
    endif
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

SKIP_STATIC_ZSTD =
ifdef STATIC_ZSTD
  BASEFLAGS += -DSTATIC_ZSTD
else
  ZCSRC2 =
  SKIP_STATIC_ZSTD = echo
  OBJ = $(CSRC:.c=.o) $(CCSRC:.cc=.o)
  OBJ2 = $(notdir $(OBJ))
  LINKFLAGS += -lzstd
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
  CXXFLAGS += -stdlib=libc++
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

CFLAGS += ${BASEFLAGS} ${CWARN2} ${CINCLUDE2}
ZCFLAGS += ${BASEFLAGS} ${CWARN2} ${ZSTD_INCLUDE2}
CXXFLAGS += ${BASEFLAGS} ${CXXWARN2}

ifdef FORCE_32BIT
  CXXFLAGS += -Wno-sign-compare
endif

all: plink2 pgen_compress

plink2: $(CSRC2) $(ZCSRC2) $(CCSRC2) ../plink2_cpu.cc ../cuda/plink2_matrix_cuda.cu
	$(CC) $(CFLAGS) $(CSRC2) -c
	$(SKIP_STATIC_ZSTD) gcc $(ZCFLAGS) $(ZCSRC2) -c
	nvcc -cudart shared ../cuda/plink2_matrix_cuda.cu -c
	$(CXX) $(CXXFLAGS) $(CCSRC2) -c
	$(CXX) $(CPUCHECK_FLAGS) ../plink2_cpu.cc -c
	$(CXX) $(OBJ2) plink2_cpu.o plink2_matrix_cuda.o $(ARCH32) -o plink2 $(BLASFLAGS) $(LINKFLAGS)

pgen_compress: $(PGCSRC2)
	$(CXX) $(CXXFLAGS) $(PGCSRC2) -o pgen_compress

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink2
	rm -f pgen_compress
