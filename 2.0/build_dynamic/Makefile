# Linux/macOS Makefile for PLINK 2.00.
#
# Compilation options (leave blank after "=" to disable, put "= 1" to enable):
#   Do not use AVX2 instructions: NO_AVX2
#   Do not use SSE4.2 instructions: NO_SSE42
#   Print clear error message if SSE42/AVX2 needed but missing: CPU_CHECK
#   Do not link to LAPACK: NO_LAPACK
#   Use only -O2 optimization for zstd (may be necessary for gcc 4.x): ZSTD_O2
#   Statically link zlib: STATIC_ZLIB
#   Statically link zstd: STATIC_ZSTD
#   Link to MKL with 64-bit indexes (dynamically): DYNAMIC_MKL
#     (this also requires MKLROOT and MKL_IOMP5_DIR to be defined, and
#     LD_LIBRARY_PATH to include the appropriate directories)
#   Link statically to AOCL with 64-bit indexes: STATIC_AOCL
#     (this supports older OSes than dynamic linking)
#   32-bit binary (also sets STATIC_ZLIB and ZSTD_O2):
#     FORCE_32BIT (warning: you may need to add a zconf.h symlink to make that
#     work)
#   Debug symbols: change DEBUG to -g (still incompatible with -flto?)
NO_AVX2 = 1
NO_SSE42 = 1
CPU_CHECK = 1
NO_LAPACK =
ZSTD_O2 = 1
STATIC_ZLIB =
STATIC_ZSTD = 1
DYNAMIC_MKL =
MKLROOT = /opt/intel/oneapi/mkl/2025.1
MKL_IOMP5_DIR = /opt/intel/oneapi/compiler/2025.1/lib
STATIC_AOCL =
AOCLROOT = /home/admin/aocl/5.1.0/gcc
FORCE_32BIT =
DEBUG = -flto
CC ?= gcc
CXX ?= g++

BASEFLAGS= -flto -DZSTD_MULTITHREAD
# ***** end configuration *****

BASEFLAGS += ${DEBUG}

include ../Makefile.src

LINKFLAGS=-lm -lpthread ${DEBUG}
ZLIB=
ARCH32=
CPUCHECK_FLAGS = ${DEBUG}

ifdef FORCE_32BIT
  # this should work on Debian 10 with the right packages installed.
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
        CPUCHECK_FLAGS = -O2 -flto -DCPU_CHECK_SSE42 ${CXXWARN2}
      endif
    endif
  else
    BASEFLAGS += -mavx2 -mbmi -mbmi2 -mfma -mlzcnt
    ifdef CPU_CHECK
      BASEFLAGS += -DCPU_CHECK_AVX2
      CPUCHECK_FLAGS = -O2 -flto -DCPU_CHECK_AVX2 ${CXXWARN2}
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

############################
# PLINK2 Linux build update:

# Old BLASFLAGS=-llapack -lblas -lcblas -latlas fails on Ubuntu 22.04/24.04 (lcblas/Atlas gone).

# Fixed: BLASFLAGS=-llapacke -llapack -lopenblas

# Install: libopenblas-dev liblapack-dev liblapacke-dev zlib1g-dev

# macOS build unaffected (still uses Accelerate)

# Works cleanly in Docker / GitHub Actions

## Old:
# BLASFLAGS=-llapack -lblas -lcblas -latlas

## Update:
BLASFLAGS=-llapacke -llapack -lopenblas

###########################

ifdef STATIC_ZLIB
  BASEFLAGS += -DSTATIC_ZLIB
  LINKFLAGS += -L. ../../zlib-${ZLIB_VER}/libz.a
else
  LINKFLAGS += -lz
endif

SKIP_STATIC_ZSTD =
ifndef STATIC_ZSTD
  BASEFLAGS += -DIGNORE_BUNDLED_ZSTD
  ZCSRC2 =
  ZSSRC2 =
  SKIP_STATIC_ZSTD = echo
  OBJ = $(CSRC:.c=.o) $(CCSRC:.cc=.o)
  OBJL = $(notdir $(OBJ))
  LINKFLAGS += -lzstd
endif

LIBTOOL=ar
STATIC=rcs
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  ifdef FORCE_32BIT
    $(error 32-bit OS X builds are not supported)
  endif
  ifdef DYNAMIC_MKL
    $(error MKL is not currently supported on OS X)
  endif
  ifdef STATIC_AOCL
    $(error AOCL is not currently supported on OS X)
  endif
  BLASFLAGS=-framework Accelerate
  CFLAGS += -I/usr/local/include
  CXXFLAGS += -I/usr/local/include -stdlib=libc++ -DACCELERATE_NEW_LAPACK
  LINKFLAGS += -L/usr/local/lib
  LIBTOOL=libtool
  STATIC=-static
else
  ifdef DYNAMIC_MKL
    ifdef NO_LAPACK
      $(error DYNAMIC_MKL and NO_LAPACK conflict)
    endif
    ifdef FORCE_32BIT
      $(error DYNAMIC_MKL + FORCE_32BIT not supported)
    endif
    ifdef STATIC_AOCL
      $(error DYNAMIC_MKL + STATIC_AOCL does not make sense)
    endif
    BASEFLAGS += -DDYNAMIC_MKL -DLAPACK_ILP64 -I${MKLROOT}/include
    BLASFLAGS = -L${MKLROOT}/lib/intel64 -L${MKL_IOMP5_DIR} -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
    LINKFLAGS += -ldl
  endif
  ifdef STATIC_AOCL
    ifdef NO_LAPACK
      $(error STATIC_AOCL and NO_LAPACK conflict)
    endif
    ifdef FORCE_32BIT
      $(error STATIC_AOCL + FORCE_32BIT not supported)
    endif
    BASEFLAGS += -DUSE_AOCL -DLAPACK_ILP64 -I${AOCLROOT}/include
    BLASFLAGS = -L. ${AOCLROOT}/lib/libflame.a -L. ${AOCLROOT}/lib/libblis.a -L. ${AOCLROOT}/lib/libau_cpuid.a -lgomp
    LINKFLAGS += -ldl
  endif
endif

ifdef NO_LAPACK
  BASEFLAGS += -DNOLAPACK
  BLASFLAGS=
endif

CFLAGS += ${BASEFLAGS} ${CWARN2} ${CINCLUDE2}
ZCFLAGS += ${BASEFLAGS} ${CWARN2} ${ZSTD_INCLUDE2}
CXXFLAGS += ${BASEFLAGS} ${CXXWARN2} ${CXXINCLUDE2}

ifdef FORCE_32BIT
  CXXFLAGS += -Wno-sign-compare
endif

SFX ?= ""

all: plink2$(SFX) pgen_compress$(SFX)

plink2$(SFX): $(CSRC2) $(ZCSRC2) $(ZSSRC2) $(CCSRC2) ../plink2_cpu.cc
	$(CC) $(CFLAGS) $(CSRC2) -c
	$(SKIP_STATIC_ZSTD) $(CC) $(ZCFLAGS) $(ZCSRC2) $(ZSSRC2) -c
	$(CXX) $(CXXFLAGS) $(CCSRC2) -c
	$(CXX) $(CPUCHECK_FLAGS) ../plink2_cpu.cc -c
	$(CXX) $(OBJL) plink2_cpu.o $(ARCH32) -o $@ $(BLASFLAGS) $(LINKFLAGS)

pgen_compress$(SFX): $(PGCSRC2)
	$(CXX) $(CXXFLAGS) $(PGCSRC2) -o $@

pgenlib.a: $(PGENLIB_CCSRC2)
	$(CXX) $(CXXFLAGS) $(PGENLIB_CCSRC2) -c
	$(LIBTOOL) $(STATIC) -o $@ $(PGENLIB_OBJL)

plink2lib.a: $(CSRC2) $(ZCSRC2) $(ZSSRC2) $(PLINK2LIB_CCSRC2)
	$(CC) $(CFLAGS) $(CSRC2) -c
	$(CC) $(ZCFLAGS) $(ZCSRC2) $(ZSSRC2) -c
	$(CXX) $(CXXFLAGS) $(PLINK2LIB_CCSRC2) -c
	$(LIBTOOL) $(STATIC) -o $@ $(PLINK2LIB_OBJL)

static_pgenlib_test: pgenlib.a ../pgen_compress.cc
	$(CXX) $(CXXFLAGS) -o $@ ../pgen_compress.cc -L. pgenlib.a

static_plink2lib_test: plink2lib.a ../pgen_compress.cc
	$(CXX) $(CXXFLAGS) -o $@ ../pgen_compress.cc -L. plink2lib.a

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink2
	rm -f pgen_compress
	rm -f pgenlib.a
	rm -f plink2lib.a
	rm -f static_pgenlib_test
	rm -f static_plink2lib_test
