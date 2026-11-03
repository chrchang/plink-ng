# Linux/MacOS Makefile for vcf_subset.

NO_AVX2 = 1
NO_SSE42 = 1
ZSTD_O2 = 1
STATIC_ZLIB =
STATIC_ZSTD ?= 1
FORCE_32BIT =
DEBUG =
CC ?= gcc
CXX ?= g++

BASEFLAGS=
# ***** end configuration *****

BASEFLAGS += ${DEBUG}

include ../../Makefile.src
OBJ = $(CSRC:.c=.o) $(ZCSRC:.c=.o) $(ZSSRC:.S=.o) $(TEXTLIB_CCSRC:.cc=.o)
OBJL = $(notdir $(OBJ))

CFLAGS ?=
CXXFLAGS ?=
LDFLAGS ?=
LDLIBS=-lpthread ${DEBUG}
ZLIB=
ARCH32=

ifdef FORCE_32BIT
  # this should work on Debian 10 with the right packages installed.
  STATIC_ZLIB = 1
  ZSTD_O2 = 1
  BASEFLAGS += -msse -mfpmath=sse
  ARCH32 = -m32 -march=i686
  CXXFLAGS += -std=c++17
else
  ifdef NO_AVX2
    ifndef NO_SSE42
      BASEFLAGS += -msse4.2
    endif
  else
    BASEFLAGS += -mavx2 -mbmi -mbmi2 -mfma -mlzcnt
  endif
  CXXFLAGS += -std=c++17
endif
BASEFLAGS += ${ARCH32}

CFLAGS += -O2 -std=gnu99
# zstd appears to be seriously targeted at -O3; see 26 Jul 2016 entry at
# cbloom.com/rants.html
ZCFLAGS = ${CFLAGS}
ifndef ZSTD_O2
  ZCFLAGS += -O3
endif
# this actually needs to be named "CXXFLAGS"
CXXFLAGS += -O2

ifdef STATIC_ZLIB
  BASEFLAGS += -DSTATIC_ZLIB
  LDFLAGS += -L.
  LDLIBS += ../../zlib-${ZLIB_VER}/libz.a
else
  LDLIBS += -lz
endif

SKIP_STATIC_ZSTD =
ifndef STATIC_ZSTD
  BASEFLAGS += -DIGNORE_BUNDLED_ZSTD
  ZCSRC3 =
  ZSSRC3 =
  SKIP_STATIC_ZSTD = echo
  TEXTLIB_OBJ = $(CSRC:.c=.o) $(TEXTLIB_CCSRC:.cc=.o)
  TEXTLIB_OBJL = $(notdir $(TEXTLIB_OBJ))
  LDLIBS += -lzstd
endif

LIBTOOL=ar
STATIC=rcs
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  ifdef FORCE_32BIT
    $(error 32-bit MacOS builds are not supported)
  endif
  CFLAGS += -I/usr/local/include
  CXXFLAGS += -I/usr/local/include -stdlib=libc++
  LDFLAGS += -L/usr/local/lib
  LIBTOOL=libtool
  STATIC=-static
endif

CFLAGS += ${BASEFLAGS} ${CWARN2} ${CINCLUDE2}
ZCFLAGS += ${BASEFLAGS} ${CWARN2} ${ZSTD_INCLUDE2}
CXXFLAGS += ${BASEFLAGS} ${CXXWARN2} ${CXXINCLUDE2}

ifdef FORCE_32BIT
  CXXFLAGS += -Wno-sign-compare
endif

SFX ?= ""

all: vcf_subset$(SFX)

vcf_subset$(SFX): $(CSRC3) $(ZCSRC3) $(ZSSRC3) $(TEXTLIB_CCSRC3) vcf_subset.cc
	$(CC) $(CFLAGS) $(CSRC3) -c
	$(SKIP_STATIC_ZSTD) $(CC) $(ZCFLAGS) $(ZCSRC3) $(ZSSRC3) -c
	$(CXX) $(CXXFLAGS) $(TEXTLIB_CCSRC3) vcf_subset.cc -c
	$(CXX) $(LDFLAGS) $(TEXTLIB_OBJL) vcf_subset.o $(ARCH32) -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f ../../include/*.o
	rm -f vcf_subset
