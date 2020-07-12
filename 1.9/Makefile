# General-purpose Makefile for PLINK 1.90
#
# Compilation options:
#   Do not link to LAPACK                    NO_LAPACK

# Leave blank after "=" to disable; put "= 1" to enable
# (when enabled, "#define NOLAPACK" must be uncommented in plink_common.h)
NO_LAPACK =

# Variable that allows additional linker flags (e.g., "-L/path/to/libs") to be
# passed into the linker; to use this, invoke 'make LDFLAGS_EXTRA="<opts>"'.
LDFLAGS_EXTRA =

.PHONY: clean install

# should autodetect system
SYS = UNIX
ifdef SystemRoot
  SYS = WIN
else
  UNAME := $(shell uname)
  ifeq ($(UNAME), Darwin)
    SYS = MAC
  endif
endif

# Allow these to be overridden by make arguments or env variables, so people
# don't have to edit the Makefile to build in a different environment.
BIN ?=		plink
CC ?=		gcc
CXX ?=		g++
CFLAGS ?=	-Wall -O2
CXXFLAGS ?=	-Wall -O2
MACFLAGS ?=	-mmacosx-version-min=10.9

PREFIX ?=	/usr/local
DESTDIR ?=	.
INSTALL ?=	install
STRIP ?=	strip

ifeq ($(SYS), MAC)
  GCC_GTEQ_43 := $(shell expr `g++ -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40300)
  ifeq "$(GCC_GTEQ_43)" "1"
    CFLAGS ?= -Wall -O2 -flax-vector-conversions
  endif
  CFLAGS += $(MACFLAGS)
  CXXFLAGS += $(MACFLAGS)
  BLASFLAGS ?= -framework Accelerate
  LDFLAGS ?= -ldl
  ZLIB ?= ../zlib-1.2.11/libz.1.2.11.dylib
endif

ifeq ($(SYS), WIN)
# Note that, unlike the Linux and Mac build processes, this STATICALLY links
# LAPACK, since we have not gotten around to trying dynamically-linked LAPACK
# on Windows.
# If you don't already have LAPACK built, you'll probably want to turn on
# NO_LAPACK.
  BLASFLAGS ?= -L. lapack/liblapack.a -L. lapack/librefblas.a
  LDFLAGS ?= -lm -static-libgcc
  ZLIB ?= ../zlib-1.2.11/libz.a
endif

# These must appear after the MAC/WIN-specific ?= statements.
BLASFLAGS ?=	-L/usr/lib64/atlas -llapack -lblas -lcblas -latlas
LDFLAGS ?=	-lm -lpthread -ldl
ZLIB ?=		../zlib-1.2.11/libz.so.1.2.11

ifdef NO_LAPACK
  BLASFLAGS=
endif

OBJS = plink.o plink_assoc.o plink_calc.o plink_cluster.o plink_cnv.o plink_common.o plink_data.o plink_dosage.o plink_family.o plink_filter.o plink_glm.o plink_help.o plink_homozyg.o plink_lasso.o plink_ld.o plink_matrix.o plink_misc.o plink_perm.o plink_rserve.o plink_set.o plink_stats.o SFMT.o dcdflib.o pigz.o yarn.o Rconnection.o hfile.o bgzf.o

# In the event that you are still concurrently using PLINK 1.07, we suggest
# renaming that binary to "plink107", and this one to "plink" or "plink1".

all:	$(BIN)

plink: $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS_EXTRA) $(BLASFLAGS) $(LDFLAGS) -L. $(ZLIB) -o $@

plinkw: $(OBJS)
	gfortran -O2 $(OBJS) $(LDFLAGS_EXTRA) -Wl,-Bstatic $(BLASFLAGS) $(LDFLAGS) -L. $(ZLIB) -o $@

install:
	$(INSTALL) -d $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -c $(BIN) $(DESTDIR)$(PREFIX)/bin

install-strip: install
	$(STRIP) $(DESTDIR)$(PREFIX)/bin/$(BIN)

clean:
	rm -f $(OBJS) plink plinkw

# Pattern-based rules for compiling object (.o) files; basically identical to
# GNU make's built-in rules, except we explicitly use "g++" instead of $(CC).

# Compiling C files with a C++ compiler is deprecated, but the code needs
# to be cleaned up before we can switch to cc.  E.g. plink_rserve.c
# includes a C++ header and exposed functions will need to be declared with
# extern "C".
%.o: %.c
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

%.o: %.cc
	$(CXX) -x c++ -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
