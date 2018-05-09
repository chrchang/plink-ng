
###########################################################################
# Note: This Makefile is intended for development.  For production
# builds, use build_dynamic/Makefile.
###########################################################################

BASEFLAGS=-mavx2 -mbmi -mbmi2 -mlzcnt -DNDEBUG -DZSTD_MULTITHREAD
# BASEFLAGS=-mavx2 -mbmi -mbmi2 -mlzcnt -DZSTD_MULTITHREAD
# BASEFLAGS=-msse4.2 -DZSTD_MULTITHREAD

include Makefile.src

# Respect the environment:
# Use defaults below only if not set in env or make arguments
CC		?= gcc
CXX		?= g++
CFLAGS		?= -O3
CXXFLAGS	?= -O2
ZLIB		?= ../zlib-1.2.11/libz-64.a
BLASFLAGS64	?= -llapack -lf77blas -latlas
ARCH32		?=

# Mandatory flags added to defaults or env settings
CFLAGS		+= -std=gnu99 $(BASEFLAGS) $(CWARN) $(CINCLUDE) $(ZSTD_INCLUDE)
CXXFLAGS	+= -std=c++14 $(BASEFLAGS) $(CXXWARN)
LDFLAGS		+= -lm -lpthread -L. $(ZLIB)

# Installation defaults
MKDIR		?= mkdir
INSTALL		?= install
STRIP_CMD	?= strip
PREFIX		?= /usr/local
DESTDIR		?= .

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  BLASFLAGS=-framework Accelerate
  BLASFLAGS64=-framework Accelerate
  LDFLAGS=-ldl -lpthread -lz
endif

%.o: %.c
	$(CC) -c $(CFLAGS) $(ARCH32) -o $@ $<

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(ARCH32) -o $@ $<

all: plink2 pgen_compress

# for clean build, "make clean" first
# Run mkdir for both plink2 and pgen_compress as we don't know which
# target will run first
plink2: $(OBJ)
	$(MKDIR) -p bin
	$(CXX) $(ARCH32) $(OBJ) -o bin/plink2 $(BLASFLAGS64) $(LDFLAGS)

# basic pgenlib_internal.h usage example; also needed for tests
pgen_compress: plink2_base.o pgenlib_internal.o pgen_compress.o
	$(MKDIR) -p bin
	$(CXX) plink2_base.o pgenlib_internal.o pgen_compress.o \
		-o bin/pgen_compress

.PHONY: install-strip install clean

install-strip: install
	$(STRIP_CMD) $(DESTDIR)$(PREFIX)/bin/*

install: all
	$(MKDIR) -p $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -c bin/* $(DESTDIR)$(PREFIX)/bin

clean:
	rm -f $(CLEAN)
