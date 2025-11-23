# Mac Universal (amd64+arm64) Makefile for PLINK 1.9

.PHONY: clean install

SYS = MAC

include ../Makefile.src
# Allow these to be overridden by make arguments or env variables, so people
# don't have to edit the Makefile to build in a different environment.
CC ?=		gcc
CXX ?=		g++
CCFLAGS ?=	-Wall -Wstrict-aliasing -O2
CXXFLAGS ?=	-Wall -Wstrict-aliasing -O2 -I../../2.0/simde "-I../../zlib-"$(ZLIB_VER)
AMDFLAGS ?=     --target=x86_64-apple-macosx10.9
ARMFLAGS ?=     --target=arm64-apple-macos11

PREFIX ?=	/usr/local
DESTDIR ?=	.
INSTALL ?=	install
STRIP ?=	strip

ifeq ($(SYS), MAC)
  GCC_GTEQ_43 := $(shell expr `g++ -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40300)
  ifeq "$(GCC_GTEQ_43)" "1"
    CFLAGS ?= -Wall -Wstrict-aliasing -O2 -flax-vector-conversions
  endif
  CFLAGS += $(MACFLAGS)
  CXXFLAGS += $(MACFLAGS)
  BLASFLAGS ?= -framework Accelerate
  LDFLAGS ?= -ldl
  ZLIB_AMD ?= -L. "../../zlib-"$(ZLIB_VER)"/libz-amd64.a"
  ZLIB_ARM ?= -L. "../../zlib-"$(ZLIB_VER)"/libz-arm64.a"
endif

CCSRC = plink.c plink_assoc.c plink_calc.c plink_cluster.c plink_cnv.c plink_common.c plink_data.c plink_dosage.c plink_family.c plink_filter.c plink_glm.c plink_help.c plink_homozyg.c plink_lasso.c plink_ld.c plink_matrix.c plink_misc.c plink_perm.c plink_rserve.c plink_set.c plink_stats.c SFMT.c dcdflib.c pigz.c yarn.c Rconnection.cc hfile.c bgzf.c
CCSRC2 = $(foreach fname,$(CCSRC),../$(fname))

CCOBJ = $(CCSRC:.cc=.o)

all: plink prettify

plink_amd64: $(CCSRC2)
	$(CXX) $(CXXFLAGS) $(AMDFLAGS) $(CCSRC2) -o plink_amd64 $(BLASFLAGS) $(LDFLAGS) $(ZLIB_AMD)

plink_arm64: $(CCSRC2)
	$(CXX) $(CXXFLAGS) $(ARMFLAGS) $(CCSRC2) -o plink_arm64 $(BLASFLAGS) $(LDFLAGS) $(ZLIB_ARM)

plink: plink_amd64 plink_arm64
	lipo -create -output plink plink_amd64 plink_arm64

prettify_amd64: ../prettify.c
	$(CC) $(CCFLAGS) $(AMDFLAGS) ../prettify.c -o prettify_amd64

prettify_arm64: ../prettify.c
	$(CC) $(CCFLAGS) $(ARMFLAGS) ../prettify.c -o prettify_arm64

prettify: prettify_amd64 prettify_arm64
	lipo -create -output prettify prettify_amd64 prettify_arm64

clean:
	rm -f plink plink_amd64 plink_arm64 prettify prettify_amd64 prettify_arm64
