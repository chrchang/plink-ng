# MinGW/MinGW-w64 Makefile for PLINK 2.00.
#
# Compilation options (leave blank after "=" to disable, put "= 1" to enable):
#   Do not use SSE 4.2 instructions: NO_SSE42
#   Do not link to OpenBLAS: NO_OPENBLAS
#   Use only -O2 optimization for zstd: ZSTD_O2
NO_SSE42 = 1
NO_OPENBLAS =
ZSTD_O2 = 1

OPENBLAS_ROOT = ../../openblas
ZLIB_STATIC = ../../zlib-1.2.11/libz.a

BASEFLAGS=-Wall -Wextra
# ***** end configuration *****

BASEFLAGS += -DSTATIC_ZLIB -fno-exceptions
LINKFLAGS=-lm -static-libgcc -L. ${ZLIB_STATIC}
ifdef NO_OPENBLAS
  BASEFLAGS += -DNOLAPACK
  BLASFLAGS=
else
  BASEFLAGS += -DUSE_OPENBLAS
  BLASFLAGS=-Wl,-Bstatic -L. ${OPENBLAS_ROOT}/lib/libopenblas.a
endif

CFLAGS=-O2 -std=gnu99
CXXFLAGS=-O2 -std=gnu++11

ifndef NO_SSE42
  BASEFLAGS += -msse4.2
endif

ifdef ZSTD_O2
  ZCFLAGS=-O2 -std=gnu99
else
  ZCFLAGS=-O3 -std=gnu99
endif

BASEFLAGS += -I${OPENBLAS_ROOT}/include

CFLAGS += ${BASEFLAGS}
ZCFLAGS += ${BASEFLAGS}
CXXFLAGS += ${BASEFLAGS}

ZSTD_INCLUDE = -I../zstd/lib -I../zstd/lib/common -I../zstd/zlibWrapper
ZCFLAGS += ${ZSTD_INCLUDE}

CSRC = ../SFMT.c ../hfile.c ../bgzf.c
ZCSRC = ../zstd/zlibWrapper/zstd_zlibwrapper.c ../zstd/zlibWrapper/gzclose.c ../zstd/zlibWrapper/gzlib.c ../zstd/zlibWrapper/gzread.c ../zstd/zlibWrapper/gzwrite.c ../zstd/lib/common/entropy_common.c ../zstd/lib/common/zstd_common.c ../zstd/lib/common/error_private.c ../zstd/lib/common/xxhash.c ../zstd/lib/common/fse_decompress.c ../zstd/lib/compress/fse_compress.c ../zstd/lib/compress/huf_compress.c ../zstd/lib/compress/zstd_compress.c ../zstd/lib/decompress/huf_decompress.c ../zstd/lib/decompress/zstd_decompress.c
CPPSRC = ../pgenlib_internal.cpp ../plink2.cpp ../plink2_adjust.cpp ../plink2_common.cpp ../plink2_compress_stream.cpp ../plink2_data.cpp ../plink2_decompress.cpp ../plink2_filter.cpp ../plink2_glm.cpp ../plink2_help.cpp ../plink2_ld.cpp ../plink2_matrix.cpp ../plink2_matrix_calc.cpp ../plink2_misc.cpp ../plink2_psam.cpp ../plink2_pvar.cpp ../plink2_random.cpp ../plink2_set.cpp ../plink2_stats.cpp

OBJ = SFMT.o hfile.o bgzf.o zstd_zlibwrapper.o gzclose.o gzlib.o gzread.o gzwrite.o entropy_common.o zstd_common.o error_private.o xxhash.o fse_decompress.o fse_compress.o huf_compress.o zstd_compress.o huf_decompress.o zstd_decompress.o pgenlib_internal.o plink2.o plink2_adjust.o plink2_common.o plink2_compress_stream.o plink2_data.o plink2_decompress.o plink2_filter.o plink2_glm.o plink2_help.o plink2_ld.o plink2_matrix.o plink2_matrix_calc.o plink2_misc.o plink2_psam.o plink2_pvar.o plink2_random.o plink2_set.o plink2_stats.o

all: plink2 pgen_compress

plink2: $(CSRC) $(ZCSRC) $(CPPSRC)
	gcc $(CFLAGS) $(CSRC) -c
	gcc $(ZCFLAGS) $(ZCSRC) -c
	g++ $(CXXFLAGS) $(CPPSRC) -c
	gfortran $(OBJ) -o plink2 $(BLASFLAGS) $(LINKFLAGS)

pgen_compress: ../pgenlib_internal.cpp ../pgen_compress.cpp
	g++ $(CXXFLAGS) ../pgenlib_internal.cpp ../pgen_compress.cpp -o pgen_compress

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink2.exe
	rm -f pgen_compress.exe
