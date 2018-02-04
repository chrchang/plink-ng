BASEFLAGS=-Wall -Wextra -mavx2 -mbmi -mbmi2 -mlzcnt -DSTATIC_ZLIB -DZSTD_MULTITHREAD
# BASEFLAGS=-Wall -Wextra -mavx2 -mbmi -mbmi2 -mlzcnt -DNDEBUG -DSTATIC_ZLIB -DZSTD_MULTITHREAD
CFLAGS=-O3 -std=gnu99
# this actually needs to be named "CXXFLAGS"
CXXFLAGS=-O2 -std=c++14 -Wold-style-cast

CFLAGS += ${BASEFLAGS}
CXXFLAGS += ${BASEFLAGS}

ZLIB=../zlib-1.2.11/libz-64.a
BLASFLAGS64=-llapack -lf77blas -latlas
LINKFLAGS=-lm -lpthread -L. ${ZLIB}
ARCH32=

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  BLASFLAGS=-framework Accelerate
  BLASFLAGS64=-framework Accelerate
  LINKFLAGS=-ldl -lpthread -L. ${ZLIB}
endif

ZSTD_INCLUDE = -Izstd/lib -Izstd/lib/common -Izstd/zlibWrapper
CFLAGS += ${ZSTD_INCLUDE}

CSRC = SFMT.c hfile.c bgzf.c zstd/zlibWrapper/zstd_zlibwrapper.c zstd/zlibWrapper/gzclose.c zstd/zlibWrapper/gzlib.c zstd/zlibWrapper/gzread.c zstd/zlibWrapper/gzwrite.c zstd/lib/common/entropy_common.c zstd/lib/common/zstd_common.c zstd/lib/common/error_private.c zstd/lib/common/xxhash.c zstd/lib/common/fse_decompress.c zstd/lib/common/pool.c zstd/lib/compress/fse_compress.c zstd/lib/compress/huf_compress.c zstd/lib/compress/zstd_double_fast.c zstd/lib/compress/zstd_fast.c zstd/lib/compress/zstd_lazy.c zstd/lib/compress/zstd_ldm.c zstd/lib/compress/zstd_opt.c zstd/lib/compress/zstd_compress.c zstd/lib/compress/zstdmt_compress.c zstd/lib/decompress/huf_decompress.c zstd/lib/decompress/zstd_decompress.c
CPPSRC = plink2_base.cpp pgenlib_internal.cpp plink2.cpp plink2_adjust.cpp plink2_cmdline.cpp plink2_common.cpp plink2_compress_stream.cpp plink2_data.cpp plink2_decompress.cpp plink2_export.cpp plink2_filter.cpp plink2_glm.cpp plink2_help.cpp plink2_import.cpp plink2_ld.cpp plink2_matrix.cpp plink2_matrix_calc.cpp plink2_misc.cpp plink2_psam.cpp plink2_pvar.cpp plink2_random.cpp plink2_set.cpp plink2_stats.cpp

COBJ = $(CSRC:.c=.o)
CPPOBJ = $(CPPSRC:.cpp=.o)

%.o: %.c
	gcc -c $(CFLAGS) $(ARCH32) -o $@ $<
%.o: %.cpp
	g++ -c $(CXXFLAGS) $(ARCH32) -o $@ $<

all: plink2 pgen_compress

# for clean build, "make clean" first
plink2: $(COBJ) $(CPPOBJ)
	g++ $(ARCH32) $(COBJ) $(CPPOBJ) -o build/plink2 $(BLASFLAGS64) $(LINKFLAGS)

# basic pgenlib_internal.h usage example; also needed for tests
pgen_compress: plink2_base.o pgenlib_internal.o pgen_compress.o
	g++ plink2_base.o pgenlib_internal.o pgen_compress.o -o build/pgen_compress

.PHONY: clean
clean:
	rm -f *.o
	rm -f zstd/lib/common/*.o
	rm -f zstd/lib/compress/*.o
	rm -f zstd/lib/decompress/*.o
	rm -f zstd/zlibWrapper/*.o
	rm -f build/plink2
	rm -f build/pgen_compress
