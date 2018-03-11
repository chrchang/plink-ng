BASEFLAGS=-Wall -Wextra -Wshadow -mavx2 -mbmi -mbmi2 -mlzcnt -DNDEBUG -DSTATIC_ZLIB -DZSTD_MULTITHREAD
# BASEFLAGS=-Wall -Wextra -Wshadow -mavx2 -mbmi -mbmi2 -mlzcnt -DSTATIC_ZLIB -DZSTD_MULTITHREAD
# BASEFLAGS=-Wall -Wextra -Wshadow -msse4.2 -DNDEBUG -DSTATIC_ZLIB -DZSTD_MULTITHREAD
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

include Makefile.src
CFLAGS += ${CINCLUDE} ${ZSTD_INCLUDE}

%.o: %.c
	gcc -c $(CFLAGS) $(ARCH32) -o $@ $<
%.o: %.cc
	g++ -c $(CXXFLAGS) $(ARCH32) -o $@ $<

all: plink2 pgen_compress

# for clean build, "make clean" first
plink2: $(OBJ)
	g++ $(ARCH32) $(OBJ) -o build/plink2 $(BLASFLAGS64) $(LINKFLAGS)

# basic pgenlib_internal.h usage example; also needed for tests
pgen_compress: plink2_base.o pgenlib_internal.o pgen_compress.o
	g++ plink2_base.o pgenlib_internal.o pgen_compress.o -o build/pgen_compress

.PHONY: clean
clean:
	rm -f $(CLEAN)
