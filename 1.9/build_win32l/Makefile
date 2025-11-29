OPENBLAS_ROOT = ../../../openblas-32

CXEXTRA ?= -DSTABLE_BUILD
BASEFLAGS=
# ***** end configuration *****

BASEFLAGS += ${CXEXTRA} -fno-exceptions
BLASFLAGS=-Wl,-Bstatic -L. ${OPENBLAS_ROOT}/lib/libopenblas.a

BASEFLAGS += -I${OPENBLAS_ROOT}/include

include ../Makefile.src
CCFLAGS=-O2
CXXFLAGS=-O2 -std=gnu++11 -I../../2.0/simde "-I../../zlib-"$(ZLIB_VER)
ZLIB_STATIC = "../../zlib-"$(ZLIB_VER)"/libz-32w.a"
LINKFLAGS=-lm -static-libgcc -L. ${ZLIB_STATIC}
CCFLAGS += ${BASEFLAGS} ${CWARN}
CXXFLAGS += ${BASEFLAGS} ${CXXWARN}

all: plink prettify

plink: $(CSRC2) $(CCSRC2)
	i686-w64-mingw32-g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	i686-w64-mingw32-g++ $(CXXFLAGS) $(CCSRC2) -c
	i686-w64-mingw32-gfortran -O2 $(OBJ2) -o plink $(BLASFLAGS) $(LINKFLAGS)

prettify: ../prettify.c
	i686-w64-mingw32-gcc $(CCFLAGS) ../prettify.c -c
	i686-w64-mingw32-gcc -O2 prettify.o -o prettify

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
	rm -f plink.exe
	rm -f prettify
	rm -f prettify.exe
