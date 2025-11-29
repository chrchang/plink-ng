# 32-bit Linux build on Debian 10.

BASEFLAGS=
CXEXTRA ?= -DSTABLE_BUILD

BASEFLAGS += ${CXEXTRA}
CCFLAGS=-O2
CXXFLAGS=-O2 -std=c++0x

ARCH32=-m32 -march=i686
BASEFLAGS += ${ARCH32}

LINKFLAGS=-lm -lpthread


include ../Makefile.src
ZLIB32="../../zlib-"$(ZLIB_VER)"/libz-32.a"
CCFLAGS += ${BASEFLAGS} ${CWARN}
CXXFLAGS += ${BASEFLAGS} ${CXXWARN} -Wno-sign-compare -I../../2.0/simde "-I../../zlib-"$(ZLIB_VER)

all: plink prettify

plink: $(CSRC2) $(CCSRC2)
	g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	g++ $(CXXFLAGS) $(CCSRC2) -c
	g++ -O2 $(OBJ2) $(ARCH32) -o plink -Wl,-Bstatic -L/usr/lib -llapack -lblas -lcblas -latlas -lgfortran -lquadmath -ldl -static-libgcc -static $(LINKFLAGS) -L. $(ZLIB32)

prettify: ../prettify.c
	gcc $(CCFLAGS) ../prettify.c -c
	gcc -O2 prettify.o $(ARCH32) -o prettify

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
	rm -f prettify
