# 32-bit Linux build on Scientific Linux 6.4.

BASEFLAGS=
CXXFLAGS=-O2 -std=c++0x

ARCH32=-m32 -march=i686
BASEFLAGS += ${ARCH32}

LINKFLAGS=-lm -lpthread


include ../Makefile.src
ZLIB32="../../zlib-"$(ZLIB_VER)"/libz-32.a"
CXXFLAGS += ${BASEFLAGS} ${CXXWARN2} -Wno-sign-compare -I../../2.0/simde "-I../../zlib-"$(ZLIB_VER)

plink: $(CSRC2) $(CCSRC2)
	g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	g++ $(CXXFLAGS) $(CCSRC2) -c
	g++ -O2 $(OBJ2) $(ARCH32) -o plink -Wl,-Bstatic -L/usr/lib -llapack -lblas -lcblas -latlas -lgfortran -ldl -static-libgcc -static $(LINKFLAGS) -L. $(ZLIB32)

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
