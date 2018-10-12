# 32-bit Linux build on Scientific Linux 6.4.

BASEFLAGS=
CXXFLAGS=-O2 -std=c++0x

ARCH32=-m32 -march=i686
BASEFLAGS += ${ARCH32}

ZLIB32=../../zlib-1.2.11/libz-32.a
LINKFLAGS=-lm -lpthread


include ../Makefile.src
CXXFLAGS += ${BASEFLAGS} ${CXXWARN2} -Wno-sign-compare

plink: $(CSRC2) $(CCSRC2)
	g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	g++ $(CXXFLAGS) $(CCSRC2) -c
	g++ -O2 $(OBJ2) $(ARCH32) -o plink -Wl,-Bstatic -L/usr/lib -llapack -lblas -lcblas -latlas -lgfortran -ldl -static-libgcc -static $(LINKFLAGS) -L. $(ZLIB32)

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
