OPENBLAS_ROOT = ../../../openblas-64
ZLIB_STATIC = ../../zlib-1.3.1/libz-64w.a

BASEFLAGS=
# ***** end configuration *****

BASEFLAGS += -fno-exceptions
LINKFLAGS=-lm -static-libgcc -L. ${ZLIB_STATIC}
BLASFLAGS=-Wl,-Bstatic -L. ${OPENBLAS_ROOT}/lib/libopenblas.a

BASEFLAGS += -I${OPENBLAS_ROOT}/include
CXXFLAGS=-O2 -std=gnu++11

include ../Makefile.src
CXXFLAGS += ${BASEFLAGS} ${CXXWARN} -I../../2.0/simde

plink: $(CSRC2) $(CCSRC2)
	x86_64-w64-mingw32-g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	x86_64-w64-mingw32-g++ -c $(CXXFLAGS) $(CCSRC2) -c
	x86_64-w64-mingw32-gfortran $(OBJ2) -o plink $(BLASFLAGS) $(LINKFLAGS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
	rm -f plink.exe
