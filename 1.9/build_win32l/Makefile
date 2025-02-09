OPENBLAS_ROOT = ../../../openblas-32

BASEFLAGS=
# ***** end configuration *****

BASEFLAGS += -fno-exceptions
BLASFLAGS=-Wl,-Bstatic -L. ${OPENBLAS_ROOT}/lib/libopenblas.a

BASEFLAGS += -I${OPENBLAS_ROOT}/include

include ../Makefile.src
CXXFLAGS=-O2 -std=gnu++11 -I../../2.0/simde "-I../../zlib-"$(ZLIB_VER)
ZLIB_STATIC = "../../zlib-"$(ZLIB_VER)"/libz-32w.a"
LINKFLAGS=-lm -static-libgcc -L. ${ZLIB_STATIC}
CXXFLAGS += ${BASEFLAGS} ${CXXWARN}

plink: $(CSRC2) $(CCSRC2)
	i686-w64-mingw32-g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	i686-w64-mingw32-g++ $(CXXFLAGS) $(CCSRC2) -c
	i686-w64-mingw32-gfortran $(OBJ2) -o plink $(BLASFLAGS) $(LINKFLAGS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
	rm -f plink.exe
