BASEFLAGS=-DUSE_MKL
CCFLAGS=-O2
CXXFLAGS=-O2 -I/opt/intel/oneapi/mkl/2025.1/include

LINKFLAGS=-lm -lpthread


include ../Makefile.src
ZLIB64="../../zlib-"$(ZLIB_VER)"/libz-64.a"
CCFLAGS += ${BASEFLAGS} ${CWARN}
CXXFLAGS += ${BASEFLAGS} ${CXXWARN} -I../../2.0/simde "-I../../"$(ZLIB_VER)

plink: $(CSRC2) $(CCSRC2)
	g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	g++ -c $(CXXFLAGS) $(CCSRC2) -c
	g++ -O2 $(OBJ2) -o plink -Wl,--start-group -L. /opt/intel/oneapi/mkl/2025.1/lib/intel64/libmkl_intel_lp64.a -L. /opt/intel/oneapi/mkl/2025.1/lib/intel64/libmkl_intel_thread.a -L. /opt/intel/oneapi/mkl/2025.1/lib/intel64/libmkl_core.a -Wl,--end-group -Wl,-Bstatic -L. /opt/intel/oneapi/compiler/2025.1/lib/libiomp5.a -ldl -static-libgcc -static $(LINKFLAGS) -L. $(ZLIB64)

prettify: ../prettify.c
	gcc $(CCFLAGS) ../prettify.c -c
	gcc -O2 prettify.o -o prettify

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
	rm -f prettify
