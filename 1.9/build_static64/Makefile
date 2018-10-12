BASEFLAGS=-DUSE_MKL
CXXFLAGS=-O2 -I/opt/intel/mkl/include

ZLIB64=../../zlib-1.2.11/libz-64.a
LINKFLAGS=-lm -lpthread


include ../Makefile.src
CXXFLAGS += ${BASEFLAGS} ${CXXWARN}

plink: $(CSRC2) $(CCSRC2)
	g++ -x c++ $(CXXFLAGS) $(CSRC2) -c
	g++ -c $(CXXFLAGS) $(CCSRC2) -c
	g++ $(OBJ2) -o plink -Wl,--start-group -L. /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a -L. /opt/intel/mkl/lib/intel64/libmkl_intel_thread.a -L. /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -Wl,-Bstatic -L. /opt/intel/lib/intel64/libiomp5.a -ldl -static-libgcc -static $(LINKFLAGS) -L. $(ZLIB64)

.PHONY: clean
clean:
	rm -f *.o
	rm -f plink
