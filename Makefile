CFLAGS=-Wall -O2 -lm -lpthread
BLASFLAGS=-lcblas -llapack -L/usr/lib64/atlas
ZLIB=zlib-1.2.7/libz.so.1.2.7

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
BLASFLAGS=-framework Accelerate
ZLIB=zlib-1.2.7/libz.1.2.7.dylib
endif

wdist: wdist.c
	gcc $(CFLAGS) $(BLASFLAGS) wdist.c -o wdist -L. $(ZLIB)

wdist64: wdist.c
	gcc $(CFLAGS) $(BLASFLAGS) -arch x86_64 wdist.c -o wdist -L. $(ZLIB)
