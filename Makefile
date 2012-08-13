CFLAGS=-Wall -O2 -lm -lpthread
ZLIB=zlib-1.2.7/libz.so.1.2.7

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
FFLAGS=-framework Accelerate
ZLIB=zlib-1.2.7/libz.1.2.7.dylib
endif

wdist: wdist.c
	gcc $(CFLAGS) $(FFLAGS) wdist.c -o wdist -L. $(ZLIB)

wdist64: wdist.c
	gcc $(CFLAGS) $(FFLAGS) -arch x86_64 wdist.c -o wdist -L. $(ZLIB)
