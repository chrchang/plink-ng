CFLAGS=-Wall -O2
BLASFLAGS=-L/usr/lib64/atlas -llapack -lcblas -latlas
LINKFLAGS=-lm -lpthread
ZLIB=zlib-1.2.7/libz.so.1.2.7

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
BLASFLAGS=-framework Accelerate
LINKFLAGS=
ZLIB=zlib-1.2.7/libz.1.2.7.dylib
endif

wdist: wdist.c
	g++ $(CFLAGS) wdist.c -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. $(ZLIB)

wdistc: wdist.c
	gcc $(CFLAGS) wdist.c -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. $(ZLIB)

wdists: wdist.c
	g++ $(CFLAGS) wdist.c -o wdist_linux_s -Wl,-Bstatic $(BLASFLAGS) -Wl,-Bdynamic $(LINKFLAGS) -L. $(ZLIB)

wdistd: wdist.c
	g++ $(CFLAGS) wdist.c -o wdist_linux $(BLASFLAGS) -Wl,-Bdynamic $(LINKFLAGS) -L. $(ZLIB)

wdist64: wdist.c
	g++ $(CFLAGS) -arch x86_64 wdist.c -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. zlib-1.2.7/libz.1.2.7-64.dylib

wdist64c: wdist.c
	gcc $(CFLAGS) -arch x86_64 wdist.c -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. zlib-1.2.7/libz.1.2.7-64.dylib
