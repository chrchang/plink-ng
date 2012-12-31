CFLAGS=-Wall -O2
BLASFLAGS=-L/usr/lib64/atlas -llapack -lcblas -latlas
LINKFLAGS=-lm -lpthread
ZLIB=zlib-1.2.7/libz.so.1.2.7

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
BLASFLAGS=-framework Accelerate
LINKFLAGS=
ZLIB=zlib-1.2.7/libz.a
endif

SRC = wdist.c wdist_common.c wdist_data.c wdist_dosage.c

wdist: $(SRC)
	g++ $(CFLAGS) $(SRC) -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. $(ZLIB)

wdistc: $(SRC)
	gcc $(CFLAGS) $(SRC) -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. $(ZLIB)

wdists: $(SRC)
	g++ $(CFLAGS) $(SRC) -o wdist_linux_s -Wl,-Bstatic $(BLASFLAGS) -Wl,-Bdynamic $(LINKFLAGS) -L. $(ZLIB)

wdistd: $(SRC)
	g++ $(CFLAGS) $(SRC) -o wdist_linux $(BLASFLAGS) -Wl,-Bdynamic $(LINKFLAGS) -L. $(ZLIB)

wdist64: $(SRC)
	g++ $(CFLAGS) -arch x86_64 $(SRC) -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. zlib-1.2.7/libz-64.a

wdist64c: $(SRC)
	gcc $(CFLAGS) -arch x86_64 $(SRC) -o wdist $(BLASFLAGS) $(LINKFLAGS) -L. zlib-1.2.7/libz-64.a

wdist64nl: $(SRC)
	g++ $(CFLAGS) -arch x86_64 $(SRC) -o wdist $(LINKFLAGS) -L. zlib-1.2.7/libz-64.a
