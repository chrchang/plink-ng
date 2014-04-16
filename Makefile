# This is a bit of a mess.  Work with Makefile.std instead.

CFLAGS=-Wall -O2
BLASFLAGS=-L/usr/lib64/atlas -llapack -lcblas -latlas
BLASFLAGS64=-L/usr/lib64/atlas -llapack -lcblas -latlas
LINKFLAGS=-lm -lpthread
ZLIB=zlib-1.2.8/libz.so.1.2.8
ARCH64=-arch x86_64

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
BLASFLAGS=-framework Accelerate
BLASFLAGS64=-framework Accelerate
LINKFLAGS=
ZLIB=zlib-1.2.8/libz.a
ZLIB64=zlib-1.2.8/libz-64.a
else
ifeq ($(UNAME), MINGW32_NT-6.2)
ARCH64=
BLASFLAGS=-Wl,-Bstatic -L. lapack/liblapack.a -L. lapack/librefblas.a
BLASFLAGS64=-Wl,-Bstatic -L. lapack/liblapack-64.a -L. lapack/librefblas-64.a
LINKFLAGS=-lm -static-libgcc
ZLIB=zlib-1.2.8/libz.a
ZLIB64=zlib-1.2.8/libz-64.a
endif
endif

SRC = plink.c plink_assoc.c plink_calc.c plink_cluster.c plink_cnv.c plink_common.c plink_data.c plink_family.c plink_filter.c plink_glm.c plink_help.c plink_homozyg.c plink_lasso.c plink_ld.c plink_matrix.c plink_misc.c plink_set.c plink_stats.c SFMT.c dcdflib.c pigz.c yarn.c

OBJ = $(SRC:.c=.o)

%.o: %.c
	g++ -c $(CFLAGS) $(ARCH64) -o $@ $<

plink: $(SRC)
	g++ $(CFLAGS) $(SRC) -o plink $(BLASFLAGS) $(LINKFLAGS) -L. $(ZLIB)

plinkw: $(SRC)
	g++ $(CFLAGS) $(SRC) -c
	gfortran -O2 $(OBJ) -o plink $(BLASFLAGS) $(LINKFLAGS) -L. $(ZLIB)

plinkc: $(SRC)
	gcc $(CFLAGS) $(SRC) -o plink $(BLASFLAGS64) $(LINKFLAGS) -L. $(ZLIB64)

plinks: $(SRC)
	g++ $(CFLAGS) $(SRC) -o plink_linux_s -Wl,-Bstatic $(BLASFLAGS) -Wl,-Bdynamic $(LINKFLAGS) -L. $(ZLIB)

plinkd: $(SRC)
	g++ $(CFLAGS) $(SRC) -o plink_linux $(BLASFLAGS) -Wl,-Bdynamic $(LINKFLAGS) -L. $(ZLIB)

plinknl: $(SRC)
	g++ $(CFLAGS) $(SRC) -o plink $(LINKFLAGS) -Wl,-Bstatic -L. $(ZLIB)

plink64: $(OBJ)
	g++ $(CFLAGS) $(ARCH64) $(OBJ) -o plink $(BLASFLAGS64) $(LINKFLAGS) -L. $(ZLIB64)
# for clean build, "make clobber" first

plink64w: $(SRC)
	g++ $(CFLAGS) $(ARCH64) $(SRC) -c
	gfortran -O2 $(OBJ) -o plink64 $(BLASFLAGS64) $(LINKFLAGS) -L. $(ZLIB64)

plink64c: $(SRC)
	gcc $(CFLAGS) $(ARCH64) $(SRC) -o plink $(BLASFLAGS64) $(LINKFLAGS) -L. $(ZLIB64)

plink64nl: $(SRC)
	g++ $(CFLAGS) $(ARCH64) $(SRC) -o plink $(LINKFLAGS) -L. $(ZLIB64)

pigz_test: pigz_test.c pigz.c yarn.c
	g++ -Wall -arch x86_64 -O2 pigz_test.c pigz.c yarn.c -o pigz_test -L. $(ZLIB64)

prettify: prettify.c
	gcc -Wall -arch x86_64 -O2 prettify.c -o prettify

prettifyw: prettify.c
	gcc -Wall -O2 prettify.c -o prettify

dbl2txt: dbl2txt.c
	gcc -Wall -arch x86_64 -O2 dbl2txt.c -o dbl2txt

clobber:
	rm -f *.o
