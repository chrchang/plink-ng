# Builds and runs fmtcheck, which checks ptxt's number formatter (FormatNum()
# in 2.0/utils/ptxt/ptxt.cc) against Ryu and printf.
#
#   make          build fmtcheck
#   make check    build and run it (N=... sets values per pattern/precision)

CC ?= cc
CXX ?= c++
ROOT = ../../..
PTXT_CC = $(ROOT)/utils/ptxt/ptxt.cc
N ?= 3000000

all: fmtcheck

# Everything from LayoutG() through the end of FormatNum().
fmt_extracted.h: $(PTXT_CC)
	awk '/^\/\/ Lays out a decimal significand/ {p = 1} p {print} p && /^static uint32_t FormatNum\(double d/ {f = 1} f && /^}$$/ {exit}' $(PTXT_CC) > $@
	grep -q '^static uint32_t FormatNum(double d' $@ || { rm -f $@; echo "fmt_extracted.h: markers not found in ptxt.cc" >&2; exit 1; }

fmtcheck: fmtcheck.cc fmt_extracted.h
	$(CC) -O2 -std=gnu99 -I$(ROOT) -c $(ROOT)/ryu/d2s.c $(ROOT)/ryu/d2fixed.c
	$(CXX) -std=c++17 -O2 -I$(ROOT)/include -I$(ROOT) -I. fmtcheck.cc $(ROOT)/include/plink2_base.cc d2s.o d2fixed.o -o $@

.PHONY: check clean
check: fmtcheck
	./fmtcheck $(N)

clean:
	rm -f fmtcheck fmt_extracted.h *.o
