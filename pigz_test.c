#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>
#include "pigz.h"

FILE* g_inf;

uint32_t emitn(uint32_t overflow_ct, unsigned char* readbuf) {
  return fread(readbuf, 1, 131072, g_inf);
}

int main(int argc, char** argv) {
  int32_t ii;
  char* out_fname;
  if (argc != 3) {
    printf("Usage: pigz_test [input filename] [output filename]\nCompresses a single file using Mark Adler's parallel gzip implementation.\n");
    return 1;
  }
  ii = sysconf(_SC_NPROCESSORS_ONLN);
  if (ii < 1) {
    ii = 1;
  } else if (ii > 8) {
    ii--;
  }
  pigz_init(ii);
  g_inf = fopen(argv[1], "rb");
  if (!g_inf) {
    printf("Error: Failed to open %s.\n", argv[1]);
    return 2;
  }
  out_fname = argv[2];
  parallel_compress(out_fname, emitn);
  return 0;
}
