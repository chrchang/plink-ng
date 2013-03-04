#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>

FILE* g_inf;

size_t readn(unsigned char* readbuf) {
  return fread(readbuf, 1, 131072, g_inf);
}

void parallel_compress(char* out_fname, size_t(* readn)(unsigned char*));

void pigz_defaults(int setprocs);

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
  pigz_defaults(ii);
  g_inf = fopen(argv[1], "rb");
  if (!g_inf) {
    printf("Error: Failed to open %s.\n", argv[1]);
    return 2;
  }
  out_fname = argv[2];
  parallel_compress(out_fname, readn);
  return 0;
}
