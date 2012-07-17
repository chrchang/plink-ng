#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void help_msg() {
  printf("Generates .map and .ped files to stress-test WDIST.\n\n");
  printf("gen_test [number of people] [number of markers]\n");
}

int main(int argc, char* argv[]) {
  int PEOPLE;
  int MARKERS;
  char* buf;
  if (argc != 3) {
    help_msg();
    return 1;
  }
  PEOPLE = atoi(argv[1]);
  MARKERS = atoi(argv[2]);
  if ((PEOPLE < 2) || (MARKERS < 1)) {
    help_msg();
    return 1;
  }
  buf = (char*)malloc((MARKERS * 4 + 21) * sizeof(char));
  memset(buf, 32, MARKERS * 4 + 20);
  buf[MARKERS * 4 + 20] = '\n';
  srand(time(NULL));
  FILE* mapfile = fopen("wdist.map", "w");
  FILE* pedfile = fopen("wdist.ped", "wb");
  int ii;
  int jj;
  for (ii = 0; ii < MARKERS; ii += 1) {
    fprintf(mapfile, "1\trs%d\t0\t1000\n", ii);
  }
  fclose(mapfile);
  for (ii = 0; ii < PEOPLE; ii += 1) {
    sprintf(buf, "1 %d 0 0 1 1 ", ii + 1000000000);
    for (jj = 0; jj < (MARKERS * 2); jj += 1) {
      buf[jj * 2 + 21] = (rand() % 2) + '1';
    }
    fwrite(buf, 1, MARKERS * 4 + 21, pedfile);
  }
  fclose(pedfile);
  return 0;
}
