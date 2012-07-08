#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void help_msg() {
  printf("Generates .map and .ped files to stress-test WDIST.\n\n");
  printf("gen_test [number of people] [number of markers]\n");
}

int main(int argc, char* argv[]) {
  int PEOPLE;
  int MARKERS;
  if (argc != 3) {
    help_msg();
    return 1;
  }
  PEOPLE = atoi(argv[1]);
  MARKERS = atoi(argv[2]);
  if ((PEOPLE < 2) || (PEOPLE > 22528) || (MARKERS < 1)) {
    help_msg();
    return 1;
  }
  srand(time(NULL));
  FILE* mapfile = fopen("wdist.map", "w");
  FILE* pedfile = fopen("wdist.ped", "w");
  int ii;
  int jj;
  int kk;
  for (ii = 0; ii < MARKERS; ii += 1) {
    fprintf(mapfile, "1\trs%d\t0\t1000\n", ii);
  }
  fclose(mapfile);
  for (ii = 0; ii < PEOPLE; ii += 1) {
    fprintf(pedfile, "1 1 0 0 1 1");
    for (jj = 0; jj < (MARKERS * 2); jj += 1) {
      kk = (rand() % 2) + 1;
      fprintf(pedfile, " %d", kk);
    }
    fprintf(pedfile, "\n");
  }
  fclose(pedfile);
  return 0;
}
