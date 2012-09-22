#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


double drand()   /* just a uniform distribution for zscore, (0..1] */
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

double zscore()  /* normal distribution, mean 0, sd 1 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

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
  buf = (char*)malloc((MARKERS * 4 + 30) * sizeof(char));
  //  memset(buf, 32, MARKERS * 4 + 20);
  // buf[MARKERS * 4 + 20] = '\n';
  srand(time(NULL));
  FILE* mapfile = fopen("wdist.map", "w");
  FILE* pedfile = fopen("wdist.ped", "wb");
  int ii;
  int jj;
  int kk;
  for (ii = 0; ii < MARKERS; ii += 1) {
    fprintf(mapfile, "1\trs%d\t0\t1000\n", ii);
  }
  fclose(mapfile);
  for (ii = 0; ii < PEOPLE; ii += 1) {
    sprintf(buf, "1 %d 0 0 1 %f ", ii + 1000000000, zscore());
    for (jj = 0; jj < (MARKERS * 2); jj += 1) {
      sprintf(&buf[strlen(buf)], "%d ", rand() %2 + 1); 
    }
    fprintf(pedfile, "%s \n", buf);
  }
  fclose(pedfile);
  return 0;
}
