// Deterministic benchmark/test data generator for ptxt.
//   gen vcf NSAMPLE NVAR   VCF-like file with GT-only genotype columns
//   gen ss NLINE           GWAS summary-statistics TSV (with header)
// Random values are drawn into locals in a fixed order, so the output is the
// same with every compiler.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long long g_state = 88172645463325252ULL;

// xorshift64
static unsigned long long Rnd(void) {
  g_state ^= g_state << 13;
  g_state ^= g_state >> 7;
  g_state ^= g_state << 17;
  return g_state;
}

int main(int argc, char** argv) {
  static char buf[1 << 20];
  const char* nt = "ACGT";
  if (argc < 3) {
    fprintf(stderr, "usage: gen vcf NSAMPLE NVAR | gen ss NLINE\n");
    return 1;
  }
  setvbuf(stdout, buf, _IOFBF, sizeof buf);
  if (!strcmp(argv[1], "vcf")) {
    if (argc < 4) {
      fprintf(stderr, "usage: gen vcf NSAMPLE NVAR\n");
      return 1;
    }
    const int ns = atoi(argv[2]);
    const int nv = atoi(argv[3]);
    printf("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (int i = 0; i < ns; i++) printf("\tS%d", i);
    putchar('\n');
    char* line = malloc(ns * 4 + 1);
    for (int v = 0; v < nv; v++) {
      const int chr = 1 + v * 22 / nv;
      const char ref = nt[Rnd() & 3];
      const char alt = nt[Rnd() & 3];
      const int af_info = (int)(Rnd() % 1000);
      const int r2_info = (int)(Rnd() % 100);
      printf("chr%d\t%d\trs%d\t%c\t%c\t.\tPASS\tAF=0.%03d;INFO=0.9%02d\tGT", chr, 1000 + v * 13, v, ref, alt, af_info, r2_info);
      const unsigned long long af = Rnd() % 64;
      char* p = line;
      for (int i = 0; i < ns; i++) {
        *p++ = '\t';
        *p++ = ((Rnd() & 63) < af) ? '1' : '0';
        *p++ = '|';
        *p++ = ((Rnd() & 63) < af) ? '1' : '0';
      }
      *p = 0;
      fputs(line, stdout);
      putchar('\n');
    }
    free(line);
  } else {
    const long n = atol(argv[2]);
    printf("CHR\tPOS\tSNP\tA1\tA2\tA1_FREQ\tN\tBETA\tSE\tT_STAT\tP\tINFO\n");
    for (long i = 0; i < n; i++) {
      double pval = (double)(Rnd() % 1000000000) / 1e9;
      if ((Rnd() & 1023) == 0) pval *= 1e-8;
      const char a1 = nt[Rnd() & 3];
      const char a2 = nt[Rnd() & 3];
      const double freq = (Rnd() % 10000) / 1e4;
      const int nobs = 400000 + (int)(Rnd() % 1000);
      const double beta = ((double)(Rnd() % 20000) - 10000) / 1e5;
      const double se = (Rnd() % 1000) / 1e5 + 0.001;
      const double tstat = ((double)(Rnd() % 20000) - 10000) / 1e3;
      const double info = (Rnd() % 1000) / 1e3;
      printf("%ld\t%ld\trs%ld\t%c\t%c\t%.4f\t%d\t%.5f\t%.5f\t%.4f\t%.4g\t%.3f\n", 1 + i * 22 / n, 1000 + i * 7, i, a1, a2, freq, nobs, beta, se, tstat, pval, info);
    }
  }
  return 0;
}
