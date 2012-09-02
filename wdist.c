// WDIST weighted genetic distance calculator
// Copyright (C) 2012  Christopher Chang  chrchang@alumni.caltech.edu

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// uncomment this to build without CBLAS/CLAPACK
// #define NOLAPACK

// TODO:
// distance MAF histograms

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

#if __LP64__
#include <emmintrin.h>
#endif

#ifndef NOLAPACK
#include <cblas.h>
#ifdef __APPLE__
#include <vecLib/clapack.h>
#else
#include <clapack.h>
#endif
#endif

#include "zlib-1.2.7/zlib.h"

#define PI 3.141592653589793
#define CACHELINE 64 // assumed number of bytes per cache line, for alignment
#define CACHELINE_DBL (CACHELINE / sizeof(double))

#define RET_SUCCESS 0
#define RET_HELP 1
#define RET_NOMEM 2
#define RET_OPENFAIL 3
#define RET_INVALID_FORMAT 4
#define RET_CALC_NOT_YET_SUPPORTED 5
#define RET_INVALID_CMDLINE 6
#define RET_WRITE_FAIL 7
#define RET_READ_FAIL 8
#define RET_MOREHELP 9
#define RET_THREAD_CREATE_FAIL 10

#define FILTER_KEEP 1
#define FILTER_REMOVE 2
#define FILTER_CUSTOM 3

#define CALC_RELATIONSHIP_NO_VAR_STD 1
#define CALC_RELATIONSHIP_SQ 2
#define CALC_RELATIONSHIP_SQ0 4
#define CALC_RELATIONSHIP_GZ 8
#define CALC_RELATIONSHIP_BIN 16
#define CALC_RELATIONSHIP_NULL 24
#define CALC_RELATIONSHIP_GRM 32
#define CALC_RELATIONSHIP_MASK 63
#define CALC_IBC 64
#define CALC_DISTANCE_SQ 128
#define CALC_DISTANCE_SQ0 256
#define CALC_DISTANCE_GZ 512
#define CALC_DISTANCE_BIN 1024
#define CALC_DISTANCE_NULL 1536
#define CALC_DISTANCE_GZBIN_MASK 1536
#define CALC_DISTANCE_MASK 1792
#define CALC_GROUPDIST 2048
#define CALC_REGRESS_DISTANCE 4096
#define CALC_UNRELATED_HERITABILITY 8192
#define CALC_UNRELATED_HERITABILITY_STRICT 16384
#define CALC_FREQ 32768
#define CALC_FREQ_GCTA 65536
#define CALC_REL_CUTOFF 131072
#define CALC_WRITE_SNPLIST 262144

#define _FILE_OFFSET_BITS 64
#define MAX_THREADS 63
#define MAX_THREADS_P1 64
#define PEDBUFBASE 256
#define FNAMESIZE 2048
#define MALLOC_DEFAULT_MB 2176
// size of generic text line load buffer.  .ped lines can of course be longer
#define MAXLINELEN 131072

#define MAX_EM_ACCEL 100.0

// default jackknife iterations
#define ITERS_DEFAULT 100000

#define DEFAULT_CHROM_MASK 0x027fffff

// Number of snp-major .bed lines to read at once for distance calc if exponent
// is zero.  Currently assumed to be a multiple of 192, and no larger than
// 1920, by the popcount_..._multibyte functions.  (The optimal value depends
// on both system-specific properties such as cache sizes, as well as the
// number of individuals in the current calculation, so in principle it's best
// to select this value at runtime.  But 960 usually works well in practice in
// my experience.)
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)

#if __LP64__
#define BITCT 64
// number of snp-major .bed lines to read at once for distance calc if exponent
// is nonzero.
#define MULTIPLEX_DIST_EXP 64
// number of snp-major .bed lines to read at once for relationship calc
#define MULTIPLEX_REL 60
#else
#define BITCT 32
#define MULTIPLEX_DIST_EXP 28
#define MULTIPLEX_REL 30
#endif

#define BITCT2 (BITCT / 2)

#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - 1)))
#define CACHEALIGN_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1)))

#define MIN(aa, bb) ((aa) > (bb))? (bb) : (aa)

const char info_str[] =
#ifdef NOLAPACK
  "WDIST weighted genetic distance calculator, v0.5.11 (2 September 2012)\n"
#else
  "WDIST weighted genetic distance calculator, v0.7.5 (2 September 2012)\n"
#endif
  "Christopher Chang (chrchang@alumni.caltech.edu), BGI Cognitive Genomics Lab\n\n"
  "wdist [flags...]\n";
const char errstr_append[] = "\nRun 'wdist --help | more' for more information.\n";
const char errstr_map_format[] = "Error: Improperly formatted .map file.\n";
const char errstr_fam_format[] = "Error: Improperly formatted .fam file.\n";
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
char tbuf[MAXLINELEN];
#ifndef __LP64__
unsigned char popcount[65536];
#endif

// manually managed, very large stack
unsigned char* wkspace;
unsigned char* wkspace_base;
long long malloc_size_mb = MALLOC_DEFAULT_MB;
long long wkspace_left;

unsigned char* wkspace_alloc(long long size) {
  unsigned char* retval;
  if (wkspace_left < size) {
    return NULL;
  }
  size = CACHEALIGN(size);
  retval = wkspace_base;
  wkspace_base += size;
  wkspace_left -= size;
  return retval;
}

void wkspace_reset(unsigned char* new_base) {
  long long freed_bytes = wkspace_base - new_base;
  wkspace_base = new_base;
  wkspace_left += freed_bytes;
}

char** subst_argv = NULL;

int dispmsg(int retval) {
  switch(retval) {
  case RET_HELP:
    printf("%s%s", info_str, errstr_append);
    break;
  case RET_NOMEM:
    printf("Error: Out of memory.  Try the --memory flag.\n");
    break;
  case RET_WRITE_FAIL:
    printf("Error: File write failure.\n");
    break;
  case RET_READ_FAIL:
    printf("Error: File read failure.\n");
    break;
  case RET_MOREHELP:
    printf(
"%s\n"
"In the command line flag definitions that follow,\n"
"  * [square brackets] denote a required parameter, where the text between the\n"
"    brackets describes its nature.\n"
"  * <angle brackets> denote an optional modifier (or if '|' is present, a set\n"
"    of mutually exclusive optional modifiers).  Use the EXACT text in the\n"
"    definition, e.g. '--distance square0'.\n"
"  * {curly braces} denote an optional parameter, where the text between the\n"
"    braces describes its nature.\n\n"
"Each run must invoke at least one of the following calculations:\n\n"
"  --freq <gcta>\n"
"    Generates the same allele frequency report as PLINK --freq (or GCTA --freq,\n"
"    if the 'gcta' modifier is present).\n\n"
"  --distance <square | square0> <gz | bin>\n"
"    Writes a lower-triangular tab-delimited table of (weighted) genetic\n"
"    distances to {output prefix}.dist, and a list of the corresponding\n"
"    family/individual IDs to {output prefix}.dist.id.\n"
"    The first row of the .dist file contains a single number describing the\n"
"    distance between the first two genotypes, the second row has the {genotype\n"
"    1-genotype 3} and {genotype 2-genotype 3} distances in that order, etc.\n"
"    If the 'square' or 'square0' modifier is present, a square matrix is\n"
"    written instead; 'square0' fills the upper right triangle with zeroes.\n"
"    If the 'gz' modifier is present, a compressed .dist.gz file is written\n"
"    instead of a plain text file.\n"
"    If the 'bin' modifier is present, a binary (square) matrix of\n"
"    double-precision floating point values, suitable for loading from R, is\n"
"    instead written to {output prefix}.dist.bin.  This can be combined with\n"
"    'square0' if you still want the upper right zeroed out.\n\n"
"  --ibc\n"
"    This yields the same output as GCTA --ibc.\n\n"
"  --make-rel <no-var-std> <square | square0> <gz | bin> <ibc1 | ibc2>\n"
"    Writes a lower-triangular variance-standardized relationship matrix to\n"
"    {output prefix}.rel, and corresponding IDs to {output prefix}.rel.id.\n"
"    The 'no-var-std' modifier causes GEMMA's \"centered relationship matrix\" to\n"
"    be calculated instead.\n"
"    'square', 'square0', 'gz', and 'bin' act as they do on --distance.\n"
"    By default, the diagonal elements in the variance-standardized matrix are\n"
"    based on --ibc's Fhat3.  Use the 'ibc1' or 'ibc2' modifiers to base them on\n"
"    Fhat1 or Fhat2 instead.\n"
"  --make-grm <no-var-std> <no-gz> <ibc1 | ibc2>\n"
"    Writes the relationship matrix in GCTA's .grm format instead.  Since GCTA\n"
"    normally compresses the .grm file, WDIST does the same unless the 'no-gz'\n"
"    modifier is present.\n"
"  --make-grm-bin <no-var-std> <square0> <ibc1 | ibc2>\n"
"    This is identical to --make-rel bin, except the output file extension is\n"
"    .grm.bin.\n\n"
"  --groupdist {iters} {d}\n"
"    Two-group genetic distance analysis, using delete-d jackknife with the\n"
"    requested number of iterations to estimate standard errors.  Binary\n"
"    phenotype required.\n"
"    If only one parameter is provided, d defaults to {number of people}^0.6\n"
"    rounded down.  With no parameters, 100k iterations are run.\n\n"
"  --regress-distance {iters} {d}\n"
"    Regresses genetic distances on average phenotypes, using delete-d\n"
"    jackknife for standard errors.  Scalar phenotype required.  Defaults for\n"
"    iters and d are the same as for --groupdist.\n\n"
#ifndef NOLAPACK
"  --unrelated-heritability <strict> {tol} {initial covg} {initial cove}\n"
"    REML estimate of additive heritability, iterating with an accelerated\n"
"    variant of the EM algorithm until the rate of change of the log likelihood\n"
"    function is less than tol.  Scalar phenotype required.\n"
"    The 'strict' modifier forces regular EM to be used.  tol defaults to\n"
"    10^{-7}, genetic covariance prior defaults to 0.45, and residual\n"
"    covariance prior defaults to (1 - covg).\n"
"    For more details, see Vattikuti S, Guo J, Chow CC (2012) Heritability and\n"
"    Genetic Correlations Explained by Common SNPs for Metabolic Syndrome\n"
"    Traits.  PLoS Genet 8(3): e1002637.  doi:10.1371/journal.pgen.1002637\n\n"
#endif
"  --write-snplist\n"
"    Write a .snplist file listing the names of all SNPs that pass the filters\n"
"    you've specified.\n\n"
"The following other flags are supported.\n"
"  --script [fname] : Include command-line options from file.\n"
"  --file [prefix]  : Specify prefix for .ped and .map files.  (When this flag\n"
"                     isn't present, the prefix is assumed to be 'wdist'.)\n"
"  --ped [filename] : Specify name of .ped file.\n"
"  --map [filename] : Specify name of .map file.\n"
"  --make-bed       : Make .bed, .bim, and .fam files.  Note: this is ALWAYS ON\n"
"                     for now, since the program core currently only handles\n"
"                     binary files.\n"
"  --no-fid         : .ped file does not contain column 1 (family ID).\n"
"  --no-parents     : .ped file does not contain columns 3-4 (parents).\n"
"  --no-sex         : .ped file does not contain column 5 (sex).\n"
"  --no-pheno       : .ped file does not contain column 6 (phenotype).\n"
"  --liability      : .ped file does contain liability (column 7).\n"
"  --bfile {prefix} : Specify .bed/.bim/.fam prefix (default 'wdist').\n"
"  --bed [filename] : Specify .bed file.\n"
"  --bim [filename] : Specify .bim file.\n"
"  --fam [filename] : Specify .fam file.\n"
"  --out [prefix]   : Specify prefix for output files.\n"
"  --silent         : Suppress output to console.\n"
"  --pheno [fname]  : Specify alternate phenotype.\n"
"  --mpheno [col]   : Specify phenotype column number in --pheno file.\n"
"  --pheno-name [c] : If phenotype file has a header row, use column with the\n"
"                     given name.\n"
"  --prune          : Remove individuals with missing phenotypes.\n"
"  --1              : .ped file affection phenotypes are interpreted as\n"
"                     0 = unaffected, 1 = affected (instead of 0 = missing,\n"
"                     1 = unaffected, 2 = affected)\n"
// --map3 implicitly supported via autodetection
// --compound-genotypes automatically supported
"  --autosome       : Include markers on all autosomes (chromosomes 1-22), and\n"
"                     no others.\n"
"  --chr [num...]   : Only consider markers on the given chromosome(s).  Valid\n"
"                     choices are 0 (unplaced), 1-22, and XY (pseudo-autosomal\n"
"                     region of X).  Separate multiple chromosomes with spaces,\n"
"                     e.g. '--chr 1 3 xy'.\n"
"                     The X, Y, and MT chromosome values are currently\n"
"                     unsupported, but this will change in the future.\n"
"  --maf {val}      : Minor allele frequency minimum threshold (default 0.01).\n"
"                     Note that the default threshold is only applied if --maf\n"
"                     is used without an accompanying value; if you do not\n"
"                     invoke --maf, WDIST does not perform any MAF-based SNP\n"
"                     filtering at all.  The other filters work the same way.\n"
"  --geno {val}     : Maximum per-SNP missing (default 0.1).\n"
"  --mind {val}     : Maximum per-person missing (default 0.1).\n"
"  --hwe {val}      : Minimum Hardy-Weinberg disequilibrium p-value (exact),\n"
"                     default 0.001.  This is checked after all other forms of\n"
"                     filtering except for --rel-cutoff.\n"
"  --rel-cutoff {v} : Exclude individuals until no remaining pairs have\n"
"  --grm-cutoff {v}   relatedness greater than the given cutoff value (default\n"
"                     0.025).  Note that maximizing the remaining sample size is\n"
"                     equivalent to the NP-hard maximum independent set problem,\n"
"                     so we use a greedy algorithm instead of guaranteeing\n"
"                     optimality.  (Use the --make-rel and --keep/--remove flags\n"
"                     if you want to try to do better.)\n"
"  --rseed [val]    : Set random number seed (relevant for jackknife standard\n"
"                     error estimation).\n"
"  --memory [val]   : Size, in MB, of initial malloc attempt.\n"
"  --threads [val]  : Maximum number of concurrent threads.\n"
"  --extract [file] : Only include SNPs in the given list.\n"
"  --exclude [file] : Exclude all SNPs in the given list.\n"
"  --keep [fname]   : Only include individuals in the given list.\n"
"  --remove [fname] : Exclude all individuals in the given list.\n"
"  --exponent [val] : When computing genetic distances, each marker has a weight\n"
"                     of (2q(1-q))^{-val}, where q is the observed MAF after\n"
"                     applying --keep/--remove/--filter/--prune and before\n"
"                     applying any other filters.  (Use --update-freq if you\n"
"                     want to explicitly specify some or all of the MAFs.)\n"
"  --update-freq [filename]  : Loads MAFs from the given PLINK- or GCTA-style\n"
"                              frequency file, instead of just setting them to\n"
"                              frequencies observed in the .ped/.bed file.\n"
"                              (This can be important when performing multiple\n"
"                              rounds of filtering, because MAFs affect the\n"
"                              behavior of some filters, and observed MAFs\n"
"                              change when you filter people out.)\n"
"  --filter [filename] [val] : Filter individuals (see PLINK documentation).\n"
"  --mfilter [col]           : Specify column number in --filter file.\n"
"  --missing-genotype [char] : Code for missing genotype.\n"
"  --missing-phenotype [val] : Code for missing phenotype.\n"
"  --make-pheno [file] [val] : Specify binary phenotype, where cases have the\n"
"                              given value.  If the value is '*', all\n"
"                              individuals present in the phenotype file are\n"
"                              affected (and other individuals in the .ped/.fam\n"
"                              are unaffected).\n"
"  --tail-pheno [Ltop] {Hbt} : Form 'low' (<= Ltop, unaffected) and 'high'\n"
"                              (greater than Hbt, affected) groups from scalar\n"
"                              phenotype data.  If Hbt is unspecified, it is set\n"
"                              equal to Ltop.  Central phenotype values are\n"
"                              treated as missing.\n\n"
         , info_str);
    break;
  }
  if (subst_argv) {
    free(subst_argv);
  }
  return retval;
}

// (copied from PLINK helper.cpp, modified to perform threshold test and avoid
// repeated allocation of the same memory)
//
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton

double* het_probs = NULL;

int SNPHWE_t(int obs_hets, int obs_hom1, int obs_hom2, double thresh) {
  if ((thresh == 0.0) || !(obs_hom1 + obs_hom2 + obs_hets)) {
    return 0;
  }
  
  // if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) 
  //   {
  //     error("Internal error: negative count in HWE test: "
  //           +int2str(obs_hets)+" "
  //           +int2str(obs_hom1)+" "
  //           +int2str(obs_hom2));
  //   }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int rare_copies2 = rare_copies / 2;
  int genotypes   = obs_hets + obs_homc + obs_homr;
  // avoid integer overflow when calculating mid
  long long rare_copies_ll = rare_copies;
  long long genotypes_ll = genotypes;

  int i;
  for (i = 0; i <= rare_copies2; i++) {
    het_probs[i] = 0.0;
  }

  /* start at midpoint */
  int mid = (int)((rare_copies_ll * (2 * genotypes_ll - rare_copies_ll)) / (2 * genotypes_ll));
  
  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies ^ mid) & 1) {
    mid++;
  }
  int mid2 = mid / 2;
  
  long long curr_hets = mid;

  het_probs[mid2] = 1.0;
  double sum = 1.0;
  long long curr_homr = rare_copies2 - mid2;
  int curr_homc = genotypes - mid - curr_homr;
  double* hptr = &(het_probs[mid2 - 1]);
  double dv = 1.0;
  while (curr_hets > 1) {
    /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
    *hptr = dv * (double)(curr_hets * (curr_hets - 1)) / (double)(4 * (++curr_homr) * (++curr_homc));
    curr_hets -= 2;
    dv = *hptr--;
    sum += dv;
  }

  curr_hets = mid;
  curr_homr = rare_copies2 - mid2;
  curr_homc = genotypes - mid - curr_homr;
  hptr = &(het_probs[mid2 + 1]);
  dv = 1.0;
  while (curr_hets <= rare_copies - 2) {
    curr_hets += 2;
    *hptr = dv * (double)(4 * (curr_homr--) * (curr_homc--)) / (double)(curr_hets * (curr_hets - 1));
    dv = *hptr++;
    sum += dv;
  }

  thresh *= sum;

  double p_hwe = 0.0;
  hptr = &(het_probs[rare_copies2]);
  sum = het_probs[obs_hets / 2];
  /*  p-value calculation for p_hwe  */
  while (*hptr <= sum) {
    p_hwe += *hptr--;
    if (p_hwe > thresh) {
      return 0;
    }
  }
  hptr = het_probs;
  while (*hptr <= sum) {
    p_hwe += *hptr++;
  }
  return (!(p_hwe > thresh));
}


// A C-program for MT19937, with initialization improved 2002/1/26.
// Coded by Takuji Nishimura and Makoto Matsumoto.

// Before using, initialize the state by using init_genrand(seed)  
// or init_by_array(init_key, key_length).

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.                          

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:

//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.

//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.

//   3. The names of its contributors may not be used to endorse or promote 
//      products derived from this software without specific prior written 
//      permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// Any feedback is very welcome.
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
// email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

// Period parameters
#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   // constant vector a
#define UPPER_MASK 0x80000000UL // most significant w-r bits
#define LOWER_MASK 0x7fffffffUL // least significant r bits

static unsigned long mt[MT_N]; // the array for the state vector
static int mti=MT_N+1; // mti==N+1 means mt[N] is not initialized

// initializes mt[MT_N] with a seed
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<MT_N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
        // In the previous versions, MSBs of the seed affect
        // only MSBs of the array mt[].
        // 2002/01/09 modified by Makoto Matsumoto
        mt[mti] &= 0xffffffffUL;
        // for >32 bit machines
    }
}

// see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
// if you want to add init_by_array back

// generates a random number on [0,0xffffffff]-interval
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1

    if (mti >= MT_N) { // generate MT_N words at one time
        int kk;

        if (mti == MT_N+1)   // if init_genrand() has not been called,
	  init_genrand(5489UL); // a default initial seed is used

        for (kk=0;kk<MT_N-MT_M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<MT_N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}


// back to our regular program

inline int is_space_or_eoln(char cc) {
  return ((cc == ' ') || (cc == '\t') || (cc == '\n') || (cc == '\0'));
}

int strlen_se(char* ss) {
  int val = 0;
  while (!is_space_or_eoln(*ss++)) {
    val++;
  }
  return val;
}

char* id_buf = NULL;

int int_cmp(const void* aa, const void* bb) {
  return *((const int*)aa) - *((const int*)bb);
}

int double_cmp(const void* aa, const void* bb) {
  double cc = *((const double*)aa) - *((const double*)bb);
  if (cc > 0.0) {
    return 1;
  } else if (cc < 0.0) {
    return -1;
  } else {
    return 0;
  }
}

double get_dmedian_i(int* sorted_arr, int len) {
  if (len) {
    if (len % 2) {
      return (double)sorted_arr[len / 2];
    } else {
      return ((double)sorted_arr[len / 2] + (double)sorted_arr[(len / 2) - 1]) * 0.5;
    }
  } else {
    return 0.0;
  }
}

double get_dmedian(double* sorted_arr, int len) {
  if (len) {
    if (len % 2) {
      return sorted_arr[len / 2];
    } else {
      return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5;
    }
  } else {
    return 0.0;
  }
}

int strcmp_casted(const void* s1, const void* s2) {
  return strcmp((char*)s1, (char*)s2);
}

int strcmp_deref(const void* s1, const void* s2) {
  return strcmp(*(char**)s1, *(char**)s2);
}

// alas, qsort_r not available on some Linux distributions
int qsort_ext(char* main_arr, int arr_length, int item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, int secondary_item_len) {
  char* proxy_arr;
  int proxy_len = secondary_item_len + sizeof(void*);
  int ii;
  if (!arr_length) {
    return 0;
  }
  if (proxy_len < item_length) {
    proxy_len = item_length;
  }
  proxy_arr = (char*)malloc(arr_length * proxy_len);
  if (!proxy_arr) {
    return -1;
  }
  for (ii = 0; ii < arr_length; ii++) {
    *(char**)(&(proxy_arr[ii * proxy_len])) = &(main_arr[ii * item_length]);
    memcpy(&(proxy_arr[ii * proxy_len + sizeof(void*)]), &(secondary_arr[ii * secondary_item_len]), secondary_item_len);
  }
  qsort(proxy_arr, arr_length, proxy_len, comparator_deref);
  for (ii = 0; ii < arr_length; ii++) {
    memcpy(&(secondary_arr[ii * secondary_item_len]), &(proxy_arr[ii * proxy_len + sizeof(void*)]), secondary_item_len);
    memcpy(&(proxy_arr[ii * proxy_len]), *(char**)(&(proxy_arr[ii * proxy_len])), item_length);
  }
  for (ii = 0; ii < arr_length; ii++) {
    memcpy(&(main_arr[ii * item_length]), &(proxy_arr[ii * proxy_len]), item_length);
  }
  free(proxy_arr);
  return 0;
}

int bsearch_str(char* id_buf, char* lptr, int max_id_len, int min_idx, int max_idx) {
  int mid_idx;
  int ii;
  if (max_idx < min_idx) {
    return -1;
  }
  mid_idx = (min_idx + max_idx) / 2;
  ii = strcmp(id_buf, &(lptr[mid_idx * max_id_len]));
  if (ii) {
    if (ii < 0) {
      return bsearch_str(id_buf, lptr, max_id_len, min_idx, mid_idx - 1);
    } else {
      return bsearch_str(id_buf, lptr, max_id_len, mid_idx + 1, max_idx);
    }
  } else {
    return mid_idx;
  }
}

int bsearch_fam_indiv(char* id_buf, char* lptr, int max_id_len, int filter_lines, char* fam_id, char* indiv_id) {
  int ii;
  int jj;
  char* sptr;
  if (!filter_lines) {
    return -1;
  }
  ii = 0;
  jj = 0;
  sptr = fam_id;
  while (!is_space_or_eoln(*sptr)) {
    sptr++;
    ii++;
  }
  sptr = indiv_id;
  while (!is_space_or_eoln(*sptr)) {
    sptr++;
    jj++;
  }
  if (ii + jj + 2 > max_id_len) {
    return -1;
  }
  memcpy(id_buf, fam_id, ii);
  id_buf[ii] = ' ';
  memcpy(&(id_buf[ii + 1]), indiv_id, jj);
  id_buf[ii + jj + 1] = '\0';
  return bsearch_str(id_buf, lptr, max_id_len, 0, filter_lines - 1);
}

inline int is_contained(char* id_buf, char* lptr, int max_id_len, int filter_lines, char* fam_id, char* indiv_id) {
  return (bsearch_fam_indiv(id_buf, lptr, max_id_len, filter_lines, fam_id, indiv_id) != -1);
}

int marker_code(char* sptr) {
  // does not require string to be null-terminated, and does not perform
  // exhaustive error-checking
  if (sptr[1] > ' ') {
    if ((sptr[0] == 'X') || (sptr[0] == 'x')) {
      return 25; // XY
    } else if ((sptr[0] == 'M') || (sptr[0] == 'm')) {
      return 26; // MT
    } else {
      return ((sptr[0] - '0') * 10 + (sptr[1] - '0'));
    }
  } else if (*sptr < 'X') {
    return (sptr[0] - '0');
  } else if ((*sptr == 'X') || (*sptr == 'x')) {
    return 23;
  } else if ((*sptr == 'Y') || (*sptr == 'y')) {
    return 24;
  } else {
    return -1;
  }
}

void cur_item(char* buf, char* sptr) {
  // no bounds-checking
  while ((*sptr != ' ') && (*sptr != '\t')) {
    *buf++ = *sptr++;
  }
  *buf = '\0';
}

void cur_item_n(char* buf, char* sptr) {
  while ((*sptr != ' ') && (*sptr != '\t') && (*sptr != '\n') && (*sptr != '\0')) {
    *buf++ = *sptr++;
  }
  *buf = '\0';
}

char* next_item(char* sptr) {
  if (!sptr) {
    return NULL;
  }
  while ((*sptr != ' ') && (*sptr != '\t')) {
    if (!(*sptr)) {
      return NULL;
    }
    sptr++;
  }
  while ((*sptr == ' ') || (*sptr == '\t')) {
    sptr++;
  }
  return sptr;
}

void fill_weights(double* weights, double* mafs, double exponent) {
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  double wtarr[MULTIPLEX_DIST_EXP / 2];
  double* wt;
#if __LP64__
  double twt[5];
  double twtf;
  __m128d* wpairs = (__m128d*)weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
#else
  int pp;
  int qq;
  double twt[7];
#endif
  for (ii = 0; ii < MULTIPLEX_DIST_EXP / 2; ii++) {
    wtarr[ii] = pow(2 * mafs[ii] * (1.0 - mafs[ii]), -exponent);
  }
  for (oo = 0; oo < 2; oo++) {
    wt = &(wtarr[7 * oo]);
#if __LP64__
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(wt[0] * 2, wt[0]);
#endif
    twt[0] = 0;
    for (ii = 0; ii < 4; ii += 1) {
      // bizarrely, turning the ii == 2 case into a memcpy doesn't actually
      // seem to speed this up
      if (ii & 1) {
	twt[0] += wt[6];
      }
      twt[1] = twt[0];
      for (jj = 0; jj < 4; jj += 1) {
	if (jj & 1) {
	  twt[1] += wt[5];
	}
	twt[2] = twt[1];
	for (kk = 0; kk < 4; kk += 1) {
	  if (kk & 1) {
	    twt[2] += wt[4];
	  }
	  twt[3] = twt[2];
	  for (mm = 0; mm < 4; mm++) {
	    if (mm & 1) {
	      twt[3] += wt[3];
	    }
	    twt[4] = twt[3];
	    for (nn = 0; nn < 4; nn++) {
	      if (nn & 1) {
		twt[4] += wt[2];
	      }
#if __LP64__
	      twtf = twt[4];
	      vpen = _mm_set1_pd(twtf);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
	      twtf += wt[1];
	      vpen = _mm_set1_pd(twtf);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
	      *wpairs = *(wpairs - 2);
	      wpairs++;
	      *wpairs = *(wpairs - 2);
	      wpairs++;
	      vpen = _mm_set1_pd(twtf + wt[1]);
	      *wpairs++ = _mm_add_pd(vpen, vfinal1);
	      *wpairs++ = _mm_add_pd(vpen, vfinal2);
#else
              twt[5] = twt[4];
              for (pp = 0; pp < 4; pp++) {
                if (pp & 1) {
                  twt[5] += wt[1];
                }
                twt[6] = twt[5];
                for (qq = 0; qq < 4; qq++) {
                  if (qq & 1) {
                    twt[6] += wt[0];
                  }
                  *weights++ = twt[6];
                }
              }
#endif
	    }
          }
	}
      }
    }
  }
#if __LP64__
  for (oo = 0; oo < 3; oo++) {
    wt = &(wtarr[14 + 6 * oo]);
    vfinal1 = _mm_set_pd(wt[0], 0.0);
    vfinal2 = _mm_set_pd(2 * wt[0], wt[0]);
    twt[0] = 0;
    for (ii = 0; ii < 4; ii += 1) {
      if (ii & 1) {
        twt[0] += wt[5];
      }
      twt[1] = twt[0];
      for (jj = 0; jj < 4; jj += 1) {
        if (jj & 1) {
          twt[1] += wt[4];
        }
        twt[2] = twt[1];
	for (kk = 0; kk < 4; kk += 1) {
          if (kk & 1) {
            twt[2] += wt[3];
          }
          twt[3] = twt[2];
          for (mm = 0; mm < 4; mm++) {
	    if (mm & 1) {
	      twt[3] += wt[2];
	    }
	    twtf = twt[3];
	    vpen = _mm_set1_pd(twtf);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
	    twtf += wt[1];
	    vpen = _mm_set1_pd(twtf);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
	    *wpairs = *(wpairs - 2);
	    wpairs++;
	    *wpairs = *(wpairs - 2);
	    wpairs++;
	    vpen = _mm_set1_pd(twtf + wt[1]);
	    *wpairs++ = _mm_add_pd(vpen, vfinal1);
	    *wpairs++ = _mm_add_pd(vpen, vfinal2);
          }
	}
      }
    }
  }
#endif
}

// Strangely, these functions optimize better than memset(arr, 0, x) under gcc.
inline void fill_long_zero(long* larr, size_t size) {
  long* lptr = &(larr[size]);
  while (larr < lptr) {
    *larr++ = 0;
  }
}

inline void fill_long_one(unsigned long* larr, size_t size) {
  unsigned long* lptr = &(larr[size]);
  while (larr < lptr) {
#if __LP64__
    *larr++ = ~(0LLU);
#else
    *larr++ = ~0;
#endif
  }
}

inline void fill_int_zero(int* iarr, size_t size) {
#if __LP64__
  fill_long_zero((long*)iarr, size >> 1);
  if (size & 1) {
    iarr[size - 1] = 0;
  }
#else
  fill_long_zero((long*)iarr, size);
#endif
}

inline void fill_int_one(unsigned int* iarr, size_t size) {
#if __LP64__
  fill_long_one((unsigned long*)iarr, size >> 1);
  if (size & 1) {
    iarr[size - 1] = -1;
  }
#else
  fill_long_one((unsigned long*)iarr, size);
#endif
}

inline void fill_double_zero(double* darr, size_t size) {
  double* dptr = &(darr[size]);
  while (darr < dptr) {
    *darr++ = 0.0;
  }
}

void fill_weights_r(double* weights, double* mafs, int var_std) {
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  // 20 markers to process in quintuplets, for 64-bit; 10, for 32-bit.
  // Each quintuplet of markers requires 40 wtarr entries, and induces
  // 2^15 writes to weights[].
  double wtarr_raw[BITCT2 * 5 + 1];
  double* wtarr = wtarr_raw;
  double twt;
  double twt2;
  double twt3;
  double twt4;
  double* wtptr;
  double mean;
  double mean_m1;
  double mean_m2;
  double mult = 1.0;
  double aux;
#if __LP64__
  __m128d* wpairs = (__m128d*)weights;
  __m128d vpen;
  __m128d vfinal1;
  __m128d vfinal2;
  __m128d vfinal3;
  __m128d vfinal4;
#endif
  if (((unsigned long)wtarr) & 15) {
    // force 16-byte alignment; can't do this at compile-time since stack
    // pointer has no 16-byte align guarantee
    wtarr++;
  }
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if (((mafs[ii] > 0.00000001) && (mafs[ii] < 0.99999999)) || (!var_std)) {
      if (mafs[ii] < 0.5) {
	mean = 2 * mafs[ii];
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * (1.0 - mafs[ii]));
        }
        aux = mean * mult;
	wtarr[ii * 8] = mean * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean_m2 * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m2 * mean_m1 * mult;
	wtarr[ii * 8 + 6] = mean_m2 * mean_m2 * mult;
      } else {
	mean = 2 * (1.0 - mafs[ii]);
	mean_m1 = mean - 1.0;
	mean_m2 = mean - 2.0;
        if (var_std) {
	  mult = 1 / (mean * mafs[ii]);
        }
        aux = mean_m2 * mult;
	wtarr[ii * 8] = mean_m2 * aux;
        wtarr[ii * 8 + 1] = 0;
	wtarr[ii * 8 + 2] = mean_m1 * aux;
	wtarr[ii * 8 + 3] = mean * aux;
	wtarr[ii * 8 + 4] = mean_m1 * mean_m1 * mult;
	wtarr[ii * 8 + 5] = mean_m1 * mean * mult;
	wtarr[ii * 8 + 6] = mean * mean * mult;
      }
    } else {
      if (mafs[ii] < 0.5) {
        wtarr[ii * 8] = 0;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = -1;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = INFINITY;
        wtarr[ii * 8 + 6] = INFINITY;
      } else {
        wtarr[ii * 8] = INFINITY;
        wtarr[ii * 8 + 1] = 0;
        wtarr[ii * 8 + 2] = INFINITY;
        wtarr[ii * 8 + 3] = -2;
        wtarr[ii * 8 + 4] = INFINITY;
        wtarr[ii * 8 + 5] = -1;
        wtarr[ii * 8 + 6] = 0;
      }
    }
    wtarr[ii * 8 + 7] = 0;
  }
  for (nn = 0; nn < BITCT / 16; nn++) {
    wtptr = &(wtarr[40 * nn]);
#if __LP64__
    vfinal1 = _mm_load_pd(wtptr);
    vfinal2 = _mm_load_pd(&(wtptr[2]));
    vfinal3 = _mm_load_pd(&(wtptr[4]));
    vfinal4 = _mm_load_pd(&(wtptr[6]));
#endif
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 32];
      for (jj = 0; jj < 8; jj++) {
        twt2 = twt + wtptr[jj + 24];
        for (kk = 0; kk < 8; kk++) {
          twt3 = twt2 + wtptr[kk + 16];
          for (mm = 0; mm < 8; mm++) {
            twt4 = twt3 + wtptr[mm + 8];
#if __LP64__
            vpen = _mm_set1_pd(twt4);
            *wpairs++ = _mm_add_pd(vpen, vfinal1);
            *wpairs++ = _mm_add_pd(vpen, vfinal2);
            *wpairs++ = _mm_add_pd(vpen, vfinal3);
            *wpairs++ = _mm_add_pd(vpen, vfinal4);
#else
            for (nn = 0; nn < 8; nn++) {
              *weights++ = twt4 + wtptr[nn];
            }
#endif
          }
        }
      }
    }
  }
}

double calc_wt_mean(double exponent, int lhi, int lli, int hhi) {
  double lcount = (double)lli + ((double)lhi * 0.5);
  long long tot = lhi + lli + hhi;
  double dtot = (double)tot;
  long long subcount = lli; // avoid 32-bit integer overflow
  double weight;
  if ((lcount == 0.0) || (dtot - lcount == 0.0)) {
    return 0.0;
  }
  weight = pow(2 * lcount * (dtot - lcount) / (dtot * dtot), -exponent);
  subcount = lhi * (subcount + hhi) + 2 * subcount * hhi;
  return (subcount * weight) / (double)(tot * (tot - 1) / 2);
}

double calc_wt_mean_maf(double exponent, double maf) {
  // assume Hardy-Weinberg equilibrium
  // homozygote frequencies: maf^2, (1-maf)^2
  // heterozygote frequency: 2maf(1-maf)
  double ll_freq = maf * maf;
  double lh_freq = 2 * maf * (1.0 - maf);
  double hh_freq = (1.0 - maf) * (1.0 - maf);
  double weight;
  if (lh_freq == 0.0) {
    return 0.0;
  }
  weight = pow(lh_freq, -exponent);
  return (lh_freq * (ll_freq + lh_freq) + 2 * ll_freq * hh_freq) * weight;
}

// ----- multithread globals -----
int indiv_ct;
int thread_ct;
double* rel_dists = NULL;
int* rel_missing = NULL;
int* idists;
double* dists = NULL;
char* pheno_c = NULL;
double* pheno_d = NULL;
unsigned char* ped_geno = NULL;
unsigned long* glptr;
double weights[2048 * BITCT];
unsigned int* weights_i = (unsigned int*)weights;
int thread_start[MAX_THREADS_P1];
int thread_start0[MAX_THREADS_P1];
int thread_start0_na[MAX_THREADS_P1];
double reg_tot_xy;
double reg_tot_x;
double reg_tot_y;
double reg_tot_xx;
int low_ct;
int high_ct;
int jackknife_iters;
int jackknife_d;
double calc_result[MAX_THREADS];
double calc_result2[MAX_THREADS];
double calc_result3[MAX_THREADS];
double calc_result4[MAX_THREADS];
double calc_result5[MAX_THREADS];
double calc_result6[MAX_THREADS];
double calc_result7[MAX_THREADS];
double calc_result8[MAX_THREADS];
double calc_result9[MAX_THREADS];
unsigned long* masks;
unsigned long* mmasks;
double* marker_weights;
unsigned int* marker_weights_i;
unsigned int* missing_tot_weights;
unsigned int* indiv_missing;

void update_rel_ibc(double* rel_ibc, unsigned long* geno, double* mafs, int ibc_type) {
  // first calculate weight array, then loop
  int ii;
  int jj;
  int kk;
  double twt;
  double* wtptr;
  double mult = 1.0;
  unsigned long ulii;
  double weights[BITCT * 10];
  double* weights1 = &(weights[64]);
  double* weights2 = &(weights[128]);
  double* weights3 = &(weights[192]);
  double* weights4 = &(weights[256]);
#if __LP64__
  double* weights5 = &(weights[320]);
  double* weights6 = &(weights[384]);
  double* weights7 = &(weights[448]);
  double* weights8 = &(weights[512]);
  double* weights9 = &(weights[576]);
#endif
  double wtarr[BITCT2 * 5];
  fill_double_zero(wtarr, BITCT2 * 5);
  double *wptr = weights;
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if ((mafs[ii] > 0.00000001) && (mafs[ii] < 0.99999999)) {
      if (ibc_type) {
        if (ibc_type == 2) {
          wtarr[ii * 8] = 2;
          wtarr[ii * 8 + 2] = 2.0 - 1.0 / (2 * mafs[ii] * (1.0 - mafs[ii]));
          wtarr[ii * 8 + 3] = 2;
        } else {
          twt = 2 * mafs[ii];
          if (ibc_type == 1) {
            mult = 1 / (twt * (1.0 - mafs[ii]));
          }
          wtarr[ii * 8] = twt * twt * mult;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt) * mult;
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt) * mult;
        }
      } else {
        twt = 1.0 - mafs[ii];
        mult = 1 / (mafs[ii] * twt);
        wtarr[ii * 8] = 1.0 + mafs[ii] * mafs[ii] * mult;
        wtarr[ii * 8 + 3] = 1.0 + twt * twt * mult;
      }
    } else {
      if (ibc_type) {
        if (ibc_type == -1) {
          twt = 2 * mafs[ii];
          wtarr[ii * 8] = twt * twt;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt);
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt);
        } else if (ibc_type == 1) {
	  wtarr[ii * 8 + 2] = INFINITY;
          if (mafs[ii] > 0.5) {
            wtarr[ii * 8] = INFINITY;
            wtarr[ii * 8 + 3] = 0;
          } else {
            wtarr[ii * 8] = 0;
            wtarr[ii * 8 + 3] = INFINITY;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 2] = -INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      } else {
        if (mafs[ii] > 0.5) {
          wtarr[ii * 8] = INFINITY;
          wtarr[ii * 8 + 3] = 1;
        } else {
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 3] = INFINITY;
        }
      }
    }
  }
  for (kk = 0; kk < (BITCT * 5) / 32; kk++) {
    wtptr = &(wtarr[16 * kk]);
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 8];
      for (jj = 0; jj < 8; jj++) {
        *wptr++ = twt + wtptr[jj];
      }
    }
  }
  for (ii = 0; ii < indiv_ct; ii++) {
    ulii = *geno++;
#if __LP64__
    *rel_ibc += weights9[ulii >> 54] + weights8[(ulii >> 48) & 63] + weights7[(ulii >> 42) & 63] + weights6[(ulii >> 36) & 63] + weights5[(ulii >> 30) & 63] + weights4[(ulii >> 24) & 63] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#else
    *rel_ibc += weights4[ulii >> 24] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#endif
    rel_ibc++;
  }
}

void exclude(unsigned char* exclude_arr, int loc, int* exclude_ct) {
  int maj = loc / 8;
  int min = 1 << (loc % 8);
  if (!(exclude_arr[maj] & min)) {
    exclude_arr[maj] |= min;
    *exclude_ct += 1;
  }
}

inline int excluded(unsigned char* exclude_arr, int loc) {
  return exclude_arr[loc / 8] & (1 << (loc % 8));
}

void exclude_multi(unsigned char* exclude_arr, int* new_excl, int indiv_ct, int* exclude_ct) {
  int ii;
  int true_loc = 0;
  for (ii = 0; ii < indiv_ct; ii++) {
    while (excluded(exclude_arr, true_loc)) {
      true_loc++;
    }
    if (new_excl[ii] == -1) {
      exclude(exclude_arr, true_loc, exclude_ct);
    }
    true_loc++;
  }
}

int marker_id_sort(char** sorted_ids_ptr, int** id_map_ptr, int unfiltered_marker_ct, unsigned char* marker_exclude, int marker_exclude_ct, char* marker_id, unsigned long max_marker_id_len) {
  // note that this allocates from wkspace
  int ii;
  int jj;
  *sorted_ids_ptr = (char*)wkspace_alloc((unfiltered_marker_ct - marker_exclude_ct) * max_marker_id_len);
  if (!(*sorted_ids_ptr)) {
    return RET_NOMEM;
  }
  *id_map_ptr = (int*)wkspace_alloc((unfiltered_marker_ct - marker_exclude_ct) * sizeof(int));
  if (!(*id_map_ptr)) {
    return RET_NOMEM;
  }
  jj = 0;
  for (ii = 0; ii < unfiltered_marker_ct; ii++) {
    if (excluded(marker_exclude, ii)) {
      continue;
    }
    memcpy(&((*sorted_ids_ptr)[jj * max_marker_id_len]), &(marker_id[ii * max_marker_id_len]), max_marker_id_len);
    (*id_map_ptr)[jj++] = ii;
  }
  qsort_ext(*sorted_ids_ptr, unfiltered_marker_ct - marker_exclude_ct, max_marker_id_len, strcmp_deref, (char*)(*id_map_ptr), sizeof(int));
  return 0;
}

int extract_exclude_markers(char* fname, char* sorted_ids, int sorted_ids_len, unsigned long max_marker_id_len, int* id_map, int unfiltered_marker_ct, unsigned char* marker_exclude, int* marker_exclude_ct_ptr, int do_exclude) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned char* marker_exclude_new = NULL;
  int marker_include_ct = 0;
  int ii;
  int jj;

  if (!do_exclude) {
    marker_exclude_new = wkspace_alloc((unfiltered_marker_ct + 7) / 8);
    if (!marker_exclude_new) {
      return RET_NOMEM;
    }
    memset(marker_exclude_new, 255, (unfiltered_marker_ct + 7) / 8);
  }
  infile = fopen(fname, "r");
  if (!infile) {
    printf("Error: Failed to open %s.\n", fname);
    wkspace_reset(wkspace_mark);
    return RET_OPENFAIL;
  }
  while (fgets(tbuf, MAXLINELEN, infile) != NULL) {
    cur_item_n(tbuf, tbuf);
    ii = bsearch_str(tbuf, sorted_ids, max_marker_id_len, 0, sorted_ids_len - 1);
    if (ii != -1) {
      jj = id_map[ii];
      if (do_exclude) {
        exclude(marker_exclude, jj, marker_exclude_ct_ptr);
      } else {
        if (!excluded(marker_exclude, jj)) {
          marker_exclude_new[jj / 8] &= ~(1 << (jj % 8));
          marker_include_ct++;
        }
      }
    }
  }
  if (!do_exclude) {
    memcpy(marker_exclude, marker_exclude_new, (unfiltered_marker_ct + 7) / 8);
    *marker_exclude_ct_ptr = unfiltered_marker_ct - marker_include_ct;
    wkspace_reset(wkspace_mark);
  }

  if (infile) {
    fclose(infile);
  }
  return 0;
}

int update_freq(char* freqname, FILE** freqfile_ptr, int unfiltered_marker_ct, unsigned char* marker_exclude, int marker_exclude_ct, char* marker_id, unsigned long max_marker_id_len, int* marker_chrom, char* marker_alleles, double* mafs) {
  unsigned char* wkspace_mark;
  char* sorted_ids;
  int* id_map;
  int ii;
  int jj;
  char cc;
  char* bufptr;
  char* bufptr2;
  double maf;
  *freqfile_ptr = fopen(freqname, "r");
  if (!(*freqfile_ptr)) {
    printf("Error: Failed to open %s.\n", freqname);
    return RET_OPENFAIL;
  }
  if (fgets(tbuf, MAXLINELEN, *freqfile_ptr) == NULL) {
    printf("Error: Empty --update-freq file.\n");
    return RET_INVALID_FORMAT;
  }
  wkspace_mark = wkspace_base;
  marker_id_sort(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_id, max_marker_id_len);
  if (!memcmp(tbuf, " CHR  ", 6)) {
    while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL) {
      bufptr = next_item(tbuf); // just initial spaces
      jj = marker_code(bufptr);
      bufptr = next_item(bufptr); // now at beginning of SNP name
      bufptr2 = next_item(bufptr);
      cur_item_n(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        if ((jj == marker_chrom[ii]) || (jj == 0) || (marker_chrom[ii] == 0)) {
          cc = *bufptr2;
          bufptr2 = next_item(bufptr2);
          bufptr = next_item(bufptr2);
          sscanf(bufptr, "%lg", &maf);
          if (marker_alleles[ii * 2] == cc) {
            if (marker_alleles[ii * 2 + 1] == *bufptr2) {
              mafs[ii] = 1.0 - maf;
            }
          } else if (marker_alleles[ii * 2 + 1] == cc) {
            if (marker_alleles[ii * 2] == *bufptr2) {
              mafs[ii] = maf;
            }
          }
        }
      }
    }
  } else {
    do {
      bufptr = next_item(tbuf);
      cur_item_n(tbuf, tbuf); // destructive read
      ii = bsearch_str(tbuf, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        cc = *bufptr;
        bufptr = next_item(bufptr);
        sscanf(bufptr, "%lg", &maf);
        if (marker_alleles[ii * 2] == cc) {
          mafs[ii] = 1.0 - maf;
        } else if (marker_alleles[ii * 2 + 1] == cc) {
          mafs[ii] = maf;
        }
      }
    } while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL);
  }
  fclose(*freqfile_ptr);
  *freqfile_ptr = NULL;
  wkspace_reset(wkspace_mark);
  return 0;
}

int is_missing(char* bufptr, int missing_pheno, int missing_pheno_len, int affection_01) {
  if ((atoi(bufptr) == missing_pheno) && is_space_or_eoln(bufptr[missing_pheno_len])) {
    return 1;
  } else if ((!affection_01) && (*bufptr == '0') && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

int eval_affection(char* bufptr, int missing_pheno, int missing_pheno_len, int affection_01) {
  if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
    return 1;
  } else if (((*bufptr == '0') || (*bufptr == '1') || ((*bufptr == '2') && (!affection_01))) && is_space_or_eoln(bufptr[1])) {
    return 1;
  }
  return 0;
}

void collapse_phenoc(char* pheno_c, unsigned char* indiv_exclude, int unfiltered_indiv_ct) {
  int ii = 0;
  int jj;
  while ((ii < unfiltered_indiv_ct) && (!excluded(indiv_exclude, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < unfiltered_indiv_ct) {
    if (!excluded(indiv_exclude, ii)) {
      pheno_c[jj++] = pheno_c[ii];
    }
  }
}

void collapse_phenod(double* pheno_d, unsigned char* indiv_exclude, int unfiltered_indiv_ct) {
  int ii = 0;
  int jj;
  while ((ii < unfiltered_indiv_ct) && (!excluded(indiv_exclude, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < unfiltered_indiv_ct) {
    if (!excluded(indiv_exclude, ii)) {
      pheno_d[jj++] = pheno_d[ii];
    }
  }
}

void collapse_copy_phenod(double *target, double* pheno_d, unsigned char* indiv_exclude, int unfiltered_indiv_ct) {
  int ii = 0;
  double* target_end = &(target[indiv_ct]);
  while (target < target_end) {
    while (excluded(indiv_exclude, ii)) {
      ii++;
    }
    *target++ = pheno_d[ii++];
  }
}

void pick_d(unsigned char* cbuf, unsigned int ct, unsigned int dd) {
  unsigned int ii;
  unsigned int jj;
  unsigned int kk;
  memset(cbuf, 0, ct);
  kk = 1073741824 % ct;
  kk = (kk * 4) % ct;
  for (ii = 0; ii < dd; ii++) {
    do {
      do {
        jj = genrand_int32();
      } while (jj < kk);
      jj %= ct;
    } while (cbuf[jj]);
    cbuf[jj] = 1;
  }
}

void pick_d_small(unsigned char* tmp_cbuf, int* ibuf, unsigned int ct, unsigned int dd) {
  int ii;
  int jj = 0;
  pick_d(tmp_cbuf, ct, dd);
  for (ii = 0; ii < ct; ii++) {
    if (tmp_cbuf[ii]) {
      ibuf[jj++] = ii;
    }
  }
  ibuf[jj] = ct;
}

#if __LP64__
typedef union {
  __m128i vi;
  unsigned long u8[2];
} __uni16;

// SSE2 implementations of Lauradoux/Walisch popcount, combined with xor to
// handle Hamming distance, and masking to handle missingness.
// Note that the size of the popcounted buffer is a hardcoded constant
// (specifically, (MULTIPLEX_DIST / BITCT) * 16 bytes).  The current code
// assumes (MULTIPLEX / BITCT) is a multiple of 3, and no greater than 30.
static inline unsigned int popcount_xor_1mask_multibyte(__m128i** xor1p, __m128i* xor2, __m128i** maskp) {
  const __m128i m1 = {0x5555555555555555LLU, 0x5555555555555555LLU};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  const __m128i m16 = {0x0000ffff0000ffffLLU, 0x0000ffff0000ffffLLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  // while loop seems to optimize better than for here
  while (xor2 < xor2_end) {
    count1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    count2 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    half1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), *((*maskp)++));
    half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
    half1 = _mm_and_si128(half1, m1);
    // Two bits can represent values from 0-3, so make each pair in count1 and
    // count2 store a partial bitcount covering themselves AND another bit from
    // elsewhere.
    count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
    count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
    count1 = _mm_add_epi64(count1, half1);
    count2 = _mm_add_epi64(count2, half2);
    // Four bits represent 0-15, so we can add four 0-3 partial bitcounts
    // together.
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    // Accumulator stores sixteen 0-255 counts in parallel.
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  }
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  // can get away with this with MULTIPLEX_DIST <= 960, since the 8-bit counts
  // are guaranteed to be <= 120, thus adding two together does not overflow
  // 255.
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return (unsigned int)(acc.u8[0] + acc.u8[1]);
}

static inline unsigned int popcount_xor_2mask_multibyte(__m128i** xor1p, __m128i* xor2, __m128i** mask1p, __m128i* mask2) {
  const __m128i m1 = {0x5555555555555555LLU, 0x5555555555555555LLU};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  const __m128i m16 = {0x0000ffff0000ffffLLU, 0x0000ffff0000ffffLLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  while (xor2 < xor2_end) {
    count1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    count2 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    half1 = _mm_and_si128(_mm_xor_si128(*((*xor1p)++), *xor2++), _mm_and_si128(*((*mask1p)++), *mask2++));
    half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
    half1 = _mm_and_si128(half1, m1);
    count1 = _mm_sub_epi64(count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
    count2 = _mm_sub_epi64(count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
    count1 = _mm_add_epi64(count1, half1);
    count2 = _mm_add_epi64(count2, half2);
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  }
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return (unsigned int)(acc.u8[0] + acc.u8[1]);
}
#else
// Simple 32-bit popcount implementation (using a 16-bit lookup table) for
// testing purposes.  If there is a serious need for a fast 32-bit version, one
// can make minor modifications to the code above, or revert to standard
// Walisch/Lauradoux if SSE2 is unavailable.
static inline unsigned int popcount_xor_1mask_multibyte(unsigned long** xor1p, unsigned long* xor2, unsigned long** maskp) {
  unsigned long* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  unsigned int bit_count = 0;
  unsigned long ulii;
  while (xor2 < xor2_end) {
    ulii = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    bit_count += popcount[ulii >> 16] + popcount[ulii & 65535];
  }
  return bit_count;
}

static inline unsigned int popcount_xor_2mask_multibyte(unsigned long** xor1p, unsigned long* xor2, unsigned long** mask1p, unsigned long* mask2) {
  unsigned long* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  unsigned int bit_count = 0;
  unsigned long ulii;
  while (xor2 < xor2_end) {
    ulii = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    bit_count += popcount[ulii >> 16] + popcount[ulii & 65535];
  }
  return bit_count;
}
#endif

void incr_dists_i(int* idists, unsigned long* geno, int tidx) {
#if __LP64__
  __m128i* glptr;
  __m128i* glptr2;
  __m128i* mptr;
  __m128i* mcptr_start;
  unsigned long* lptr;
#else
  unsigned long* glptr;
  unsigned long* glptr2;
  unsigned long* mptr;
  unsigned long* mcptr_start;
#endif
  int ii;
  int jj;
  unsigned long mask_fixed;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    jj = ii * (MULTIPLEX_2DIST / BITCT);
#if __LP64__
    glptr = (__m128i*)geno;
    glptr2 = (__m128i*)(&(geno[jj]));
    lptr = &(masks[jj]);
    mcptr_start = (__m128i*)lptr;
    mask_fixed = *lptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *lptr++;
    }
    mptr = (__m128i*)masks;
#else
    glptr = geno;
    glptr2 = &(geno[jj]);
    mcptr_start = &(masks[jj]);
    mptr = mcptr_start;
    mask_fixed = *mptr++;
    for (jj = 0; jj < MULTIPLEX_2DIST / BITCT - 1; jj++) {
      mask_fixed &= *mptr++;
    }
    mptr = masks;
#endif
    if (~mask_fixed) {
      while (glptr < glptr2) {
	*idists += popcount_xor_2mask_multibyte(&glptr, glptr2, &mptr, mcptr_start);
	idists++;
      }
    } else {
      while (glptr < glptr2) {
	*idists += popcount_xor_1mask_multibyte(&glptr, glptr2, &mptr);
	idists++;
      }
    }
  }
}

void* calc_idist_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  incr_dists_i(&(idists[(ii * (ii - 1)) / 2]), (unsigned long*)ped_geno, (int)tidx);
  return NULL;
}

void incr_dists(double* dists, unsigned long* geno, int tidx) {
  unsigned long* glptr;
  unsigned long ulii;
  unsigned long mask_fixed;
  unsigned long uljj;
  unsigned long* mptr;
  double* weights1 = &(weights[16384]);
#if __LP64__
  double* weights2 = &(weights[32768]);
  double* weights3 = &(weights[36864]);
  double* weights4 = &(weights[40960]);
#endif
  int ii;
  int jj;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    glptr = geno;
    ulii = geno[ii];
    mptr = masks;
    mask_fixed = masks[ii];
#if __LP64__
    if (mask_fixed == 0xffffffffffffffffLLU) {
      for (jj = 0; jj < ii; jj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + weights[uljj & 16383];
	dists++;
      }
    } else {
      for (jj = 0; jj < ii; jj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
        *dists += weights4[uljj >> 52] + weights3[(uljj >> 40) & 4095] + weights2[(uljj >> 28) & 4095] + weights1[(uljj >> 14) & 16383] + weights[uljj & 16383];
	dists++;
      }
    }
#else
    if (mask_fixed == 0x0fffffff) {
      for (jj = 0; jj < ii; jj++) {
	uljj = (*glptr++ ^ ulii) & (*mptr++);
	*dists += weights1[uljj >> 14] + weights[uljj & 16383];
	dists++;
      }
    } else {
      for (jj = 0; jj < ii; jj++) {
	uljj = (*glptr++ ^ ulii) & (mask_fixed & (*mptr++));
	*dists += weights1[uljj >> 14] + weights[uljj & 16383];
	dists++;
      }
    }
#endif
  }
}

void* calc_dist_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  incr_dists(&(dists[(ii * (ii - 1)) / 2]), (unsigned long*)ped_geno, (int)tidx);
  return NULL;
}

void decr_dist_missing(unsigned int* mtw, int tidx) {
  unsigned long* glptr;
  unsigned long ulii;
  unsigned long uljj;
  int ii;
  int jj;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    glptr = mmasks;
    ulii = mmasks[ii];
    if (ulii) {
      for (jj = 0; jj < ii; jj++) {
	uljj = (*glptr++) & ulii;
        while (uljj) {
          mtw[jj] += weights_i[__builtin_ctzl(uljj)];
          uljj &= uljj - 1;
        }
      }
    }
    mtw = &(mtw[ii]);
  }
}

void* calc_distm_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  decr_dist_missing(&(missing_tot_weights[(ii * (ii - 1)) / 2]), (int)tidx);
  return NULL;
}

void incr_dists_r(double* dists, unsigned long* geno, unsigned long* masks, int tidx, double* weights) {
  unsigned long* glptr;
  unsigned long* maskptr;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long basemask;
  double* weights1 = &(weights[32768]);
#if __LP64__
  double* weights2 = &(weights[65536]);
  double* weights3 = &(weights[98304]);
#endif
  int ii;
  int jj;
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    glptr = geno;
    ulii = geno[ii];
    maskptr = masks;
    basemask = masks[ii];
    if (!basemask) {
      for (jj = 0; jj < ii; jj++) {
	uljj = ((*glptr++) + ulii) | (*maskptr++);
#if __LP64__
	*dists += weights3[uljj >> 45] + weights2[(uljj >> 30) & 32767] + weights1[(uljj >> 15) & 32767] + weights[uljj & 32767];
#else
	*dists += weights1[uljj >> 15] + weights[uljj & 32767];
#endif
	dists++;
      }
    } else {
      for (jj = 0; jj < ii; jj++) {
        uljj = ((*glptr++) + ulii) | ((*maskptr++) | basemask);
#if __LP64__
	*dists += weights3[uljj >> 45] + weights2[(uljj >> 30) & 32767] + weights1[(uljj >> 15) & 32767] + weights[uljj & 32767];
#else
	*dists += weights1[uljj >> 15] + weights[uljj & 32767];
#endif
	dists++;
      }
    }
  }
}

void* calc_rel_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  incr_dists_r(&(rel_dists[(ii * (ii - 1)) / 2]), (unsigned long*)ped_geno, masks, (int)tidx, weights);
  return NULL;
}

void incr_dists_rm(int* idists, int tidx) {
  // count missing intersection, optimized for sparsity
  unsigned long* glptr;
  unsigned long ulii;
  unsigned long uljj;
  int ii;
  int jj;
  for (ii = thread_start0[tidx]; ii < thread_start0[tidx + 1]; ii++) {
    glptr = mmasks;
    ulii = mmasks[ii];
    if (ulii) {
      for (jj = 0; jj < ii; jj++) {
        uljj = (*glptr++) & ulii;
        while (uljj) {
          idists[jj] += 1;
          uljj &= uljj - 1;
        }
      }
    }
    idists = &(idists[ii]);
  }
}

void* calc_relm_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start0[tidx];
  incr_dists_rm(&(rel_missing[(ii * (ii - 1)) / 2]), (int)tidx);
  return NULL;
}

int triangle_bsearch(long long cur_prod, int modif, int lbound, int ubound) {
  int center;
  long long vv;
  if (lbound == ubound) {
    return lbound;
  }
  center = (lbound + ubound) / 2;
  vv = center;
  vv = vv * (vv + modif);
  if (vv < cur_prod) {
    return triangle_bsearch(cur_prod, modif, center + 1, ubound);
  } else {
    return triangle_bsearch(cur_prod, modif, lbound, center);
  }
}

// set align to 2 for no alignment
void triangle_fill(int* target_arr, int ct, int pieces, int start, int align) {
  long long ct_tr = (long long)ct;
  long long cur_prod = 0;
  int modif = 1 - start * 2;
  int cur_piece = 1;
  int lbound = start;
  int ii;
  int align_m1;
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = start;
  target_arr[pieces] = ct;
  ct_tr = (ct_tr * (ct_tr + modif)) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_bsearch(cur_prod, modif, start, ct);
    ii = (lbound - start) & align_m1;
    if ((ii) && (ii != align_m1)) {
      lbound = start + ((lbound - start) | align_m1);
    }
    if (lbound > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

void groupdist_jack(int* ibuf, double* returns) {
  int* iptr = ibuf;
  int* jptr;
  int ii;
  double* dptr = dists;
  char* cptr;
  char* cptr2;
  char cc;
  char cc2;
  double neg_tot_aa = 0.0;
  double neg_tot_au = 0.0;
  double neg_tot_uu = 0.0;
  int neg_a = 0;
  int neg_u = 0;
  double dxx;
  for (ii = 0; ii < indiv_ct; ii++) {
    cc = pheno_c[ii];
    if (cc == -1) {
      continue;
    }
    dptr = &(dists[(ii * (ii - 1)) / 2]);
    if (ii == *iptr) {
      cptr = pheno_c;
      cptr2 = &(pheno_c[ii]);
      if (cc) {
        neg_a++;
        while (cptr < cptr2) {
          cc2 = *cptr++;
          dxx = *dptr++;
          if (cc2 == 1) {
            neg_tot_aa += dxx;
          } else if (!cc2) {
            neg_tot_au += dxx;
          }
        }
      } else {
        neg_u++;
        while (cptr < cptr2) {
          cc2 = *cptr++;
          dxx = *dptr++;
          if (cc2 == 1) {
            neg_tot_au += dxx;
          } else if (!cc2) {
            neg_tot_uu += dxx;
          }
        }
      }
      iptr++;
    } else {
      jptr = ibuf;
      while (jptr < iptr) {
        cc2 = pheno_c[*jptr];
        dxx = dptr[*jptr];
	if ((cc == 1) && (cc2 == 1)) {
	  neg_tot_aa += dxx;
	} else if ((cc == 0) && (cc2 == 0)) {
	  neg_tot_uu += dxx;
	} else {
	  neg_tot_au += dxx;
	}
        jptr++;
      }
    }
  }
  returns[0] = (reg_tot_x - neg_tot_aa) / (double)(((high_ct - neg_a) * (high_ct - neg_a - 1)) / 2);
  returns[1] = (reg_tot_xy - neg_tot_au) / (double)((high_ct - neg_a) * (low_ct - neg_u));
  returns[2] = (reg_tot_y - neg_tot_uu) / (double)(((low_ct - neg_u) * (low_ct - neg_u - 1)) / 2);
}

void small_remap(int* ibuf, unsigned int ct, unsigned int dd) {
  int* ibuf_end = &(ibuf[dd]);
  int missings = 0;
  int curpos = 0;
  do {
    if (pheno_c[curpos] == -1) {
      missings++;
    } else if (*ibuf == curpos - missings) {
      *ibuf = curpos;
      ibuf++;
    }
    curpos++;
  } while (ibuf < ibuf_end);
}

void* groupdist_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(ped_geno[tidx * CACHEALIGN(high_ct + low_ct + (jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(ped_geno[tidx * CACHEALIGN(high_ct + low_ct + (jackknife_d + 1) * sizeof(int)) + jackknife_d * sizeof(int)]);
  unsigned long long ulii;
  unsigned long long uljj = jackknife_iters / 100;
  double returns[3];
  double results[9];
  double new_old_diff[3];
  fill_double_zero(results, 9);
  for (ulii = 0; ulii < jackknife_iters; ulii++) {
    pick_d_small(cbuf, ibuf, high_ct + low_ct, jackknife_d);
    if (high_ct + low_ct < indiv_ct) {
      small_remap(ibuf, high_ct + low_ct, jackknife_d);
    }
    groupdist_jack(ibuf, returns);
    if (ulii > 0) {
      new_old_diff[0] = returns[0] - results[0];
      new_old_diff[1] = returns[1] - results[1];
      new_old_diff[2] = returns[2] - results[2];
      results[0] += new_old_diff[0] / (ulii + 1); // AA mean
      results[1] += new_old_diff[1] / (ulii + 1); // AU mean
      results[2] += new_old_diff[2] / (ulii + 1); // UU mean
      results[3] += (returns[0] - results[0]) * new_old_diff[0]; // AA var
      results[4] += (returns[1] - results[1]) * new_old_diff[1]; // AU var
      results[5] += (returns[2] - results[2]) * new_old_diff[2]; // UU var
      results[6] += (returns[0] - results[0]) * new_old_diff[1]; // AA-AU cov
      results[7] += (returns[0] - results[0]) * new_old_diff[2]; // AA-UU cov
      results[8] += (returns[1] - results[1]) * new_old_diff[2]; // AU-UU cov
    } else {
      results[0] += returns[0];
      results[1] += returns[1];
      results[2] += returns[2];
    }
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100) / jackknife_iters;
      printf("\r%lld%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1) * jackknife_iters) / 100;
    }
  }
  // don't write until end, to avoid false sharing
  calc_result[tidx] = results[0];
  calc_result2[tidx] = results[1];
  calc_result3[tidx] = results[2];
  calc_result4[tidx] = results[3];
  calc_result5[tidx] = results[4];
  calc_result6[tidx] = results[5];
  calc_result7[tidx] = results[6];
  calc_result8[tidx] = results[7];
  calc_result9[tidx] = results[8];
  return NULL;
}

double regress_jack(int* ibuf) {
  int* iptr = ibuf;
  int* jptr;
  int ii;
  int jj;
  double* dptr = dists;
  double* dptr2;
  double* dptr3;
  double neg_tot_xy = 0.0;
  double neg_tot_x = 0.0;
  double neg_tot_y = 0.0;
  double neg_tot_xx = 0.0;
  double dxx;
  double dxx1;
  double dyy;
  if (*iptr == 0) {
    iptr++;
  }
  for (ii = 1; ii < indiv_ct; ii++) {
    dxx1 = pheno_d[ii];
    if (ii == *iptr) {
      dptr2 = pheno_d;
      dptr3 = &(pheno_d[ii]);
      do {
	dxx = (dxx1 + (*dptr2++)) * 0.5;
	dyy = *dptr++;
	neg_tot_xy += dxx * dyy;
	neg_tot_x += dxx;
	neg_tot_y += dyy;
	neg_tot_xx += dxx * dxx;
      } while (dptr2 < dptr3);
      iptr++;
    } else {
      jptr = ibuf;
      do {
        jj = *jptr++;
        dxx = (dxx1 + pheno_d[jj]) * 0.5;
        dyy = dptr[jj];
	neg_tot_xy += dxx * dyy;
	neg_tot_x += dxx;
	neg_tot_y += dyy;
	neg_tot_xx += dxx * dxx;
      } while (jptr < iptr);
      dptr = &(dists[(ii * (ii + 1)) / 2]);
    }
  }
  dxx = reg_tot_x - neg_tot_x;
  dyy = indiv_ct - jackknife_d;
  dyy = dyy * (dyy - 1.0) * 0.5;
  return ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_y - neg_tot_y) / dyy) / ((reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

void* regress_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(ped_geno[tidx * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(ped_geno[tidx * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int)) + jackknife_d * sizeof(int)]);
  unsigned long long ulii;
  unsigned long long uljj = jackknife_iters / 100;
  double sum = 0.0;
  double sum_sq = 0.0;
  double dxx;
  for (ulii = 0; ulii < jackknife_iters; ulii++) {
    pick_d_small(cbuf, ibuf, indiv_ct, jackknife_d);
    dxx = regress_jack(ibuf);
    sum += dxx;
    sum_sq += dxx * dxx;
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100) / jackknife_iters;
      printf("\r%lld%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1) * jackknife_iters) / 100;
    }
  }
  calc_result[tidx] = sum;
  calc_result2[tidx] = sum_sq;
  return NULL;
}

#ifndef NOLAPACK
// Replaces matrix[][] with mult_val * matrix[][] + add_val * I.
// Multithreading doesn't help here.
void matrix_const_mult_add(double* matrix, double mult_val, double add_val) {
  int ii;
  int loop_end = indiv_ct - 1;
  int jj;
  double* dptr = matrix;
#if __LP64__
  __m128d* vptr;
  __m128d v_mult_val = _mm_set1_pd(mult_val);
#endif
  for (ii = 0; ii < loop_end; ii++) {
    *dptr = (*dptr) * mult_val + add_val;
    dptr++;
#if __LP64__
    if ((unsigned long)dptr & 8) {
      *dptr *= mult_val;
      dptr++;
      jj = 1;
    } else {
      jj = 0;
    }
    vptr = (__m128d*)dptr;
    while (jj < loop_end) {
      *vptr = _mm_mul_pd(*vptr, v_mult_val);
      vptr++;
      jj += 2;
    }
    dptr = (double*)vptr;
    if (jj < indiv_ct) {
      *dptr *= mult_val;
      dptr++;
    }
#else
    for (jj = 0; jj < indiv_ct; jj++) {
      *dptr *= mult_val;
      dptr++;
    }
#endif
  }
  *dptr = (*dptr) * mult_val + add_val;
}

// sums[idx] = matrix[idx][1] + matrix[idx][2] + ....  matrix assumed to be
// symmetric, and only FORTRAN upper right is referenced.
void matrix_row_sum_ur(double* sums, double* matrix) {
  int ii;
  double* dptr;
  double acc;
  double* sptr_end;
  double* sptr;
  fill_double_zero(sums, indiv_ct);
  for (ii = 0; ii < indiv_ct; ii++) {
    dptr = &(matrix[ii * indiv_ct]);
    acc = 0.0;
    sptr_end = &(sums[ii]);
    sptr = sums;
    while (sptr < sptr_end) {
      acc += *dptr;
      *sptr += *dptr++;
      sptr++;
    }
    *sptr += acc + *dptr;
  }
}

// one-trait REML via EM.
//
// wkbase is assumed to have space for three cache-aligned
// indiv_ct * indiv_ct double matrices plus three more rows.  The unpacked
// relationship matrix is stored in the SECOND slot.
void reml_em_one_trait(double* wkbase, double* pheno, double* covg_ref, double* cove_ref, double tol, int strict) {
  double ll_change;
  long long mat_offset = indiv_ct;
  double* rel_dists;
#ifdef __LP64__
  int* irow;
  int lwork;
#else
  long int* irow;
  long int lwork;
#endif
  double* row;
  double* row2;
  double* work;
  double* dptr;
  double* dptr2;
  double* matrix_pvg;
#ifdef __APPLE__
#ifdef __LP64__
  int info;
#else
  long int indiv_ct_li = indiv_ct;
  long int info;
#endif
#endif
  double dxx;
  double dyy;
  double max_jump;
  double dzz;
  double dlg;
  double dle;
  double covg_cur_change = 1.0;
  double cove_cur_change = 1.0;
  double covg_last_change;
  double cove_last_change;
  double indiv_ct_d = 1 / (double)indiv_ct;
  int ii;
  int jj;
  mat_offset = CACHEALIGN_DBL(mat_offset * mat_offset);
  rel_dists = &(wkbase[mat_offset]);
  row = &(wkbase[mat_offset * 3]);
#ifdef __LP64__
  irow = (int*)row;
#else
  irow = (long int*)row;
#endif
  row2 = &(row[indiv_ct]);
  work = &(wkbase[mat_offset * 2]);
  lwork = mat_offset;
  matrix_pvg = work;
  if (!lwork) {
    lwork = CACHELINE_DBL;
  }
  fill_double_zero(matrix_pvg, mat_offset);
  fill_double_zero(row, indiv_ct);
  fill_double_zero(row2, indiv_ct);
  do {
    memcpy(wkbase, rel_dists, mat_offset * sizeof(double));
    matrix_const_mult_add(wkbase, *covg_ref, *cove_ref);
#ifdef __APPLE__
#if __LP64__
    dgetrf_(&indiv_ct, &indiv_ct, wkbase, &indiv_ct, irow, &info);
    dgetri_(&indiv_ct, wkbase, &indiv_ct, irow, work, &lwork, &info);
#else
    dgetrf_(&indiv_ct_li, &indiv_ct_li, wkbase, &indiv_ct_li, irow, &info);
    dgetri_(&indiv_ct_li, wkbase, &indiv_ct_li, irow, work, &lwork, &info);
#endif
#else
    clapack_dgetrf(CblasColMajor, indiv_ct, indiv_ct, wkbase, indiv_ct, irow);
    clapack_dgetri(CblasColMajor, indiv_ct, wkbase, indiv_ct, irow);
#endif
    matrix_row_sum_ur(row, wkbase);
    dxx = 0.0;
    dptr = row;
    dptr2 = &(row[indiv_ct]);
    while (dptr < dptr2) {
      dxx += *dptr++;
    }
    dxx = -1 / dxx;
    cblas_dger(CblasColMajor, indiv_ct, indiv_ct, dxx, row, 1, row, 1, wkbase, indiv_ct);
    // unfortunately, cblas_dsymm is much worse than cblas_dgemm on OS X
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, indiv_ct, indiv_ct, indiv_ct, 1.0, wkbase, indiv_ct, rel_dists, indiv_ct, 0.0, matrix_pvg, indiv_ct);
    dlg = 0.0;
    dle = 0.0;
    jj = indiv_ct + 1;
    for (ii = 0; ii < indiv_ct; ii++) {
      dlg -= matrix_pvg[ii * jj];
      dle -= wkbase[ii * jj];
    }
    cblas_dsymv(CblasColMajor, CblasUpper, indiv_ct, 1.0, wkbase, indiv_ct, pheno, 1, 0.0, row2, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, indiv_ct, 1.0, matrix_pvg, indiv_ct, row2, 1, 0.0, row, 1);
    dlg += cblas_ddot(indiv_ct, pheno, 1, row, 1);
    cblas_dsymv(CblasColMajor, CblasUpper, indiv_ct, 1.0, wkbase, indiv_ct, row2, 1, 0.0, row, 1);
    dle += cblas_ddot(indiv_ct, pheno, 1, row, 1);
    covg_last_change = covg_cur_change;
    cove_last_change = cove_cur_change;
    covg_cur_change = (*covg_ref) * (*covg_ref) * dlg * indiv_ct_d;
    cove_cur_change = (*cove_ref) * (*cove_ref) * dle * indiv_ct_d;
    if (strict) {
      max_jump = 1.0;
    } else {
      // acceleration factor:
      // min(half covg distance to 0 or 1, cove distance to 0 or 1, pi/4 divided
      // by last angular change, 1.0 / (1 - ratio of last two step lengths),
      // MAX_EM_ACCEL)
      dxx = atan2(covg_last_change, cove_last_change) - atan2(covg_cur_change, cove_cur_change);
      if (dxx < 0.0) {
	dxx = -dxx;
      }
      if (dxx > PI) {
	dxx = 2 * PI - dxx;
      }
      dyy = sqrt((covg_cur_change * covg_cur_change + cove_cur_change * cove_cur_change) / (covg_last_change * covg_last_change + cove_last_change * cove_last_change));
      if (covg_cur_change < 0.0) {
	max_jump = *covg_ref * (-0.5) / covg_cur_change;
      } else {
	max_jump = (1.0 - *covg_ref) * 0.5 / covg_cur_change;
      }
      if (cove_cur_change < 0.0) {
	dzz = *cove_ref * (-0.5) / cove_cur_change;
      } else {
	dzz = (1.0 - *cove_ref) * 0.5 / cove_cur_change;
      }
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      dzz = (PI / 4) / dxx;
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      if (dyy < 1.0) {
	dzz = 1 / (1.0 - dyy);
      }
      if (dzz < max_jump) {
	max_jump = dzz;
      }
      if (max_jump < 1.0) {
	max_jump = 1.0;
      } else if (max_jump > MAX_EM_ACCEL) {
	max_jump = MAX_EM_ACCEL;
      }
    }
    *covg_ref += covg_cur_change * max_jump;
    *cove_ref += cove_cur_change * max_jump;
    ll_change = (covg_cur_change * dlg) + (cove_cur_change * dle);
    printf("\b\b\b\b\b\b      \rcovg: %g  cove: %g  EM step log likelihood change: %g", *covg_ref, *cove_ref, ll_change);
    fflush(stdout);
  } while (ll_change > tol);
  printf("\n");
}
#endif // NOLAPACK

inline int distance_req(int calculation_type) {
  return ((calculation_type & CALC_DISTANCE_MASK) || (calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE));
}

inline int relationship_req(int calculation_type) {
  return ((calculation_type & CALC_RELATIONSHIP_MASK) || (calculation_type & CALC_UNRELATED_HERITABILITY) || (calculation_type & CALC_REL_CUTOFF));
}

inline int relationship_or_ibc_req(int calculation_type) {
  return (relationship_req(calculation_type) || (calculation_type & CALC_IBC));
}

int wdist(char* pedname, char* mapname, char* famname, char* phenoname, char* extractname, char* excludename, char* filtername, char* freqname, char* makepheno_str, int filter_type, char* filterval, int mfilter_col, int make_bed, int ped_col_1, int ped_col_34, int ped_col_5, int ped_col_6, int ped_col_7, char missing_geno, int missing_pheno, int mpheno_col, char* phenoname_str, int prune, int affection_01, unsigned int chrom_mask, double exponent, double min_maf, double geno_thresh, double mind_thresh, double hwe_thresh, double grm_cutoff, int tail_pheno, double tail_bottom, double tail_top, char* outname, int calculation_type, int groupdist_iters, int groupdist_d, int regress_iters, int regress_d, double unrelated_herit_tol, double unrelated_herit_covg, double unrelated_herit_cove, int ibc_type) {
  FILE* pedfile = NULL;
  FILE* mapfile = NULL;
  FILE* famfile = NULL;
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  FILE* phenofile = NULL;
  FILE* filterfile = NULL;
  FILE* bedtmpfile = NULL;
  FILE* bimtmpfile = NULL;
  FILE* famtmpfile = NULL;
  FILE* freqfile = NULL;
  int unfiltered_marker_ct = 0;
  int unfiltered_marker_ct4 = 0;
  int marker_ct;
  char* outname_end;
  unsigned char* pedbuf = NULL;
  unsigned char* marker_exclude = NULL;
  long long* line_locs = NULL;
  int max_people;
  int pedbuflen;
  unsigned long max_marker_id_len = 0;
  double* mafs = NULL;
  int ped_recalc_len = 0;
  char* fgets_return;
  int unfiltered_indiv_ct = 0;
  int unfiltered_indiv_ct4 = 0;
  unsigned char* indiv_exclude = NULL;
  int indiv_exclude_ct = 0;
  int ii;
  int jj = 0;
  int kk = 0;
  int mm;
  int nn = 0;
  int oo = 0;
  int pp;
  int qq;
  unsigned int uii = 0;
  unsigned int ujj;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long ulkk;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  long long llxx;
  int* marker_chrom = NULL;
  char* marker_id = NULL;
  // int* marker_pos = NULL;
  char* marker_alleles = NULL;
  char* marker_alleles_tmp = NULL;
  int* marker_allele_cts = NULL;
  int* person_missing_cts = NULL;
  char* bufptr;
  int retval;
  int map_cols = 3;
  int affection = 0;
  double* phenor_d = NULL;
  char* phenor_c = NULL;
  char* person_id = NULL;
  int max_person_id_len = 4;
  unsigned char* wkspace_mark = NULL;
  unsigned int* giptr;
  unsigned int* giptr2;
  unsigned int* giptr3;
  char* cptr = NULL;
  char* cptr2;
  // char* cptr3;
  unsigned char* gptr;
  unsigned char* gptr2;
  unsigned long* glptr2;
  unsigned long* glptr3;
  int* iwptr;
  int* iptr;
  long long dists_alloc = 0;
  int geno_window_size;
  char cc;
  double* dist_ptr = NULL;
  double* dptr2;
  double* dptr3 = NULL;
  double* dptr4 = NULL;
  double* rel_ibc;
  long long last_tell;
  int binary_files = 0;
  int maf_int_thresh;
  int geno_int_thresh;
  int mind_int_thresh;
  int marker_exclude_ct = 0;
  int pheno_lines = 0;
  int makepheno_all = 0;
  int filter_lines = 0;
  int snp_major = 0;
  int max_pid_len = 0;
  int max_id_len = 4;
  char* pid_list = NULL;
  char* id_list = NULL;
  double maf_buf[MULTIPLEX_DIST];
  unsigned int wtbuf[MULTIPLEX_DIST];
  double missing_phenod = (double)missing_pheno;
  int missing_pheno_len = 1;
  int var_std = 1;
  int* hwe_ll;
  int* hwe_lh;
  int* hwe_hh;
  int* hwe_u_ll;
  int* hwe_u_lh;
  int* hwe_u_hh;
  int* hwe_a_ll;
  int* hwe_a_lh;
  int* hwe_a_hh;
  int hwe_lli;
  int hwe_lhi;
  int hwe_hhi;
  int hwe_u_lli;
  int hwe_u_lhi;
  int hwe_u_hhi;
  int hwe_a_lli;
  int hwe_a_lhi;
  int hwe_a_hhi;
  int multiplex = 0;
  int bin_pheno;
  double* ll_pool;
  double* lh_pool;
  double* hh_pool;
  double* ll_poolp;
  double* lh_poolp;
  double* hh_poolp;
  double dhh_ssq;
  double dhl_ssq;
  double dll_ssq;
  double dhh_sd;
  double dhl_sd;
  double dll_sd;
  int ll_size;
  int lh_size;
  int hh_size;
  double ll_med;
  double lh_med;
  double hh_med;
  pthread_t threads[MAX_THREADS];
  int exp0 = (exponent == 0.0);

  ii = missing_pheno;
  if (ii < 0) {
    ii = -ii;
    missing_pheno_len++;
  }
  while (ii > 9) {
    ii /= 10;
    missing_pheno_len++;
  }
  if (calculation_type & CALC_RELATIONSHIP_NO_VAR_STD) {
    var_std = 0;
    ibc_type = -1;
  }

  retval = RET_OPENFAIL;
  pedfile = fopen(pedname, "r");
  if (!pedfile) {
    printf("Error: Failed to open %s.\n", pedname);
    goto wdist_ret_0;
  }
  mapfile = fopen(mapname, "r");
  if (!mapfile) {
    printf("Error: Failed to open %s.\n", mapname);
    goto wdist_ret_0;
  }
  if (famname[0]) {
    binary_files = 1;
    famfile = fopen(famname, "r");
    if (!famfile) {
      printf("Error: Failed to open %s.\n", famname);
      goto wdist_ret_0;
    }
  }
  if (phenoname[0]) {
    phenofile = fopen(phenoname, "r");
    if (!phenofile) {
      printf("Error: Failed to open %s.\n", famname);
      goto wdist_ret_0;
    }
  }
  outname_end = outname;
  while (*outname_end) {
    outname_end++;
  }
  retval = RET_NOMEM;

  // note that this actually only allocates up to half of the main workspace to
  // the distance array, because of the difference between n(n-1)/2 and n^2.
  if (exp0) {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(int));
  } else {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(double));
  }

  // ----- .map/.bim load, first pass -----
  tbuf[MAXLINELEN - 6] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  while (fgets(tbuf, MAXLINELEN - 5, mapfile) != NULL) {
    if (tbuf[0] > ' ') {
      bufptr = next_item(tbuf);
      cptr = next_item(bufptr);
      if ((unsigned long)(cptr - bufptr) > max_marker_id_len) {
	max_marker_id_len = (unsigned long)(cptr - bufptr);
      }
      if (!unfiltered_marker_ct) {
	bufptr = next_item(cptr);
        if (binary_files) {
          bufptr = next_item(bufptr);
          bufptr = next_item(bufptr);
        }
	if (!bufptr) {
	  retval = RET_INVALID_FORMAT;
	  printf(errstr_map_format);
	  goto wdist_ret_1;
	}
	if (*bufptr > ' ') {
	  map_cols = 4;
	}
      }
      unfiltered_marker_ct += 1;
    }
    if (!tbuf[MAXLINELEN - 6]) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Excessively long line in .map/.bim file (max %d chars).\n", MAXLINELEN - 8);
      goto wdist_ret_1;
    }
  }
  if (!unfiltered_marker_ct) {
    retval = RET_INVALID_FORMAT;
    printf("Error: No markers in .map/.bim file.");
    goto wdist_ret_0;
  }
  rewind(mapfile);
  marker_exclude = (unsigned char*)calloc(sizeof(char), ((unfiltered_marker_ct + 7) / 8));
  if (!marker_exclude) {
    goto wdist_ret_1;
  }
  if (binary_files || freqname[0]) {
    marker_alleles = (char*)calloc(sizeof(char), unfiltered_marker_ct * 2);
    if (!marker_alleles) {
      goto wdist_ret_1;
    }
    marker_allele_cts = (int*)calloc(sizeof(int), unfiltered_marker_ct * 2);
    if (!marker_allele_cts) {
      goto wdist_ret_1;
    }
  }
  if (binary_files) {
    unfiltered_marker_ct4 = (unfiltered_marker_ct + 3) / 4;
    pedbuflen = unfiltered_marker_ct4;
  } else {
    pedbuflen = unfiltered_marker_ct * 4 + PEDBUFBASE;
  }
  pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
  if (!pedbuf) {
    goto wdist_ret_1;
  }
  marker_chrom = (int*)malloc(unfiltered_marker_ct * sizeof(int));
  if (!marker_chrom) {
    goto wdist_ret_1;
  }
  marker_id = (char*)malloc(unfiltered_marker_ct * max_marker_id_len);
  if (!marker_id) {
    goto wdist_ret_1;
  }
  mafs = (double*)malloc(unfiltered_marker_ct * sizeof(double));
  if (!mafs) {
    goto wdist_ret_1;
  }
  for (ii = 0; ii < unfiltered_marker_ct; ii++) {
    mafs[ii] = -1.0;
  }
  // marker_pos = (int*)malloc(unfiltered_marker_ct * sizeof(int));
  // if (!marker_pos) {
  //   goto wdist_ret_2;
  // }

  if (!(calculation_type & (CALC_DISTANCE_MASK | CALC_RELATIONSHIP_MASK | CALC_IBC | CALC_FREQ | CALC_WRITE_SNPLIST))) {
    prune = 1;
  }
  // ----- .map/.bim load, second pass -----
  for (ii = 0; ii < unfiltered_marker_ct; ii += 1) {
    do {
      fgets(tbuf, MAXLINELEN, mapfile);
    } while (tbuf[0] <= ' ');
    jj = marker_code(tbuf);
    if (jj == -1) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Invalid chromosome index in .map/.bim file.\n");
      goto wdist_ret_1;
    }
    marker_chrom[ii] = jj;
    if (!(chrom_mask & (1 << marker_code(tbuf)))) {
      exclude(marker_exclude, ii, &marker_exclude_ct);
    } else {
      bufptr = next_item(tbuf);
      cur_item(&(marker_id[ii * max_marker_id_len]), bufptr);
      bufptr = next_item(bufptr);
      if (map_cols == 4) {
	bufptr = next_item(bufptr);
      }
      if (*bufptr == '-') {
	if (binary_files) {
	  retval = RET_INVALID_FORMAT;
	  printf("Error: Negative marker position in .map/.bim file.\n");
	  goto wdist_ret_1;
	}
	exclude(marker_exclude, ii, &marker_exclude_ct);
      } else if (binary_files) {
        bufptr = next_item(bufptr);
        marker_alleles[ii * 2] = *bufptr;
        bufptr = next_item(bufptr);
        marker_alleles[ii * 2 + 1] = *bufptr;
      }
    }
    // marker_pos[ii] = atoi(bufptr);
  }
  if (marker_exclude_ct) {
    printf("%d markers loaded (after excluding %d).\n", unfiltered_marker_ct - marker_exclude_ct, marker_exclude_ct);
  } else {
    printf("%d markers loaded.\n", unfiltered_marker_ct - marker_exclude_ct);
  }

  // ----- phenotype file load, first pass -----
  if (phenofile) {
    if (makepheno_str) {
      if ((makepheno_str[0] == '*') && (makepheno_str[1] == '\0')) {
        makepheno_all = 1;
      } else {
        kk = strlen(makepheno_str);
      }
    }
    mm = tail_pheno; // is affection determined from file yet?
    nn = 0; // existence of header row determined from file yet?
    while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
      if (*tbuf == '\n') {
        continue;
      }
      bufptr = tbuf;
      ii = 0;
      jj = 0;
      while (!is_space_or_eoln(*bufptr)) {
        ii++;
        bufptr++;
      }
      if ((*bufptr == '\n') || (*bufptr == '\0')) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Improperly formatted phenotype file.\n");
        goto wdist_ret_1;
      }
      while ((*bufptr == ' ') || (*bufptr == '\t')) {
        bufptr++;
      }
      cptr = bufptr;
      while (!is_space_or_eoln(*bufptr)) {
        jj++;
        bufptr++;
      }
      if (ii + jj + 2 > max_pid_len) {
        max_pid_len = ii + jj + 2;
      }
      if (!nn) {
        nn = 1;
        if ((ii == 3) && (jj == 3) && (!memcmp("FID", tbuf, 3)) && (!memcmp("IID", cptr, 3))) {
          pheno_lines = -1;
          if (phenoname_str) {
            jj = strlen(phenoname_str);
            do {
              bufptr = next_item(bufptr);
              if ((!bufptr) || (*bufptr == '\n') || (*bufptr == '\0')) {
		retval = RET_INVALID_FORMAT;
		printf("Error: --pheno-name column not found.\n");
		goto wdist_ret_1;
              }
              mpheno_col++;
              ii = 0;
              cptr = bufptr;
              while (!is_space_or_eoln(*bufptr)) {
                bufptr++;
                ii++;
              }
	      if (ii != jj) {
		continue;
	      }
            } while (memcmp(cptr, phenoname_str, ii));
            phenoname_str = NULL;
          }
        } else if (phenoname_str) {
          retval = RET_INVALID_FORMAT;
          printf("Error: No header row in phenotype file.\n");
          goto wdist_ret_1;
        }
        if (!mpheno_col) {
          mpheno_col = 1;
        }
      }
      if (pheno_lines == -1) {
        pheno_lines = 0;
      } else {
        for (ii = 0; ii < mpheno_col; ii++) {
          bufptr = next_item(bufptr);
        }
        if (!bufptr || (*bufptr == '\n') || (*bufptr == '\0')) {
	  retval = RET_INVALID_FORMAT;
	  printf("Error: Improperly formatted phenotype file.\n");
	  goto wdist_ret_1;
        }
        if (makepheno_str && (!makepheno_all)) {
          if (!strncmp(makepheno_str, bufptr, kk)) {
            if (is_space_or_eoln(bufptr[kk])) {
              pheno_lines++;
            }
          }
        } else {
          if (!mm) {
            affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
            mm = 1;
          }
          if (affection) {
            if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
              pheno_lines++;
            }
          } else {
	    if (sscanf(bufptr, "%lg", &dxx) != 1) {
	      retval = RET_INVALID_FORMAT;
	      printf("Error: Improperly formatted phenotype file.\n");
	      goto wdist_ret_1;
	    }
            if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
              pheno_lines++;
            }
          }
        }
      }
      if (!tbuf[MAXLINELEN - 1]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Excessively long line in phenotype file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_1;
      }
    }
    id_buf = (char*)malloc(max_pid_len * sizeof(char));
    if (!id_buf) {
      goto wdist_ret_1;
    }
    if (pheno_lines) {
      pid_list = (char*)malloc(max_pid_len * pheno_lines * sizeof(char));
      if (!pid_list) {
	goto wdist_ret_1;
      }
      if (affection || tail_pheno) {
	phenor_c = (char*)malloc(pheno_lines * sizeof(char));
	if (!phenor_c) {
	  goto wdist_ret_1;
	}
      } else if (!makepheno_str) {
	phenor_d = (double*)malloc(pheno_lines * sizeof(double));
	if (!phenor_d) {
	  goto wdist_ret_1;
	}
      }

      rewind(phenofile);
      ii = 0;
      // ----- phenotype file load, second pass -----
      while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
	if (ii == pheno_lines) {
	  break;
	}
	if (*tbuf == '\n') {
	  continue;
	}
	cptr = &(pid_list[ii * max_pid_len]);
	bufptr = tbuf;
	while ((*bufptr != ' ') && (*bufptr != '\t')) {
	  *cptr++ = *bufptr++;
	}
	*cptr++ = ' ';
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	while (!is_space_or_eoln(*bufptr)) {
	  *cptr++ = *bufptr++;
	}
	for (jj = 0; jj < mpheno_col; jj++) {
	  bufptr = next_item(bufptr);
	}
	*cptr = '\0';
	if (makepheno_str) {
	  if (makepheno_all) {
	    ii++;
	  } else {
	    if (!strncmp(makepheno_str, bufptr, kk)) {
	      if (is_space_or_eoln(bufptr[kk])) {
		ii++;
	      }
	    }
	  }
	} else if (affection) {
	  if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
	    cc = -1;
	  } else if (affection_01) {
	    cc = *bufptr - '0';
	  } else {
	    cc = *bufptr - '1';
	  }
	  if ((cc != -1) || (!prune)) {
	    phenor_c[ii++] = cc;
	  }
	} else {
	  sscanf(bufptr, "%lg", &dxx);
	  if (!(prune && (dxx == missing_phenod))) {
	    if (tail_pheno) {
	      if (dxx == missing_phenod) {
		phenor_c[ii++] = -1;
	      } else if (dxx <= tail_bottom) {
		phenor_c[ii++] = 0;
	      } else if (dxx > tail_top) {
		phenor_c[ii++] = 1;
	      } else if (!prune) {
		phenor_c[ii++] = -1;
	      }
	    } else {
	      phenor_d[ii++] = dxx;
	    }
	  }
	}
      }
      if (affection || tail_pheno) {
	if (qsort_ext(pid_list, pheno_lines, max_pid_len, strcmp_deref, phenor_c, sizeof(char)) == -1) {
	  goto wdist_ret_1;
        }
      } else if (makepheno_str) {
        qsort(pid_list, pheno_lines, max_pid_len, strcmp_casted);
      } else {
	if (qsort_ext(pid_list, pheno_lines, max_pid_len, strcmp_deref, (char*)phenor_d, sizeof(double)) == -1) {
          goto wdist_ret_1;
        }
      }
    } else {
      printf("Note: No valid entries in phenotype file.\n");
    }
    fclose(phenofile);
    phenofile = NULL;
  }

  if (filter_type) {
    filterfile = fopen(filtername, "r");
    if (!filterfile) {
      printf("Error: Failed to open %s.\n", filtername);
      goto wdist_ret_1;
    }
  }

  if (extractname[0] || excludename[0]) {
    wkspace_mark = wkspace_base;
    ii = unfiltered_marker_ct - marker_exclude_ct;
    retval = marker_id_sort(&cptr, &iptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_id, max_marker_id_len);
    if (extractname[0]) {
      retval = extract_exclude_markers(extractname, cptr, ii, max_marker_id_len, iptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 0);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    if (excludename[0]) {
      retval = extract_exclude_markers(excludename, cptr, ii, max_marker_id_len, iptr, unfiltered_marker_ct, marker_exclude, &marker_exclude_ct, 1);
      if (retval) {
	goto wdist_ret_1;
      }
    }
    wkspace_reset(wkspace_mark);
  }
  retval = RET_NOMEM;

  // ----- filter load, first pass -----
  if (filter_type) {
    if (filter_type == FILTER_CUSTOM) {
      jj = strlen(filterval);
      if (!mfilter_col) {
        mfilter_col = 1;
      }
    }
    while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
      if (*tbuf == '\n') {
        continue;
      }
      bufptr = tbuf;
      ii = 0;
      while (!is_space_or_eoln(*bufptr)) {
        ii++;
        bufptr++;
      }
      if ((*bufptr == '\n') || (*bufptr == '\0')) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Improperly formatted filter file.\n");
        goto wdist_ret_1;
      }
      while ((*bufptr == ' ') || (*bufptr == '\t')) {
        bufptr++;
      }
      while (!is_space_or_eoln(*bufptr)) {
        ii++;
        bufptr++;
      }
      ii += 2;
      if (ii > max_id_len) {
        max_id_len = ii;
      }
      if (filter_type == FILTER_CUSTOM) {
        for (kk = 0; kk < mfilter_col; kk++) {
          bufptr = next_item(bufptr);
        }
	if ((!bufptr) || (*bufptr == '\n') || (*bufptr == '\0')) {
	  retval = RET_INVALID_FORMAT;
	  printf("Error: Improperly formatted filter file.\n");
	  goto wdist_ret_1;
        }
        if (!strncmp(filterval, bufptr, jj)) {
          filter_lines += 1;
        }
      } else {
        filter_lines += 1;
      }
      if (!tbuf[MAXLINELEN - 1]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Excessively long line in filter file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_1;
      }
    }
    if (max_pid_len && (max_id_len > max_pid_len)) {
      free(id_buf);
    }
    id_buf = (char*)malloc(max_id_len * sizeof(char));
    if (!id_buf) {
      goto wdist_ret_1;
    }
    if (filter_lines) {
      id_list = (char*)malloc(max_id_len * filter_lines * sizeof(char));
      if (!id_list) {
	goto wdist_ret_1;
      }
      rewind(filterfile);
      ii = 0;

      // ----- filter load, second pass -----
      while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
	if (ii == filter_lines) {
	  break;
	}
	if (*tbuf == '\n') {
	  continue;
	}
	cptr = &(id_list[ii * max_id_len]);
	bufptr = tbuf;
	while ((*bufptr != ' ') && (*bufptr != '\t')) {
	  *cptr++ = *bufptr++;
	}
	*cptr++ = ' ';
	while ((*bufptr == ' ') || (*bufptr == '\t')) {
	  bufptr++;
	}
	while (!is_space_or_eoln(*bufptr)) {
	  *cptr++ = *bufptr++;
	}
	if (filter_type == FILTER_CUSTOM) {
	  for (kk = 0; kk < mfilter_col; kk++) {
	    bufptr = next_item(bufptr);
	  }
	  if (!strncmp(filterval, bufptr, jj)) {
	    ii++;
	    *cptr = '\0';
	  }
	} else {
	  ii++;
	  *cptr = '\0';
	}
      }
      qsort(id_list, filter_lines, max_id_len, strcmp_casted);
    } else {
      printf("Note: No valid entries in filter file.\n");
    }
    fclose(filterfile);
    filterfile = NULL;
  }
  line_locs = (long long*)malloc(max_people * sizeof(long long));
  if (!line_locs) {
    goto wdist_ret_1;
  }

  mind_int_thresh = (int)(2 * mind_thresh * (unfiltered_marker_ct - marker_exclude_ct));
  if (binary_files) {
    nn = 0; // number of people that pass initial filter
    // ----- .fam load, first pass -----
    while (fgets(tbuf, MAXLINELEN, famfile) != NULL) {
      if (tbuf[0] > ' ') {
        if (tbuf[0] != '#') {
	  bufptr = next_item(tbuf);
          if (!ped_col_1) {
            cptr = tbuf;
          } else {
            cptr = bufptr;
          }
          ii = strlen_se(tbuf) + strlen_se(cptr) + 2;
          if (ii > max_person_id_len) {
            max_person_id_len = ii;
          }
          if (filter_type) {
            ii = is_contained(id_buf, id_list, max_id_len, filter_lines, tbuf, cptr);
            if (filter_type == FILTER_REMOVE) {
              ii = 1 - ii;
            }
          } else {
            ii = 1;
          }
          if (phenoname[0]) {
            if (ii) {
	      if (makepheno_str || (!prune)) {
		line_locs[nn++] = unfiltered_indiv_ct;
	      } else if (is_contained(id_buf, pid_list, max_pid_len, pheno_lines, tbuf, cptr)) {
		line_locs[nn++] = unfiltered_indiv_ct;
	      }
            }
          } else if (ii || (!unfiltered_indiv_ct)) {
	    if (ped_col_1) {
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_34) {
	      bufptr = next_item(bufptr);
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_5) {
	      bufptr = next_item(bufptr);
	    }
	    if (!bufptr) {
	      retval = RET_INVALID_FORMAT;
	      printf(errstr_fam_format);
	      goto wdist_ret_1;
	    }
            if (!unfiltered_indiv_ct) {
              affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
            }
            if (ii) {
	      if (affection) {
                if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
                  line_locs[nn++] = unfiltered_indiv_ct;
                }
	      } else {
                if (sscanf(bufptr, "%lg", &dxx) != 1) {
                  retval = RET_INVALID_FORMAT;
                  printf(errstr_fam_format);
                  goto wdist_ret_1;
                }
                if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
                  line_locs[nn++] = unfiltered_indiv_ct;
                }
	      }
            }
          }
          unfiltered_indiv_ct++;
        }
      }
      if (!tbuf[MAXLINELEN - 1]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Excessively long line in .fam file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_1;
      }
    }
    if (nn < 2) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Less than two valid people in .fam file after initial filter.\n");
      goto wdist_ret_1;
    } else if (nn > max_people) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Too many people in .fam file for this algorithm.\n");
      goto wdist_ret_1;
    }

    bin_pheno = (makepheno_str || affection || tail_pheno);
    if (bin_pheno) {
      pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
      if (!pheno_c) {
	goto wdist_ret_1;
      }
    } else {
      pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
      if (!pheno_d) {
	goto wdist_ret_1;
      }
    }
    person_id = (char*)malloc(unfiltered_indiv_ct * max_person_id_len * sizeof(char));
    if (!person_id) {
      goto wdist_ret_1;
    }

    indiv_exclude = (unsigned char*)calloc(sizeof(char), ((unfiltered_indiv_ct + 7) / 8));
    if (!indiv_exclude) {
      goto wdist_ret_1;
    }

    rewind(famfile);
    // ----- .fam load, second pass -----
    ii = 0; // raw line number
    jj = 0; // loaded lines
    while (fgets(tbuf, MAXLINELEN, famfile) != NULL) {
      if (tbuf[0] > ' ') {
        if (tbuf[0] != '#') {
	  if (jj == nn) {
            while (ii < unfiltered_indiv_ct) {
              exclude(indiv_exclude, ii++, &indiv_exclude_ct);
            }
	    break;
	  }
	  if (line_locs[jj] > ii) {
            exclude(indiv_exclude, ii++, &indiv_exclude_ct);
	    continue;
	  }

	  bufptr = next_item(tbuf);
	  if (!ped_col_1) {
	    cptr = tbuf;
	  } else {
	    cptr = bufptr;
	  }
	  kk = strlen_se(tbuf);
	  mm = strlen_se(cptr);
	  memcpy(&(person_id[ii * max_person_id_len]), tbuf, kk);
	  person_id[ii * max_person_id_len + kk] = '\t';
	  memcpy(&(person_id[ii * max_person_id_len + kk + 1]), cptr, mm);
	  person_id[ii * max_person_id_len + kk + mm + 1] = '\0';
          if (phenoname[0]) {
	    kk = bsearch_fam_indiv(id_buf, pid_list, max_pid_len, pheno_lines, tbuf, cptr);
            if (makepheno_str) {
              if (kk == -1) {
                pheno_c[ii] = 0;
              } else {
                pheno_c[ii] = 1;
              }
            } else if (affection || tail_pheno) {
              if (kk == -1) {
                pheno_c[ii] = -1;
              } else {
                pheno_c[ii] = phenor_c[kk];
              }
            } else if (kk == -1) {
              pheno_d[ii] = missing_phenod;
            } else {
              pheno_d[ii] = phenor_d[kk];
            }
          } else {
	    if (ped_col_1) {
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_34) {
	      bufptr = next_item(bufptr);
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_5) {
	      bufptr = next_item(bufptr);
	    }
            if (affection) {
              if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
                pheno_c[ii] = -1;
              } else if (affection_01) {
                pheno_c[ii] = *bufptr - '0';
              } else {
                pheno_c[ii] = *bufptr - '1';
              }
            } else {
              sscanf(bufptr, "%lg", &dxx);
              if (tail_pheno) {
                if (dxx == missing_phenod) {
                  pheno_c[ii] = -1;
                } if (dxx <= tail_bottom) {
                  pheno_c[ii] = 0;
                } else if (dxx > tail_top) {
                  pheno_c[ii] = 1;
                } else {
                  pheno_c[ii] = -1;
                }
              } else {
                pheno_d[ii] = dxx;
              }
            }
          }
	  ii++;
	  jj++;
        }
      }
    }

    if (freqname[0]) {
      retval = update_freq(freqname, &freqfile, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_id, max_marker_id_len, marker_chrom, marker_alleles, mafs);
      if (retval) {
	goto wdist_ret_1;
      }
    }

    maf_int_thresh = 2 * (nn - indiv_exclude_ct) - (int)((1.0 - min_maf) * (nn - indiv_exclude_ct) * 2);
    geno_int_thresh = 2 * (nn - indiv_exclude_ct) - (int)(geno_thresh * 2 * (nn - indiv_exclude_ct));

    if (fread(pedbuf, 1, 3, pedfile) < 3) {
      retval = RET_READ_FAIL;
      goto wdist_ret_1;
    }
    if ((pedbuf[0] != 'l') || (pedbuf[1] != '\x1b')) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Invalid or pre-v1.00 BED file.\n");
      goto wdist_ret_1;
    }

    if (pedbuf[2] == '\0') {
      // ----- individual-major .bed load, first pass -----
      for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
        if (excluded(indiv_exclude, ii)) {
          fseeko(pedfile, (off_t)unfiltered_marker_ct4, SEEK_CUR);
          continue;
        }
        if (fread(pedbuf, 1, unfiltered_marker_ct4, pedfile) < unfiltered_marker_ct4) {
          retval = RET_READ_FAIL;
          goto wdist_ret_1;
        }
        gptr = pedbuf - 1;
        mm = 0; // missing
        for (jj = 0; jj < unfiltered_marker_ct; jj++) {
          if (excluded(marker_exclude, jj)) {
            continue;
	  }
          kk = jj % 4;
          if (!kk) {
            nn = *(++gptr);
          } else {
            nn >>= 2;
          }
          oo = nn & 3;
          if (oo) {
            if (oo == 2) {
              marker_allele_cts[jj * 2] += 1;
              marker_allele_cts[jj * 2 + 1] += 1;
            } else if (oo == 3) {
              marker_allele_cts[jj * 2] += 2;
            } else {
              mm++;
            }
          } else {
            marker_allele_cts[jj * 2 + 1] += 2;
          }
        }
        if (mm > mind_int_thresh) {
          exclude(indiv_exclude, ii, &indiv_exclude_ct);
        }
      }
      for (ii = 0; ii < unfiltered_marker_ct; ii++) {
        if ((marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1] < geno_int_thresh) || ((mafs[ii] < 0.0) && ((marker_allele_cts[ii * 2 + 1] < maf_int_thresh) || (marker_allele_cts[ii * 2] < maf_int_thresh))) || ((mafs[ii] >= 0.0) && ((mafs[ii] < min_maf) || (mafs[ii] > 1.0 - min_maf)))) {
          exclude(marker_exclude, ii, &marker_exclude_ct);
        }
      }
      het_probs = (double*)malloc(((unfiltered_marker_ct - marker_exclude_ct) / 2 + 1) * sizeof(double));
      if (!het_probs) {
	goto wdist_ret_1;
      }
      if (wkspace_left < unfiltered_marker_ct * 9 * sizeof(int)) {
	goto wdist_ret_1;
      }
      // ----- individual-major .bed load, second pass -----
      fseeko(pedfile, 3, SEEK_SET);
      hwe_ll = (int*)wkspace;
      hwe_lh = (int*)(&wkspace[unfiltered_marker_ct * sizeof(int)]);
      hwe_hh = (int*)(&wkspace[unfiltered_marker_ct * 2 * sizeof(int)]);
      hwe_u_ll = (int*)(&wkspace[unfiltered_marker_ct * 3 * sizeof(int)]);
      hwe_u_lh = (int*)(&wkspace[unfiltered_marker_ct * 4 * sizeof(int)]);
      hwe_u_hh = (int*)(&wkspace[unfiltered_marker_ct * 5 * sizeof(int)]);
      hwe_a_ll = (int*)(&wkspace[unfiltered_marker_ct * 6 * sizeof(int)]);
      hwe_a_lh = (int*)(&wkspace[unfiltered_marker_ct * 7 * sizeof(int)]);
      hwe_a_hh = (int*)(&wkspace[unfiltered_marker_ct * 8 * sizeof(int)]);
      fill_int_zero((int*)wkspace, unfiltered_marker_ct * 9);
      for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
	if (excluded(indiv_exclude, ii)) {
	  fseeko(pedfile, (off_t)unfiltered_marker_ct4, SEEK_CUR);
	  continue;
	}
	if (fread(pedbuf, 1, unfiltered_marker_ct4, pedfile) < unfiltered_marker_ct4) {
	  retval = RET_READ_FAIL;
	  goto wdist_ret_1;
	}
	gptr = pedbuf - 1;
	for (jj = 0; jj < unfiltered_marker_ct; jj++) {
	  kk = jj % 4;
	  if (!kk) {
	    nn = *(++gptr);
	  } else {
	    nn >>= 2;
	  }
	  oo = nn & 3;
	  if (oo) {
	    if (oo == 2) {
	      hwe_lh[jj] += 1;
	      if (bin_pheno) {
		if (pheno_c[ii] == 0) {
		  hwe_u_lh[jj] += 1;
		} else if (pheno_c[ii] == 1) {
		  hwe_a_lh[jj] += 1;
		}
	      }
	    } else if (oo == 3) {
	      hwe_hh[jj] += 1;
	      if (bin_pheno) {
		if (pheno_c[ii] == 0) {
		  hwe_u_hh[jj] += 1;
		} else if (pheno_c[ii] == 1) {
		  hwe_a_hh[jj] += 1;
		}
	      }
	    }
	  } else {
	    hwe_ll[jj] += 1;
	    if (bin_pheno) {
	      if (pheno_c[ii] == 0) {
		hwe_u_ll[jj] += 1;
	      } else if (pheno_c[ii] == 1) {
		hwe_a_ll[jj] += 1;
	      }
	    }
	  }
	}
      }
      marker_weights = (double*)wkspace_alloc((unfiltered_marker_ct - marker_exclude_ct) * sizeof(double));
      if (!marker_weights) {
	goto wdist_ret_2;
      }
      dptr2 = marker_weights;

      for (ii = 0; ii < unfiltered_marker_ct; ii++) {
	if (!excluded(marker_exclude, ii)) {
	  if (SNPHWE_t(hwe_lh[ii], hwe_ll[ii], hwe_hh[ii], hwe_thresh)) {
	    exclude(marker_exclude, ii, &marker_exclude_ct);
	  } else if (bin_pheno) {
	    if (SNPHWE_t(hwe_u_lh[ii], hwe_u_ll[ii], hwe_u_hh[ii], hwe_thresh) || SNPHWE_t(hwe_a_lh[ii], hwe_a_ll[ii], hwe_a_hh[ii], hwe_thresh)) {
	      exclude(marker_exclude, ii, &marker_exclude_ct);
	    }
	  }
	  if (!excluded(marker_exclude, ii)) {
            marker_allele_cts[ii * 2] = hwe_hh[ii] * 2 + hwe_lh[ii];
            marker_allele_cts[ii * 2 + 1] = hwe_ll[ii] * 2 + hwe_lh[ii];
            if (mafs[ii] < 0.0) {
	      *dptr2++ = calc_wt_mean(exponent, hwe_lh[ii], hwe_ll[ii], hwe_hh[ii]);
	      mafs[ii] = ((double)marker_allele_cts[ii * 2]) / ((double)(marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]));
            } else {
              *dptr2++ = calc_wt_mean_maf(exponent, mafs[ii]);
            }
	  }
	}
      }
      free(het_probs);
      het_probs = NULL;
    } else {
      // ----- snp-major .bed load, first pass -----
      snp_major = 1;
      person_missing_cts = (int*)calloc(sizeof(int), unfiltered_indiv_ct);
      if (!person_missing_cts) {
        goto wdist_ret_1;
      }
      unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
      if (pedbuflen < (4 * unfiltered_indiv_ct4)) {
        free(pedbuf);
        pedbuflen = 4 * unfiltered_indiv_ct4;
        pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
        if (!pedbuf) {
          goto wdist_ret_2;
        }
      }
      for (ii = 0; ii < unfiltered_marker_ct; ii++) {
        if (excluded(marker_exclude, ii)) {
          fseeko(pedfile, (off_t)unfiltered_indiv_ct4, SEEK_CUR);
          continue;
        }
        if (fread(pedbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
          retval = RET_READ_FAIL;
          goto wdist_ret_2;
        }
        mm = 0; // major allele ct
        nn = 0; // minor allele ct
        gptr = pedbuf - 1;
        for (jj = 0; jj < unfiltered_indiv_ct; jj++) {
          if (excluded(indiv_exclude, jj)) {
            continue;
          }
          kk = jj % 4;
          if (!kk) {
            uii = *(++gptr);
          } else {
            uii >>= 2;
          }
          pp = uii & 3;
          if (pp) {
            if (pp == 2) {
              mm++;
              nn++;
            } else if (pp == 3) {
              mm += 2;
            } else {
              person_missing_cts[jj] += 1;
            }
          } else {
            nn += 2;
          }
        }
        if ((mm + nn < geno_int_thresh) || ((mafs[ii] < 0.0) && ((mm < maf_int_thresh) || (nn < maf_int_thresh))) || ((mafs[ii] >= 0.0) && ((mafs[ii] < min_maf) || (mafs[ii] > 1.0 - min_maf)))) {
          exclude(marker_exclude, ii, &marker_exclude_ct);
        }
      }
      for (ii = 0; ii < unfiltered_indiv_ct; ii++) {
        if (person_missing_cts[ii] > mind_int_thresh) {
	  exclude(indiv_exclude, ii, &indiv_exclude_ct);
        }
      }

      het_probs = (double*)malloc(((unfiltered_marker_ct - marker_exclude_ct) / 2 + 1) * sizeof(double));
      if (!het_probs) {
	goto wdist_ret_1;
      }
      marker_weights = (double*)wkspace_alloc((unfiltered_marker_ct - marker_exclude_ct) * sizeof(double));
      if (!marker_weights) {
	goto wdist_ret_2;
      }
      dptr2 = marker_weights;

      fseeko(pedfile, 3, SEEK_SET);
      // ----- .snp-major .bed load, second pass -----
      for (ii = 0; ii < unfiltered_marker_ct; ii++) {
	if (excluded(marker_exclude, ii)) {
	  fseeko(pedfile, (off_t)unfiltered_indiv_ct4, SEEK_CUR);
	  continue;
	}
	if (fread(pedbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	  retval = RET_READ_FAIL;
	  goto wdist_ret_2;
	}
	gptr = pedbuf;
	hwe_lli = 0;
	hwe_lhi = 0;
	hwe_hhi = 0;
	hwe_u_lli = 0;
	hwe_u_lhi = 0;
	hwe_u_hhi = 0;
	hwe_a_lli = 0;
	hwe_a_lhi = 0;
	hwe_a_hhi = 0;
	for (jj = 0; jj < unfiltered_indiv_ct; jj++) {
	  kk = jj % 4;
	  if (!kk) {
	    uii = *gptr++;
	  } else {
	    uii >>= 2;
	  }
	  if (!excluded(indiv_exclude, jj)) {
	    ujj = uii & 3;
	    if (ujj) {
	      if (ujj == 2) {
		hwe_lhi += 1;
		if (bin_pheno) {
		  if (pheno_c[jj] == 0) {
		    hwe_u_lhi += 1;
		  } else if (pheno_c[jj] == 1) {
		    hwe_a_lhi += 1;
		  }
		}
	      } else if (ujj == 3) {
		hwe_hhi += 1;
		if (bin_pheno) {
		  if (pheno_c[jj] == 0) {
		    hwe_u_hhi += 1;
		  } else if (pheno_c[jj] == 1) {
		    hwe_a_hhi += 1;
		  }
		}
	      }
	    } else {
	      hwe_lli += 1;
	      if (bin_pheno) {
		if (pheno_c[jj] == 0) {
		  hwe_u_lli += 1;
		} else if (pheno_c[jj] == 1) {
		  hwe_a_lli += 1;
		}
	      }
	    }
	  }
	}
	kk = SNPHWE_t(hwe_lhi, hwe_lli, hwe_hhi, hwe_thresh);
	if (bin_pheno && (!kk)) {
	  kk = SNPHWE_t(hwe_u_lhi, hwe_u_lli, hwe_u_hhi, hwe_thresh) || SNPHWE_t(hwe_a_lhi, hwe_a_lli, hwe_a_hhi, hwe_thresh);
	}
	if (kk) {
	  exclude(marker_exclude, ii, &marker_exclude_ct);
	} else {
          marker_allele_cts[ii * 2] = hwe_hhi * 2 + hwe_lhi;
          marker_allele_cts[ii * 2 + 1] = hwe_lli * 2 + hwe_lhi;
          if (mafs[ii] < 0.0) {
	    *dptr2++ = calc_wt_mean(exponent, hwe_lhi, hwe_lli, hwe_hhi);
	    mafs[ii] = ((double)marker_allele_cts[ii * 2]) / ((double)(marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]));
          } else {
            *dptr2++ = calc_wt_mean_maf(exponent, mafs[ii]);
          }
	}
      }
      free(het_probs);
      het_probs = NULL;
    }
  } else {
    if (freqname[0]) {
      printf("--update-freq currently ignored on text file load.\n");
    }
    pedbuf[pedbuflen - 1] = ' ';
    nn = 0; // number of individuals after first-level filter
    // define exclude counters
    // ----- .ped load, first pass -----
    while (fgets((char*)pedbuf, pedbuflen, pedfile) != NULL) {
      if ((pedbuf[0] > ' ') && (pedbuf[0] != '#')) {
	bufptr = next_item((char*)pedbuf);
	if (!ped_col_1) {
	  cptr = (char*)pedbuf;
	} else {
	  cptr = bufptr;
	}
	ii = strlen_se(tbuf) + strlen_se(cptr) + 2;
	if (ii > max_person_id_len) {
	  max_person_id_len = ii;
	}
	if (filter_type) {
	  ii = is_contained(id_buf, id_list, max_id_len, filter_lines, (char*)pedbuf, cptr);
	  if (filter_type == FILTER_REMOVE) {
	    ii = 1 - ii;
	  }
	} else {
	  ii = 1;
	}
	if (phenoname[0]) {
	  if (ii) {
	    if (makepheno_str || (!prune)) {
	      line_locs[nn++] = unfiltered_indiv_ct;
	    } else {
	      if (is_contained(id_buf, pid_list, max_pid_len, pheno_lines, (char*)pedbuf, cptr)) {
		line_locs[nn++] = unfiltered_indiv_ct;
	      }
	    }
	  }
	} else if (ii || (!unfiltered_indiv_ct)) {
	  if (ped_col_1) { // actually column 2
	    bufptr = next_item(bufptr);
	  }
	  if (ped_col_34) {
	    bufptr = next_item(bufptr);
	    bufptr = next_item(bufptr);
	  }
	  if (ped_col_5) {
	    bufptr = next_item(bufptr);
	  }
	  if (!bufptr) {
	    retval = RET_INVALID_FORMAT;
	    printf(errstr_ped_format);
	    goto wdist_ret_1;
	  }
	  if (!unfiltered_indiv_ct) {
	    affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
	  }
	  if (affection) {
	    if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
	      line_locs[nn++] = unfiltered_indiv_ct;
	    }
	  } else {
	    if (sscanf(bufptr, "%lg", &dxx) != 1) {
	      retval = RET_INVALID_FORMAT;
	      printf(errstr_ped_format);
	      goto wdist_ret_1;
	    }
	    if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
	      line_locs[nn++] = unfiltered_indiv_ct;
	    }
          }
	}
        unfiltered_indiv_ct++;
      }
      if (!pedbuf[pedbuflen - 1]) {
        pedbuf[pedbuflen - 1] = ' ';
        if (pedbuf[pedbuflen - 2] == '\n') {
	  last_tell = ftello(pedfile);
          continue;
        }
	ii = 0;
	do {
	  ii += pedbuflen - 1;
	  pedbuf[pedbuflen - 1] = ' ';
	  fgets_return = fgets((char*)pedbuf, pedbuflen, pedfile);
	} while (fgets_return && !pedbuf[pedbuflen - 1] && (pedbuf[pedbuflen - 2] != '\n'));
	ii += strlen((char*)pedbuf) + 1;
	if (ii > ped_recalc_len) {
	  ped_recalc_len = ii;
	}
      }
    }
    if (nn < 2) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Less than two valid people in .ped file.\n");
      goto wdist_ret_1;
    } else if (unfiltered_indiv_ct > max_people) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Too many people in .ped file for this algorithm.\n");
      goto wdist_ret_1;
    }
    rewind(pedfile);
    if (ped_recalc_len) {
      free(pedbuf);
      pedbuflen = ped_recalc_len;
      pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
      if (!pedbuf) {
	goto wdist_ret_1;
      }
    }

    person_id = (char*)malloc(unfiltered_indiv_ct * max_person_id_len * sizeof(char));
    if (!person_id) {
      goto wdist_ret_1;
    }

    maf_int_thresh = 2 * unfiltered_indiv_ct - (int)((1.0 - min_maf) * unfiltered_indiv_ct * 2);
    geno_int_thresh = 2 * unfiltered_indiv_ct - (int)(geno_thresh * 2 * unfiltered_indiv_ct);

    marker_alleles = (char*)calloc(sizeof(char), unfiltered_marker_ct * 4);
    if (!marker_alleles) {
      goto wdist_ret_1;
    }
    marker_allele_cts = (int*)calloc(sizeof(int), unfiltered_marker_ct * 4);
    if (!marker_allele_cts) {
      goto wdist_ret_1;
    }
    bin_pheno = (makepheno_str || affection || tail_pheno);
    if (bin_pheno) {
      pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
      if (!pheno_c) {
	goto wdist_ret_1;
      }
    } else {
      pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
      if (!pheno_d) {
	goto wdist_ret_1;
      }
    }
    indiv_exclude = (unsigned char*)calloc(sizeof(char), ((unfiltered_indiv_ct + 7) / 8));
    if (!indiv_exclude) {
      goto wdist_ret_1;
    }

    // ----- .ped load, second pass -----
    kk = 0;
    for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
      do {
        last_tell = ftello(pedfile);
        fgets((char*)pedbuf, pedbuflen, pedfile);
      } while ((pedbuf[0] <= ' ') || (pedbuf[0] == '#'));
      if (line_locs[kk] > ii) {
        exclude(indiv_exclude, ii, &indiv_exclude_ct);
        continue;
      }
      bufptr = (char*)pedbuf;
      oo = strlen_se(bufptr);
      if (ped_col_1) {
	bufptr = next_item(bufptr);
      }
      mm = strlen_se(bufptr);
      memcpy(&(person_id[ii * max_person_id_len]), pedbuf, oo);
      person_id[ii * max_person_id_len + oo] = '\t';
      memcpy(&(person_id[ii * max_person_id_len + oo + 1]), bufptr, mm);
      person_id[ii * max_person_id_len + oo + mm + 1] = '\0';
      if (phenoname[0]) {
        jj = bsearch_fam_indiv(id_buf, pid_list, max_pid_len, pheno_lines, (char*)pedbuf, bufptr);
        cc = 0;
        if (makepheno_str) {
          if (jj != -1) {
            cc = 1;
          }
          pheno_c[ii] = cc;
        } else if (affection || tail_pheno) {
          if (jj == -1) {
            pheno_c[ii] = -1;
          } else {
            pheno_c[ii] = phenor_c[jj];
          }
        } else if (jj == -1) {
          pheno_d[ii] = missing_phenod;
        } else {
          pheno_d[ii] = phenor_d[jj];
        }
        bufptr = next_item(bufptr);
        if (ped_col_34) {
          bufptr = next_item(bufptr);
          bufptr = next_item(bufptr);
        }
        if (ped_col_5) {
          bufptr = next_item(bufptr);
        }
      } else {
	if (ped_col_34) {
	  bufptr = next_item(bufptr);
	  bufptr = next_item(bufptr);
	}
	if (ped_col_5) {
	  bufptr = next_item(bufptr);
	}
	if (affection) {
          if (is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01)) {
            pheno_c[ii] = -1;
          } else if (affection_01) {
            pheno_c[ii] = *bufptr - '0';
	  } else {
            pheno_c[ii] = *bufptr - '1';
	  }
	} else {
	  sscanf(bufptr, "%lg", &dxx);
          if (tail_pheno) {
            if (dxx == missing_phenod) {
              pheno_c[ii] = -1;
	    } else if (dxx <= tail_bottom) {
              pheno_c[ii] = 0;
            } else if (dxx > tail_top) {
              pheno_c[ii] = 1;
            } else {
              pheno_c[ii] = -1;
            }
          } else {
            pheno_d[ii] = dxx;
          }
	}
      }
      if (ped_col_6) {
        bufptr = next_item(bufptr);
      }
      if (ped_col_7) {
	bufptr = next_item(bufptr);
      }
      if (!bufptr) {
        retval = RET_INVALID_FORMAT;
        printf(errstr_ped_format);
        goto wdist_ret_1;
      }
      line_locs[kk] = last_tell + (bufptr - (char*)pedbuf);
      mm = 0; // number of missing
      for (jj = 0; jj < unfiltered_marker_ct; jj += 1) {
        for (mm = 0; mm < 2; mm++) {
          cc = *bufptr;
	  if (cc == '\0') {
	    retval = RET_INVALID_FORMAT;
	    printf(errstr_ped_format);
	    goto wdist_ret_1;
	  }
          if (!excluded(marker_exclude, jj)) {
	    if (cc == missing_geno) {
	      mm += 1;
	    } else if (cc == marker_alleles[jj * 4]) {
	      marker_allele_cts[jj * 4] += 1;
	    } else if (cc == marker_alleles[jj * 4 + 1]) {
	      marker_allele_cts[jj * 4 + 1] += 1;
	    } else if (cc == marker_alleles[jj * 4 + 2]) {
	      marker_allele_cts[jj * 4 + 2] += 1;
	    } else if (cc == marker_alleles[jj * 4 + 3]) {
	      marker_allele_cts[jj * 4 + 3] += 1;
	    } else if (marker_alleles[jj * 4 + 3]) {
	      retval = RET_INVALID_FORMAT;
	      printf("Error: More than 4 different allele types at marker %d.\n", jj + 1);
	      goto wdist_ret_1;
	    } else if (marker_alleles[jj * 4 + 2]) {
	      marker_alleles[jj * 4 + 3] = cc;
	      marker_allele_cts[jj * 4 + 3] = 1;
	    } else if (marker_alleles[jj * 4 + 1]) {
	      marker_alleles[jj * 4 + 2] = cc;
	      marker_allele_cts[jj * 4 + 2] = 1;
	    } else if (marker_alleles[jj * 4]) {
	      marker_alleles[jj * 4 + 1] = cc;
	      marker_allele_cts[jj * 4 + 1] = 1;
	    } else {
	      marker_alleles[jj * 4] = cc;
	      marker_allele_cts[jj * 4] = 1;
	    }
          }
	  bufptr++;
	  while ((*bufptr == ' ') || (*bufptr == '\t')) {
	    bufptr++;
	  }
        }
      }
      if (mm > mind_int_thresh) {
        exclude(indiv_exclude, ii, &indiv_exclude_ct);
      }
      kk++;
    }
    rewind(pedfile);
    for (ii = 0; ii < unfiltered_marker_ct; ii += 1) {
      for (jj = 1; jj < 3; jj += 1) {
        for (kk = jj - 1; kk >= 0; kk -= 1) {
          if (marker_allele_cts[ii * 4 + kk] < marker_allele_cts[ii * 4 + kk + 1]) {
            cc = marker_alleles[ii * 4 + kk];
            marker_alleles[ii * 4 + kk] = marker_alleles[ii * 4 + kk + 1];
            marker_alleles[ii * 4 + kk + 1] = cc;
            mm = marker_allele_cts[ii * 4 + kk];
            marker_allele_cts[ii * 4 + kk] = marker_allele_cts[ii * 4 + kk + 1];
            marker_allele_cts[ii * 4 + kk + 1] = mm;
          }
        }
      }
      mm = marker_allele_cts[ii * 4] + marker_allele_cts[ii * 4 + 1] + marker_allele_cts[ii * 4 + 2] + marker_allele_cts[ii * 4 + 3];
      if (marker_allele_cts[ii * 4] + marker_allele_cts[ii * 4 + 1] < geno_int_thresh) {
        exclude(marker_exclude, ii, &marker_exclude_ct);
      } else if (marker_allele_cts[ii * 4 + 1] < maf_int_thresh) {
        exclude(marker_exclude, ii, &marker_exclude_ct);
      } else if (marker_allele_cts[ii * 4 + 1] == marker_allele_cts[ii * 4 + 2]) {
        retval = RET_INVALID_FORMAT;
        printf("Error: Ambiguous minor allele at marker %s.\n", &(marker_id[(ii + 1) * max_marker_id_len]));
        if (!marker_allele_cts[ii * 4 + 1]) {
          printf("(Only one allele is present.  Consider using the --maf flag to avoid this\nproblem.)\n");
        }
        goto wdist_ret_1;
      }
    }
    if (indiv_exclude_ct > (unfiltered_indiv_ct - 2)) {
      retval = RET_INVALID_FORMAT;
      printf("Error: Too many people fail QC.\n");
      goto wdist_ret_1;
    } else if (marker_exclude_ct == unfiltered_marker_ct) {
      retval = RET_INVALID_FORMAT;
      printf("Error: All markers fail QC.\n");
      goto wdist_ret_1;
    }

    het_probs = (double*)malloc(((unfiltered_marker_ct - marker_exclude_ct) / 2 + 1) * sizeof(double));
    if (!het_probs) {
      goto wdist_ret_1;
    }
    marker_weights = (double*)wkspace_alloc((unfiltered_marker_ct - marker_exclude_ct) * sizeof(double));
    if (!marker_weights) {
      goto wdist_ret_2;
    }
    dptr2 = marker_weights;

    if (wkspace_left < unfiltered_marker_ct * 9 * sizeof(int)) {
      goto wdist_ret_1;
    }
    // unfortunately, this requires a third pass...
    hwe_ll = (int*)wkspace_base;
    hwe_lh = (int*)(&(wkspace_base[unfiltered_marker_ct * sizeof(int)]));
    hwe_hh = (int*)(&(wkspace_base[unfiltered_marker_ct * 2 * sizeof(int)]));
    hwe_u_ll = (int*)(&(wkspace_base[unfiltered_marker_ct * 3 * sizeof(int)]));
    hwe_u_lh = (int*)(&(wkspace_base[unfiltered_marker_ct * 4 * sizeof(int)]));
    hwe_u_hh = (int*)(&(wkspace_base[unfiltered_marker_ct * 5 * sizeof(int)]));
    hwe_a_ll = (int*)(&(wkspace_base[unfiltered_marker_ct * 6 * sizeof(int)]));
    hwe_a_lh = (int*)(&(wkspace_base[unfiltered_marker_ct * 7 * sizeof(int)]));
    hwe_a_hh = (int*)(&(wkspace_base[unfiltered_marker_ct * 8 * sizeof(int)]));
    fill_int_zero((int*)wkspace_base, unfiltered_marker_ct * 9);

    for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
      if (excluded(indiv_exclude, ii)) {
	continue;
      }
      fseeko(pedfile, line_locs[ii], SEEK_SET);
      fgets((char*)pedbuf, pedbuflen, pedfile);
      bufptr = (char*)pedbuf;
      for (jj = 0; jj < unfiltered_marker_ct; jj += 1) {
	if (excluded(marker_exclude, jj)) {
	  continue;
	}
	mm = 0;
	for (kk = 0; kk < 2; kk++) {
	  cc = *bufptr;
	  if (cc == marker_alleles[jj * 4 + 1]) {
	    mm += 1;
	  } else if (cc != marker_alleles[jj * 4]) {
	    mm = -2;
	  }
	  bufptr++;
	  while ((*bufptr == ' ') || (*bufptr == '\t')) {
	    bufptr++;
	  }
	}
	if (mm == 0) {
	  hwe_ll[jj] += 1;
	  if (bin_pheno) {
	    if (pheno_c[ii] == 0) {
	      hwe_u_ll[jj] += 1;
	    } else if (pheno_c[ii] == 1) {
	      hwe_a_ll[jj] += 1;
	    }
	  }
	} else if (mm == 1) {
	  hwe_lh[jj] += 1;
	  if (bin_pheno) {
	    if (pheno_c[ii] == 0) {
	      hwe_u_lh[jj] += 1;
	    } else if (pheno_c[ii] == 1) {
	      hwe_a_lh[jj] += 1;
	    }
	  }
	} else if (mm == 2) {
	  hwe_hh[jj] += 1;
	  if (bin_pheno) {
	    if (pheno_c[ii] == 0) {
	      hwe_u_hh[jj] += 1;
	    } else if (pheno_c[ii] == 1) {
	      hwe_a_hh[jj] += 1;
	    }
	  }
	}
      }
    }
    for (ii = 0; ii < unfiltered_marker_ct; ii++) {
      if (!excluded(marker_exclude, ii)) {
	if (SNPHWE_t(hwe_lh[ii], hwe_ll[ii], hwe_hh[ii], hwe_thresh)) {
	  exclude(marker_exclude, ii, &marker_exclude_ct);
	} else if (bin_pheno) {
	  if (SNPHWE_t(hwe_u_lh[ii], hwe_u_ll[ii], hwe_u_hh[ii], hwe_thresh) || SNPHWE_t(hwe_a_lh[ii], hwe_a_ll[ii], hwe_a_hh[ii], hwe_thresh)) {
	    exclude(marker_exclude, ii, &marker_exclude_ct);
	  }
	}
	if (!excluded(marker_exclude, ii)) {
	  *dptr2++ = calc_wt_mean(exponent, hwe_lh[ii], hwe_ll[ii], hwe_hh[ii]);
          // most common allele is counted in hwe_ll, not hwe_hh as in the
          // binary loads
          // should standardize this when implementing .ped load without forced
          // --make-bed
	  mafs[ii] = ((double)(hwe_ll[ii] * 2 + hwe_lh[ii])) / ((double)(hwe_ll[ii] + hwe_lh[ii] + hwe_hh[ii]) * 2);
	}
      }
    }
    free(het_probs);
    het_probs = NULL;
    rewind(pedfile);

    ulii = unfiltered_indiv_ct - indiv_exclude_ct;
    ulii = (ulii * (ulii - 1)) / 2;
    dists_alloc = ulii * (sizeof(int) + sizeof(double));
    
    llxx = malloc_size_mb * 1048576 - dists_alloc;
    geno_window_size = llxx / (unfiltered_indiv_ct - indiv_exclude_ct);
    unfiltered_marker_ct4 = (unfiltered_marker_ct - marker_exclude_ct + 3) / 4;
    uljj = (unfiltered_indiv_ct - indiv_exclude_ct + 3) / 4;
    if (make_bed && (uljj <= ((malloc_size_mb * 1048576) / (unfiltered_marker_ct - marker_exclude_ct)))) {
      snp_major = 1;
      printf("Writing binary files...");
      strcpy(outname_end, ".bed");
      bedtmpfile = fopen(outname, "wb");
      if (!bedtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_1;
      }
      strcpy(outname_end, ".bim");
      bimtmpfile = fopen(outname, "wb");
      if (!bimtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      strcpy(outname_end, ".fam");
      famtmpfile = fopen(outname, "wb");
      if (!famtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (3 != fwrite("l\x1b\x01", 1, 3, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_2;
      }
      if (wkspace_left < uljj * (unfiltered_marker_ct - marker_exclude_ct)) {
        goto wdist_ret_2;
      }
      ped_geno = wkspace;
      memset(ped_geno, 0, uljj * (unfiltered_marker_ct - marker_exclude_ct));
      // ----- .ped (fourth pass) -> .fam + .bed conversion, snp-major -----
      rewind(pedfile);
      ii = 0;
      oo = 0; // final person index
      last_tell = 0;

      while (fgets((char*)pedbuf, pedbuflen, pedfile) != NULL) {
        if (ii == unfiltered_indiv_ct) {
          break;
        }
        llxx = ftello(pedfile);
        if (llxx < line_locs[ii]) {
          last_tell = llxx;
          continue;
        }
        if (excluded(indiv_exclude, ii)) {
          last_tell = llxx;
          ii++;
          continue;
        }
        jj = (int)(line_locs[ii++] - last_tell);
        last_tell = llxx;
        pedbuf[jj - 1] = '\n';
        if (jj != fwrite(pedbuf, 1, jj, famtmpfile)) {
          retval = RET_WRITE_FAIL;
          printf("\n");
          goto wdist_ret_2;
        }
        bufptr = (char*)(&pedbuf[jj]);
        mm = 0; // final SNP index
        pp = uljj;
	for (jj = 0; jj < unfiltered_marker_ct; jj += 1) {
          if (excluded(marker_exclude, jj)) {
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
            continue;
          }
	  nn = 0;
	  for (kk = 0; kk < 2; kk++) {
	    cc = *bufptr;
	    if (cc == marker_alleles[jj * 4]) {
	      if (nn == 2) {
		nn = 3;
	      } else if (!nn) {
		nn = 2;
	      }
	    } else if (cc != marker_alleles[jj * 4 + 1]) {
              nn = 1;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	  }
	  ped_geno[mm * pp + (oo / 4)] |= (nn << ((oo % 4) * 2));
          mm++;
	}
	oo++;
      }
      ulii = uljj * (unfiltered_marker_ct - marker_exclude_ct);
      if (ulii != fwrite(ped_geno, 1, ulii, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_2;
      }
    } else if (unfiltered_marker_ct4 > geno_window_size) {
      if (make_bed) {
        printf("Writing individual-major binary files, since .ped is very large...");
      } else {
        printf(".ped file too large for direct read.  Converting to binary...");
      }
      if (make_bed) {
        strcpy(outname_end, ".bed");
      } else {
        strcpy(outname_end, ".bed.tmp");
      }
      bedtmpfile = fopen(outname, "wb");
      if (!bedtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_1;
      }
      if (make_bed) {
        strcpy(outname_end, ".bim");
      } else {
        strcpy(outname_end, ".bim.tmp");
      }
      bimtmpfile = fopen(outname, "wb");
      if (!bimtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (make_bed) {
        strcpy(outname_end, ".fam");
      } else {
        strcpy(outname_end, ".fam.tmp");
      }
      famtmpfile = fopen(outname, "wb");
      if (!famtmpfile) {
        retval = RET_OPENFAIL;
        printf("\nError: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (3 != fwrite("l\x1b", 1, 3, bedtmpfile)) {
	retval = RET_WRITE_FAIL;
        printf("\n");
	goto wdist_ret_2;
      }
      // ----- .ped (fourth pass) -> .fam + .bed conversion, indiv-major -----
      rewind(pedfile);
      ii = 0; // line count
      last_tell = 0; // last file location
      while (fgets((char*)pedbuf, pedbuflen, pedfile) != NULL) {
        if (ii == unfiltered_indiv_ct) {
          break;
        }
        llxx = ftello(pedfile);
        if (llxx < line_locs[ii]) {
          last_tell = llxx;
          continue;
        }
        if (excluded(indiv_exclude, ii)) {
          last_tell = llxx;
          ii++;
          continue;
        }
        jj = (int)(line_locs[ii++] - last_tell);
        last_tell = llxx;
        pedbuf[jj - 1] = '\n';
        if (jj != fwrite(pedbuf, 1, jj, famtmpfile)) {
          retval = RET_WRITE_FAIL;
          goto wdist_ret_2;
        }
        bufptr = (char*)(&pedbuf[jj]);
        memset(ped_geno, 0, unfiltered_marker_ct4);
        gptr = ped_geno;
        mm = 0;
	for (jj = 0; jj < unfiltered_marker_ct; jj += 1) {
          if (excluded(marker_exclude, jj)) {
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
            continue;
          }
	  nn = 0;
	  for (kk = 0; kk < 2; kk++) {
	    cc = *bufptr;
	    if (cc == marker_alleles[jj * 4]) {
	      if (nn == 2) {
		nn = 3;
	      } else if (!nn) {
		nn = 2;
	      }
	    } else if (cc != marker_alleles[jj * 4 + 1]) {
              nn = 1;
	    }
	    bufptr++;
	    while ((*bufptr == ' ') || (*bufptr == '\t')) {
	      bufptr++;
	    }
	  }
	  *gptr |= (nn << (mm * 2));
          mm = (mm + 1) % 4;
	  if (mm == 0) {
	    gptr++;
	  }
	}
	if (unfiltered_marker_ct4 != fwrite(ped_geno, 1, unfiltered_marker_ct4, bedtmpfile)) {
	  retval = RET_WRITE_FAIL;
	  goto wdist_ret_2;
	}
      }
    }
    if (make_bed || (unfiltered_marker_ct4 > geno_window_size)) {
      marker_alleles_tmp = (char*)malloc((unfiltered_marker_ct - marker_exclude_ct) * 2 * sizeof(char));
      if (!marker_alleles_tmp) {
        goto wdist_ret_2;
      }
      cptr = marker_alleles_tmp;
      iwptr = (int*)malloc((unfiltered_marker_ct - marker_exclude_ct) * 2 * sizeof(int));
      iptr = iwptr;
      for (ii = 0; ii < unfiltered_marker_ct; ii += 1) {
        if (excluded(marker_exclude, ii)) {
          continue;
        }
        *cptr++ = marker_alleles[ii * 4];
        *cptr++ = marker_alleles[ii * 4 + 1];
        *iptr++ = marker_allele_cts[ii * 4];
        *iptr++ = marker_allele_cts[ii * 4 + 1];
      }
      free(marker_alleles);
      marker_alleles = marker_alleles_tmp;
      marker_alleles_tmp = NULL;
      marker_allele_cts = iwptr;

      // ----- .map -> .bim conversion -----
      rewind(mapfile);
      kk = 0;
      for (ii = 0; ii < unfiltered_marker_ct; ii += 1) {
        do {
          fgets(tbuf, MAXLINELEN, mapfile);
        } while (tbuf[0] <= ' ');
        if (excluded(marker_exclude, ii)) {
          continue;
        }
        bufptr = next_item(tbuf);
        bufptr = next_item(bufptr);
        if (map_cols == 4) {
          bufptr = next_item(bufptr);
        }
        while (!is_space_or_eoln(*bufptr)) {
          bufptr++;
        }
        *bufptr++ = '\t';
        *bufptr++ = marker_alleles[kk * 2 + 1];
        *bufptr++ = '\t';
        *bufptr++ = marker_alleles[kk * 2];
        *bufptr++ = '\n';
        kk++;
        jj = (int)(bufptr - tbuf);
        if (jj != fwrite(tbuf, 1, jj, bimtmpfile)) {
	  retval = RET_WRITE_FAIL;
	  goto wdist_ret_2;
        }
      }
      printf(" done.\n");
      fclose(pedfile);
      fclose(bedtmpfile);
      fclose(bimtmpfile);
      fclose(famtmpfile);
      bedtmpfile = NULL;
      bimtmpfile = NULL;
      famtmpfile = NULL;
      if (make_bed) {
        strcpy(outname_end, ".bed");
      } else {
        strcpy(outname_end, ".bed.tmp");
      }
      pedfile = fopen(outname, "rb");
      if (!pedfile) {
        retval = RET_OPENFAIL;
        printf("Error: Failed to open %s.\n", outname);
        goto wdist_ret_2;
      }
      binary_files = 1;
      if (pheno_c) {
        collapse_phenoc(pheno_c, indiv_exclude, unfiltered_indiv_ct);
      } else if (pheno_d) {
        collapse_phenod(pheno_d, indiv_exclude, unfiltered_indiv_ct);
      }
      unfiltered_marker_ct -= marker_exclude_ct;
      marker_exclude_ct = 0;
      memset(marker_exclude, 0, (unfiltered_marker_ct + 7) / 8);
      unfiltered_indiv_ct -= indiv_exclude_ct;
      unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
      indiv_exclude_ct = 0;
      memset(indiv_exclude, 0, (unfiltered_indiv_ct + 7) / 8);
    }
  }

  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  if ((calculation_type & CALC_GROUPDIST) && (!bin_pheno)) {
    retval = RET_INVALID_CMDLINE;
    printf("Error: --groupdist calculation requires binary phenotype.\n");
    goto wdist_ret_2;
  } else if ((calculation_type & CALC_REGRESS_DISTANCE) && (!pheno_d)) {
    retval = RET_INVALID_CMDLINE;
    printf("Error: --regress-distance calculation requires scalar phenotype.\n");
    goto wdist_ret_2;
  } else if ((calculation_type & CALC_UNRELATED_HERITABILITY) && (!pheno_d)) {
    retval = RET_INVALID_CMDLINE;
    printf("Error: --unrelated-heritability requires scalar phenotype.\n");
    goto wdist_ret_2;
  }
  indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
  if (distance_req(calculation_type)) {
    // normalize marker weights to add to 2^32 - 1
    dxx = 0.0;
    dptr2 = marker_weights;
    dptr3 = &(marker_weights[marker_ct]);
    while (dptr2 < dptr3) {
      dxx += *dptr2++;
    }
    dxx = 4294967295.0 / dxx;
    dptr2 = marker_weights;
    giptr = (unsigned int*)marker_weights;
    while (dptr2 < dptr3) {
      *giptr++ = (unsigned int)((*dptr2++) * dxx + 0.5);
    }

    wkspace_reset((unsigned char*)marker_weights);
    marker_weights_i = (unsigned int*)wkspace_alloc(marker_ct * sizeof(int));
  }
  printf("%d markers and %d people pass filters and QC.\n", marker_ct, indiv_ct);
  if (thread_ct > 1) {
    printf("Using %d threads (change this with --threads).\n", thread_ct);
  }

  if (calculation_type & CALC_FREQ) {
    if (calculation_type & CALC_FREQ_GCTA) {
      strcpy(outname_end, ".freq");
      outfile = fopen(outname, "w");
      if (!outfile) {
	retval = RET_OPENFAIL;
	printf("Error: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      for (ii = 0; ii < unfiltered_marker_ct; ii++) {
        if (excluded(marker_exclude, ii)) {
          continue;
        }
        fprintf(outfile, "%s\t%c\t%g\n", &(marker_id[ii * max_marker_id_len]), marker_alleles[ii * 2], 1.0 - mafs[ii]);
      }
    } else {
      strcpy(outname_end, ".frq");
      outfile = fopen(outname, "w");
      if (!outfile) {
	retval = RET_OPENFAIL;
	printf("Error: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (max_marker_id_len < 4) {
	fprintf(outfile, " CHR  SNP   A1   A2          MAF  NCHROBS\n");
	strcpy(tbuf, "%4d  %3s    %c    %c  %11g  %7d\n");
      } else {
	sprintf(tbuf, " CHR  %%%lus   A1   A2          MAF  NCHROBS\n", max_marker_id_len - 1);
	fprintf(outfile, tbuf, "SNP");
	sprintf(tbuf, "%%4d  %%%lus    %%c    %%c  %%11g  %%7d\n", max_marker_id_len - 1);
      }
      for (ii = 0; ii < unfiltered_marker_ct; ii++) {
	if (excluded(marker_exclude, ii)) {
	  continue;
	}
	if (mafs[ii] < 0.5) {
	  fprintf(outfile, tbuf, marker_chrom[ii], &(marker_id[ii * max_marker_id_len]), marker_alleles[ii * 2 + 1], marker_alleles[ii * 2], mafs[ii], marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]);
	} else {
	  fprintf(outfile, tbuf, marker_chrom[ii], &(marker_id[ii * max_marker_id_len]), marker_alleles[ii * 2], marker_alleles[ii * 2 + 1], 1.0 - mafs[ii], marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]);
	}
      }
    }
    fclose(outfile);
    printf("Allele frequencies written to %s.\n", outname);
    outfile = NULL;
  }

  if (relationship_or_ibc_req(calculation_type)) {
    if (relationship_req(calculation_type)) {
      ulii = indiv_ct;
      ulii = (ulii * (ulii - 1)) / 2;
      dists_alloc = ulii * sizeof(double);
      rel_dists = (double*)wkspace_alloc(dists_alloc);
      if (!rel_dists) {
	goto wdist_ret_2;
      }
      fill_double_zero(rel_dists, ulii);
    }
    if (calculation_type & CALC_IBC) {
      ii = indiv_ct * 3;
    } else {
      ii = indiv_ct;
    }
    rel_ibc = (double*)wkspace_alloc(ii * sizeof(double));
    if (!rel_ibc) {
      goto wdist_ret_2;
    }
    fill_double_zero(rel_ibc, ii);
    wkspace_mark = wkspace_base;
    if (relationship_req(calculation_type)) {
      rel_missing = (int*)wkspace_alloc(dists_alloc);
      if (!rel_missing) {
	goto wdist_ret_2;
      }
      fill_int_zero(rel_missing, ulii);
    }
    indiv_missing = (unsigned int*)wkspace_alloc(indiv_ct * sizeof(int));
    if (!indiv_missing) {
      goto wdist_ret_2;
    }
    fill_int_zero((int*)indiv_missing, indiv_ct);
;
    if (binary_files && snp_major) {
      fseeko(pedfile, 3, SEEK_SET);
      ii = 0;
      pp = 0;
      ped_geno = wkspace_alloc(indiv_ct * sizeof(long));
      if (!ped_geno) {
        goto wdist_ret_2;
      }
      mmasks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
      if (!mmasks) {
	goto wdist_ret_2;
      }
      gptr = wkspace_alloc(MULTIPLEX_REL * unfiltered_indiv_ct4);
      if (!gptr) {
        goto wdist_ret_2;
      }
      masks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
      if (!masks) {
	goto wdist_ret_2;
      }

      triangle_fill(thread_start, indiv_ct, thread_ct, 1, CACHELINE_DBL);
      triangle_fill(thread_start0, indiv_ct, thread_ct, 0, CACHELINE_DBL);
      // See later comments on CALC_DISTANCE.
      // The difference is, we have to use + instead of XOR here to distinguish
      // the cases, so it's best to allow 3 bits per marker.
      // The rest is similar to CALC_DISTANCE on nonzero exponents: we
      // precalculate tables covering all possible pairwise differences between
      // two individuals on sets of 5 consecutive markers, slashing the number
      // of floating-point additions we need to perform by roughly 80%.
      while (pp < marker_ct) {
        jj = 0;
        while ((jj < MULTIPLEX_REL) && (pp < marker_ct)) {
          ulii = 0;
          while (excluded(marker_exclude, ii)) {
            ii++;
            ulii++;
          }
          if (ulii) {
            fseeko(pedfile, ulii * unfiltered_indiv_ct4, SEEK_CUR);
          }
          if (fread(&(gptr[jj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
            retval = RET_READ_FAIL;
            goto wdist_ret_2;
          }
          maf_buf[jj] = mafs[ii];
          ii++;
          jj++;
          pp++;
        }
        if (jj < MULTIPLEX_REL) {
          memset(&(gptr[jj * unfiltered_indiv_ct4]), 0, (MULTIPLEX_REL - jj) * unfiltered_indiv_ct4);
          fill_double_zero(&(maf_buf[jj]), MULTIPLEX_REL - jj);
        }
        fill_long_zero((long*)mmasks, indiv_ct);

	for (nn = 0; nn < jj; nn += MULTIPLEX_REL / 3) {
          fill_long_zero((long*)masks, indiv_ct);
          oo = 0;
          glptr2 = (unsigned long*)ped_geno;
	  for (qq = 0; oo < indiv_ct; qq++) {
	    while (excluded(indiv_exclude, qq)) {
              qq++;
            }
	    ulii = 0;
	    gptr2 = &(gptr[qq / 4 + nn * unfiltered_indiv_ct4]);
	    kk = (qq % 4) * 2;
	    for (mm = 0; mm < (MULTIPLEX_REL / 3); mm++) {
	      uljj = (gptr2[mm * unfiltered_indiv_ct4] >> kk) & 3;
	      if (uljj == 1) {
		masks[oo] |= 7LLU << (mm * 3);
		mmasks[oo] |= 1LLU << (nn + mm);
		indiv_missing[oo] += 1;
	      }
	      ulii |= uljj << (mm * 3);
	    }
	    *glptr2++ = ulii;
	    oo++;
	  }
          if (calculation_type & CALC_IBC) {
            for (oo = 0; oo < 3; oo++) {
              update_rel_ibc(&(rel_ibc[oo * indiv_ct]), (unsigned long*)ped_geno, &(maf_buf[nn]), oo);
            }
          } else {
            update_rel_ibc(rel_ibc, (unsigned long*)ped_geno, &(maf_buf[nn]), ibc_type);
          }
          if (relationship_req(calculation_type)) {
	    fill_weights_r(weights, &(maf_buf[nn]), var_std);
	    for (ulii = 1; ulii < thread_ct; ulii++) {
	      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_rel_thread, (void*)ulii)) {
		printf("Error: Could not create thread.\n");
		retval = RET_THREAD_CREATE_FAIL;
		goto wdist_ret_2;
	      }
	    }
	    incr_dists_r(rel_dists, (unsigned long*)ped_geno, masks, 0, weights);
	    for (oo = 0; oo < thread_ct - 1; oo++) {
	      pthread_join(threads[oo], NULL);
	    }
          }
	}
        if (relationship_req(calculation_type)) {
	  for (ulii = 1; ulii < thread_ct; ulii++) {
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_relm_thread, (void*)ulii)) {
	      printf("Error: Could not create thread.\n");
	      retval = RET_THREAD_CREATE_FAIL;
	      goto wdist_ret_2;
	    }
	  }
	  incr_dists_rm(rel_missing, 0);
	  for (oo = 0; oo < thread_ct - 1; oo++) {
	    pthread_join(threads[oo], NULL);
	  }
	}
        printf("\r%d markers complete.", pp);
        fflush(stdout);
      }
      if (calculation_type & (CALC_RELATIONSHIP_MASK | CALC_UNRELATED_HERITABILITY | CALC_REL_CUTOFF)) {
        printf("\rRelationship matrix calculation complete.\n");
        dist_ptr = rel_dists;
      } else {
        printf("\n");
      }
      ulii = indiv_ct;
      ulii = ulii * (ulii + 1) / 2;
      dptr2 = rel_ibc;
      if (calculation_type & CALC_IBC) {
        dptr3 = &(rel_ibc[indiv_ct]);
        dptr4 = &(rel_ibc[indiv_ct * 2]);
      }
      iwptr = rel_missing;
      for (ii = 0; ii < indiv_ct; ii++) {
        uii = marker_ct - indiv_missing[ii];
        if (calculation_type & (CALC_RELATIONSHIP_MASK | CALC_UNRELATED_HERITABILITY | CALC_REL_CUTOFF)) {
          giptr = indiv_missing;
	  for (jj = 0; jj < ii; jj++) {
	    *dist_ptr /= uii - (*giptr++) + (*iwptr++);
	    dist_ptr++;
	  }
        }
        if (calculation_type & CALC_IBC) {
          *dptr2 /= uii;
          dptr2++;
          *dptr3 /= uii;
          dptr3++;
          *dptr4 /= uii;
          dptr4++;
        } else {
          *dptr2 /= uii;
          dptr2++;
        }
      }
      if (calculation_type & CALC_REL_CUTOFF) {
        // Algorithm:
        // - Whenever there is at least one individual with exactly one
        // remaining too-close relation, prune the other side of that
        // relationship, because doing so is never suboptimal.
        // - Otherwise, there's no efficient rule that is always optimal
        // (assuming P != NP, anyway), so we use a simple heuristic: prune the
        // first individual with the largest number of remaining too-close
        // relationships.
	ii = 0; // total number of individuals excluded
        jj = 0; // number of individuals with exactly one too-close relation

        // number of too-close relations, -1 if excluded
	iptr = (int*)wkspace_alloc(indiv_ct * sizeof(int));
	fill_int_zero(iptr, indiv_ct);
        dist_ptr = rel_dists;
	for (kk = 1; kk < indiv_ct; kk++) {
	  for (mm = 0; mm < kk; mm++) {
	    if (*dist_ptr++ > grm_cutoff) {
	      iptr[kk] += 1;
	      iptr[mm] += 1;
	    }
	  }
	}
        for (kk = 0; kk < indiv_ct; kk++) {
          if (iptr[kk] == 1) {
            jj++;
          }
        }
	do {
          if (jj) {
            // there is at least one individual with exactly one too-close
            // relation left, find the first one
            kk = 0;
            while (iptr[kk] != 1) {
              kk++;
            }
            // and now find the identity of the other side
            dist_ptr = &(rel_dists[(kk * (kk - 1)) / 2]);
            for (mm = 0; mm < kk; mm++) {
              if (*dist_ptr > grm_cutoff) {
                *dist_ptr = 0.0;
                break;
              }
              dist_ptr++;
            }
            if (mm == kk) {
              do {
                mm++;
	        dist_ptr = &(rel_dists[(mm * (mm - 1)) / 2 + kk]);
              } while (*dist_ptr <= grm_cutoff);
              *dist_ptr = 0.0;
            }
            iptr[kk] = 0;
            jj--;
            if (iptr[mm] == 1) {
              jj--;
              iptr[mm] = -1;
              ii++;
              continue;
            }
          } else {
            // find identity of first individual with maximum number of
            // remaining too-close relations
            kk = 0; // highest too-close pair count
            mm = -1; // associated individual index
            for (nn = 0; nn < indiv_ct; nn++) {
              if (iptr[nn] > kk) {
                kk = iptr[nn];
                mm = nn;
              }
            }
            // no too-close relations left at all, we're done
            if (mm == -1) {
              break;
            }
          }
	  dist_ptr = &(rel_dists[(mm * (mm - 1)) / 2]);
	  for (kk = 0; kk < mm; kk++) {
	    if (*dist_ptr > grm_cutoff) {
	      *dist_ptr = 0.0;
	      iptr[kk] -= 1;
	      if (iptr[kk] < 2) {
		if (iptr[kk] == 0) {
		  jj--;
		} else if (iptr[kk] == 1) {
		  jj++;
		}
	      }
	    }
            dist_ptr++;
	  }
	  for (kk = mm + 1; kk < indiv_ct; kk++) {
	    dist_ptr = &(rel_dists[(kk * (kk - 1)) / 2 + mm]);
	    if (*dist_ptr > grm_cutoff) {
	      *dist_ptr = 0.0;
	      iptr[kk] -= 1;
	      if (iptr[kk] < 2) {
		if (iptr[kk] == 0) {
		  jj--;
		} else if (iptr[kk] == 1) {
		  jj++;
		}
	      }
	    }
	  }
	  iptr[mm] = -1;
	  ii++;
	} while (1);
	exclude_multi(indiv_exclude, iptr, indiv_ct, &indiv_exclude_ct);
        if (ii) {
	  dist_ptr = rel_dists; // write
	  dptr2 = rel_dists; // read
	  dptr3 = rel_ibc; // write
	  dptr4 = rel_ibc; // read
	  for (jj = 0; jj < indiv_ct; jj++) {
	    if (iptr[jj] != -1) {
	      if (calculation_type & CALC_IBC) {
		dptr3[indiv_ct] = dptr4[indiv_ct];
		dptr3[indiv_ct * 2] = dptr4[indiv_ct * 2];
	      }
	      *dptr3 = *dptr4++;
	      dptr3++;
	      for (kk = 0; kk < jj; kk++) {
		if (iptr[kk] != -1) {
		  *dist_ptr = *dptr2++;
		  dist_ptr++;
		} else {
		  dptr2++;
		}
	      }
	    } else {
	      dptr4++;
	      dptr2 = &(dptr2[jj]);
	    }
	  }
	  indiv_ct -= ii;
	  if (calculation_type & CALC_IBC) {
	    for (jj = 0; jj < indiv_ct; jj++) {
	      *dptr3++ = *dptr4++;
	    }
            dptr4 = &(dptr4[ii]);
            for (jj = 0; jj < indiv_ct; jj++) {
              *dptr3++ = *dptr4++;
            }
            giptr = indiv_missing;
            giptr2 = indiv_missing;
            for (jj = 0; jj < indiv_ct + ii; jj++) {
              if (iptr[jj] != -1) {
                *giptr = *giptr2++;
                giptr++;
              } else {
                giptr2++;
              }
            }
	  }
        }
        if (ii == 1) {
          printf("1 individual excluded by --rel-cutoff.\n");
        } else {
	  printf("%d individuals excluded by --rel-cutoff.\n", ii);
        }
      }

      if (calculation_type & CALC_IBC) {
	strcpy(outname_end, ".ibc");
	outfile = fopen(outname, "w");
	if (!outfile) {
          retval = RET_OPENFAIL;
	  printf("Error: Failed to open %s.\n", outname);
	  goto wdist_ret_2;
	}
        dptr2 = rel_ibc;
        dptr3 = &(rel_ibc[indiv_ct]);
        dptr4 = &(rel_ibc[indiv_ct * 2]);
        fprintf(outfile, "FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3\n");
	for (ii = 0; ii < indiv_ct; ii++) {
          fprintf(outfile, "%d\t%d\t%d\t%g\t%g\t%g\n", ii + 1, ii + 1, marker_ct - indiv_missing[ii], *dptr3++ - 1.0, *dptr4++ - 1.0, *dptr2++ - 1.0);
	}
	fclose(outfile);
        printf("%s written.\n", outname);
	outfile = NULL;
      }
      if (calculation_type & CALC_RELATIONSHIP_MASK) {
        if (calculation_type & CALC_IBC) {
          dptr2 = &(rel_ibc[ibc_type * indiv_ct]);
        } else {
          dptr2 = rel_ibc;
        }
        mm = 1;
	if ((calculation_type & CALC_RELATIONSHIP_BIN) && !(calculation_type & CALC_RELATIONSHIP_GZ)) {
          if (calculation_type & CALC_RELATIONSHIP_SQ0) {
            fill_double_zero((double*)ped_geno, indiv_ct - 1);
          }
          if (calculation_type & CALC_RELATIONSHIP_GRM) {
            strcpy(outname_end, ".grm.bin");
          } else {
            strcpy(outname_end, ".rel.bin");
          }
          outfile = fopen(outname, "wb");
          for (ii = 0; ii < indiv_ct; ii++) {
            fwrite(&(rel_dists[(ii * (ii - 1)) / 2]), 1, ii * sizeof(double), outfile);
            fwrite(dptr2++, 1, sizeof(double), outfile);
            if (calculation_type & CALC_RELATIONSHIP_SQ0) {
              fwrite(ped_geno, 1, (indiv_ct - ii - 1) * sizeof(double), outfile);
            } else {
              for (jj = ii + 1; jj < indiv_ct; jj++) {
                fwrite(&(rel_dists[(jj * (jj - 1) / 2) + ii]), 1, sizeof(double), outfile);
              }
            }
	    if ((ii + 1) * 100 >= mm * indiv_ct) {
	      mm = (ii + 1) / indiv_ct;
	      printf("\rWriting... %d%%", mm++);
	      fflush(stdout);
	    }
          }
          fclose(outfile);
          outfile = NULL;
        } else if (calculation_type & CALC_RELATIONSHIP_GRM) {
          iwptr = rel_missing;
	  if (calculation_type & CALC_RELATIONSHIP_GZ) {
	    strcpy(outname_end, ".grm.gz");
	    gz_outfile = gzopen(outname, "wb");
	    if (!gz_outfile) {
              retval = RET_OPENFAIL;
	      printf("Error: Failed to open %s.\n", outname);
	      goto wdist_ret_2;
	    }
	    dist_ptr = rel_dists;
	    for (ii = 0; ii < indiv_ct; ii += 1) {
	      if ((long long)ii * (ii + 1) * 100 >= (long long)indiv_ct * (indiv_ct + 1) * mm) {
		mm = ((long long)ii * (ii + 1) * 100) / ((long long)indiv_ct * (indiv_ct + 1));
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	      for (jj = 0; jj < ii; jj += 1) {
		kk = marker_ct - *iwptr++;
		gzprintf(gz_outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dist_ptr++);
	      }
              kk = marker_ct - indiv_missing[ii];
              gzprintf(gz_outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dptr2++);
	    }
	    gzclose(gz_outfile);
	    gz_outfile = NULL;
	  } else {
	    strcpy(outname_end, ".grm");
	    outfile = fopen(outname, "w");
	    if (!outfile) {
              retval = RET_OPENFAIL;
	      printf("Error: Failed to open %s.\n", outname);
	      goto wdist_ret_2;
	    }
	    dist_ptr = rel_dists;
	    for (ii = 0; ii < indiv_ct; ii += 1) {
	      if ((long long)ii * (ii + 1) * 100 >= (long long)indiv_ct * (indiv_ct + 1) * mm) {
		mm = ((long long)ii * (ii + 1) * 100) / ((long long)indiv_ct * (indiv_ct + 1));
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	      for (jj = 0; jj < ii; jj += 1) {
		kk = marker_ct - *iwptr++;
		fprintf(outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dist_ptr++);
	      }
              kk = marker_ct - *iwptr++;
              fprintf(outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dptr2++);
	    }
	    fclose(outfile);
	    outfile = NULL;
	  }
        } else {
	  if (calculation_type & CALC_RELATIONSHIP_SQ0) {
	    cptr2 = (char*)(&ulii);
	    for (ii = 0; ii < sizeof(long); ii += 2) {
	      cptr2[ii] = '\t';
	      cptr2[ii + 1] = '0';
	    }
	    ii = (indiv_ct * 2 + sizeof(long) - 4) / sizeof(long);
	    glptr2 = (unsigned long*)ped_geno;
	    for (jj = 0; jj < ii; jj++) {
	      *glptr2++ = ulii;
	    }
	  }
	  if ((calculation_type & CALC_RELATIONSHIP_GZ) && !(calculation_type & CALC_RELATIONSHIP_BIN)) {
	    strcpy(outname_end, ".rel.gz");
	    gz_outfile = gzopen(outname, "wb");
	    if (!gz_outfile) {
              retval = RET_OPENFAIL;
	      printf("Error: Failed to open %s.\n", outname);
	      goto wdist_ret_2;
	    }
	    dist_ptr = rel_dists;
            gzprintf(gz_outfile, "%g", *dptr2++);
            if (calculation_type & CALC_RELATIONSHIP_SQ0) {
              gzwrite(gz_outfile, ped_geno, (indiv_ct - 1) * 2);
            } else if (calculation_type & CALC_RELATIONSHIP_SQ) {
              for (jj = 1; jj < indiv_ct; jj++) {
                gzprintf(gz_outfile, "\t%g", rel_dists[(jj * (jj - 1)) / 2]);
	      }
            }
            gzprintf(gz_outfile, "\n");
	    for (ii = 1; ii < indiv_ct; ii++) {
	      gzprintf(gz_outfile, "%g", *dist_ptr++);
	      for (jj = 1; jj < ii; jj++) {
		gzprintf(gz_outfile, "\t%g", *dist_ptr++);
	      }
              gzprintf(gz_outfile, "\t%g", *dptr2++);
	      if (calculation_type & CALC_RELATIONSHIP_SQ0) {
		gzwrite(gz_outfile, ped_geno, (indiv_ct - jj - 1) * 2);
		if ((ii + 1) * 100 >= mm * indiv_ct) {
		  mm = (ii + 1) / indiv_ct;
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
	      } else {
                if (calculation_type & CALC_RELATIONSHIP_SQ) {
                  for (jj = ii + 1; jj < indiv_ct; jj++) {
                    gzprintf(gz_outfile, "\t%g", rel_dists[((jj * (jj - 1)) / 2) + ii]);
                  }
                }
		if ((long long)(ii + 1) * (ii + 2) * 100 >= (long long)indiv_ct * (indiv_ct + 1) * mm) {
		  mm = ((long long)(ii + 1) * (ii + 2) * 100) / ((long long)indiv_ct * (indiv_ct + 1));
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
              }
	      gzprintf(gz_outfile, "\n");
	    }
	    gzclose(gz_outfile);
	    gz_outfile = NULL;
	  } else {
	    strcpy(outname_end, ".rel");
	    outfile = fopen(outname, "w");
	    if (!outfile) {
              retval = RET_OPENFAIL;
	      printf("Error: Failed to open %s.\n", outname);
	      goto wdist_ret_2;
	    }
	    dist_ptr = rel_dists;
            fprintf(outfile, "%g", *dptr2++);
            if (calculation_type & CALC_RELATIONSHIP_SQ0) {
              fwrite(ped_geno, 1, (indiv_ct - 1) * 2, outfile);
            } else if (calculation_type & CALC_RELATIONSHIP_SQ) {
              for (jj = 1; jj < indiv_ct; jj++) {
                fprintf(outfile, "\t%g", rel_dists[(jj * (jj - 1)) / 2]);
              }
            }
            fprintf(outfile, "\n");
	    for (ii = 1; ii < indiv_ct; ii++) {
	      fprintf(outfile, "%g", *dist_ptr++);
	      for (jj = 1; jj < ii; jj++) {
		fprintf(outfile, "\t%g", *dist_ptr++);
	      }
              fprintf(outfile, "\t%g", *dptr2++);
	      if (calculation_type & CALC_RELATIONSHIP_SQ0) {
		fwrite(ped_geno, 1, (indiv_ct - jj - 1) * 2, outfile);
		if ((ii + 1) * 100 >= mm * indiv_ct) {
		  mm = (ii + 1) / indiv_ct;
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
	      } else {
                if (calculation_type & CALC_RELATIONSHIP_SQ) {
                  for (jj = ii + 1; jj < indiv_ct; jj++) {
                    fprintf(outfile, "\t%g", rel_dists[((jj * (jj - 1)) / 2) + ii]);
                  }
                }
		if ((long long)(ii + 1) * (ii + 2) * 100 >= (long long)indiv_ct * (indiv_ct + 1) * mm) {
		  mm = ((long long)(ii + 1) * (ii + 2) * 100) / ((long long)indiv_ct * (indiv_ct + 1));
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
              }
	      fprintf(outfile, "\n");
	    }
	    fclose(outfile);
	    outfile = NULL;
	  }
	}
	printf("\rRelationship matrix written to %s.\n", outname);
	strcpy(&(outname_end[4]), ".id");
	outfile = fopen(outname, "w");
	if (!outfile) {
          retval = RET_OPENFAIL;
	  printf("Error: Failed to open %s.\n", outname);
	  goto wdist_ret_2;
	}
	for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
	  if (!excluded(indiv_exclude, ii)) {
	    fprintf(outfile, "%s\n", &(person_id[ii * max_person_id_len]));
	  }
	}
      }
      wkspace_reset(wkspace_mark);
    } else {
      printf("Error: Relationship calculation currently doesn't support individual-major\n.bed file.\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto wdist_ret_2;
    }

#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      ulii = indiv_ct;
      ulii = CACHEALIGN_DBL(ulii * ulii);
      dptr4 = &(rel_dists[ulii]);
      ulii = ulii * 3 + CACHEALIGN_DBL(indiv_ct) * 3;
      wkspace_reset((unsigned char*)rel_dists);
      rel_dists = (double*)wkspace_alloc(ulii * sizeof(double));
      if (!rel_dists) {
        goto wdist_ret_2;
      }
      dptr2 = &(rel_dists[ulii - CACHEALIGN_DBL(indiv_ct)]);
      collapse_copy_phenod(dptr2, pheno_d, indiv_exclude, unfiltered_indiv_ct);
      dxx = 0.0;
      dyy = 0.0;
      dptr3 = dptr2;
      dist_ptr = &(dptr2[indiv_ct]);
      while (dptr3 < dist_ptr) {
        dzz = *dptr3++;
        dxx += dzz;
        dyy += dzz * dzz;
      }
      dxx /= (double)indiv_ct;
      dxx = 1 / sqrt((dyy / (double)indiv_ct) - dxx * dxx);
      dptr3 = dptr2;
      while (dptr3 < dist_ptr) {
        *dptr3 *= dxx;
        dptr3++;
      }
      if (calculation_type & CALC_IBC) {
        dptr3 = &(rel_ibc[ibc_type * indiv_ct]);
      } else {
        dptr3 = rel_ibc;
      }
      for (ulii = 0; ulii < indiv_ct; ulii++) {
	memcpy(&(dptr4[ulii * indiv_ct]), &(rel_dists[(ulii * (ulii - 1)) / 2]), ulii * sizeof(double));
        dptr4[ulii * (indiv_ct + 1)] = *dptr3++;
        for (uljj = ulii + 1; uljj < indiv_ct; uljj++) {
          dptr4[ulii * indiv_ct + uljj] = rel_dists[(uljj * (uljj - 1)) / 2 + ulii];
        }
      }
      reml_em_one_trait(rel_dists, dptr2, &unrelated_herit_covg, &unrelated_herit_cove, unrelated_herit_tol, calculation_type & CALC_UNRELATED_HERITABILITY_STRICT);
      printf("h^2 estimate: %g\n", unrelated_herit_covg);
    }
#endif
    wkspace_reset((unsigned char*)rel_dists);
  }

  if (calculation_type & CALC_WRITE_SNPLIST) {
    strcpy(outname_end, ".snplist");
    outfile = fopen(outname, "w");
    if (!outfile) {
      retval = RET_OPENFAIL;
      printf("Error: Failed to open %s.\n", outname);
      goto wdist_ret_2;
    }
    for (ii = 0; ii < unfiltered_marker_ct; ii++) {
      if (!excluded(marker_exclude, ii)) {
        fprintf(outfile, "%s\n", &(marker_id[ii * max_marker_id_len]));
      }
    }
    fclose(outfile);
    outfile = NULL;
  }

  if (distance_req(calculation_type)) {
    ulii = indiv_ct;
    ulii = (ulii * (ulii - 1)) / 2;
    dists_alloc = ulii * sizeof(double);
    // additional + CACHELINE is to fix weird-ass aliasing bug that shows up
    // with -O2 in some cases
    dists = (double*)wkspace_alloc(dists_alloc + CACHELINE);
    if (!dists) {
      goto wdist_ret_2;
    }
    wkspace_mark = wkspace_base;
    missing_tot_weights = (unsigned int*)wkspace_alloc(ulii * sizeof(int));
    if (!missing_tot_weights) {
      goto wdist_ret_2;
    }
    fill_int_zero((int*)missing_tot_weights, ulii);

    if (exp0) {
      idists = (int*)(((char*)missing_tot_weights) - CACHEALIGN(ulii * sizeof(int)));
      fill_int_zero(idists, ulii);
      masks = (unsigned long*)wkspace_alloc(indiv_ct * (MULTIPLEX_2DIST / 8));
    } else {
      fill_double_zero(dists, ulii);
      masks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
    }
    if (!masks) {
      goto wdist_ret_2;
    }
    indiv_missing = (unsigned int*)wkspace_alloc(indiv_ct * sizeof(int));
    if (!indiv_missing) {
      goto wdist_ret_2;
    }
    fill_int_zero((int*)indiv_missing, indiv_ct);
    mmasks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
    if (!mmasks) {
      goto wdist_ret_2;
    }

    if (exp0) {
      multiplex = MULTIPLEX_DIST;
      ped_geno = wkspace_alloc(indiv_ct * (MULTIPLEX_2DIST / 8));
    } else {
      multiplex = MULTIPLEX_DIST_EXP;
      ped_geno = wkspace_alloc(indiv_ct * sizeof(long));
    }
    if (!ped_geno) {
      goto wdist_ret_2;
    }

    if (binary_files) {
      if (snp_major) {
	if (pedbuflen < multiplex * unfiltered_indiv_ct4) {
	  free(pedbuf);
	  pedbuf = (unsigned char*)malloc(multiplex * unfiltered_indiv_ct4 * sizeof(char));
	  if (!pedbuf) {
	    goto wdist_ret_2;
	  }
	}
	fseeko(pedfile, 3, SEEK_SET);
        if (exp0) {
          ii = sizeof(int);
        } else {
          ii = sizeof(double);
        }
	triangle_fill(thread_start, indiv_ct, thread_ct, 1, CACHELINE / ii);
	ii = 0; // current SNP index
	pp = 0; // after subtracting out excluded
	while (pp < marker_ct) {
	  for (jj = 0; jj < multiplex; jj++) {
	    maf_buf[jj] = 0.5;
	  }
          fill_int_zero((int*)wtbuf, multiplex);
	  jj = 0; // actual SNPs read
	  // If exponent == 0.0, we just need to calculate bitwise Hamming
          // distance between two strings (throwing in some masking if there
          // are missing markers).  Since the x86 processor has (at least until
          // very recently) lacked a popcount instruction, others have spent a
          // fair bit of effort on designing efficient popcount
          // implementations.  The popcount_..._multibyte functions reimplement
          // the most efficient known 64-bit algorithm (developed by Cedric
          // Lauradoux and Kim Walisch) using SSE2 instructions, and, to the
          // best of my knowledge, is faster than all other x86_64 open source
          // implementations.  (Other existing SSE2/SSSE3 implementations are
          // based on simpler algorithms.  See
          // http://www.dalkescientific.com/writings/diary/archive/2011/11/02/
          // faster_popcount_update.html for more information.)
          //
          // For nonzero exponents, we create lookup tables for all possible
          // combinations of 6-7 markers (splitting up the 32 markers that fit
          // into a 64-bit word 7-7-6-6-6; unfortunately, an 8-8-8-8 split is
          // too hard on the L2 cache of many machines).  Floating point adds
          // are the primary bottleneck here.
	  while ((jj < multiplex) && (pp < marker_ct)) {
	    while (excluded(marker_exclude, ii)) {
	      ii++;
	      fseeko(pedfile, (off_t)unfiltered_indiv_ct4, SEEK_CUR);
	    }
	    if (fread(&(pedbuf[jj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	      retval = RET_READ_FAIL;
	      goto wdist_ret_2;
	    }
	    maf_buf[jj] = mafs[ii];
            wtbuf[jj++] = marker_weights_i[pp++];
	    ii++;
	  }
	  if (jj < multiplex) {
	    memset(&(pedbuf[jj * unfiltered_indiv_ct4]), 0, (multiplex - jj) * unfiltered_indiv_ct4);
            if (exp0) {
              fill_long_zero((long*)ped_geno, indiv_ct * (MULTIPLEX_2DIST / BITCT));
              fill_long_zero((long*)masks, indiv_ct * (MULTIPLEX_2DIST / BITCT));
	    } else {
              fill_long_zero((long*)ped_geno, indiv_ct);
              fill_long_zero((long*)masks, indiv_ct);
            }
	  }
	  if (exp0) {
            for (nn = 0; nn < jj; nn += BITCT) {
	      glptr = &(((unsigned long*)ped_geno)[nn / BITCT2]);
	      glptr2 = &(masks[nn / BITCT2]);
              glptr3 = mmasks;
              giptr = indiv_missing;
	      for (oo = 0; oo < unfiltered_indiv_ct; oo++) {
		if (!excluded(indiv_exclude, oo)) {
		  kk = (oo % 4) * 2;
		  ulii = 0;
                  ulkk = 0;
                  gptr = &(pedbuf[oo / 4 + nn * unfiltered_indiv_ct4]);
		  for (mm = 0; mm < BITCT2; mm++) {
		    uljj = (gptr[mm * unfiltered_indiv_ct4] >> kk) & 3;
		    ulii |= uljj << (mm * 2);
                    if (uljj == 1) {
                      ulkk |= 1LLU << mm;
                      *giptr += wtbuf[mm + nn];
                    }
		  }
		  // use xor to convert representation to 0 = missing,
		  // 1 or 2 = homozygote, 3 = heterozygote
#if __LP64__
		  ulii ^= 0x5555555555555555LLU;
		  *glptr++ = ulii;
		  ulii = (ulii | (ulii >> 1)) & 0x5555555555555555LLU;
#else
		  ulii ^= 0x55555555;
		  *glptr++ = ulii;
		  ulii = (ulii | (ulii >> 1)) & 0x55555555;
#endif
		  *glptr2++ = ulii * 3;
                  *glptr3 = ulkk;
		  ulii = 0;
                  ulkk = 0;
                  gptr = &(pedbuf[oo / 4 + (nn + BITCT2) * unfiltered_indiv_ct4]);
		  for (mm = 0; mm < BITCT2; mm++) {
		    uljj = (gptr[mm * unfiltered_indiv_ct4] >> kk) & 3;
		    ulii |= uljj << (mm * 2);
                    if (uljj == 1) {
                      ulkk |= 1LLU << mm;
                      *giptr += wtbuf[mm + nn + BITCT2];
                    }
		  }
#if __LP64__
		  ulii ^= 0x5555555555555555LLU;
		  *glptr = ulii;
		  ulii = (ulii | (ulii >> 1)) & 0x5555555555555555LLU;
#else
		  ulii ^= 0x55555555;
		  *glptr = ulii;
		  ulii = (ulii | (ulii >> 1)) & 0x55555555;
#endif
		  *glptr2 = ulii * 3;
                  *glptr3++ |= ulkk << BITCT2;
                  glptr = &(glptr[(MULTIPLEX_2DIST / BITCT) - 1]);
                  glptr2 = &(glptr2[(MULTIPLEX_2DIST / BITCT) - 1]);
                  giptr++;
		}
	      }
	      weights_i = &(wtbuf[nn]);

	      for (ulii = 1; ulii < thread_ct; ulii++) {
		if (pthread_create(&(threads[ulii - 1]), NULL, &calc_distm_thread, (void*)ulii)) {
		  printf("Error: Could not create thread.\n");
		  retval = RET_THREAD_CREATE_FAIL;
		  goto wdist_ret_2;
		}
	      }
	      decr_dist_missing(missing_tot_weights, 0);
	      for (oo = 0; oo < thread_ct - 1; oo++) {
		pthread_join(threads[oo], NULL);
	      }
            }
	    for (ulii = 1; ulii < thread_ct; ulii++) {
	      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_idist_thread, (void*)ulii)) {
		printf("Error: Could not create thread.\n");
		retval = RET_THREAD_CREATE_FAIL;
		goto wdist_ret_2;
	      }
	    }
	    incr_dists_i(idists, (unsigned long*)ped_geno, 0);
	    for (oo = 0; oo < thread_ct - 1; oo++) {
	      pthread_join(threads[oo], NULL);
	    }
	  } else {
            fill_long_zero((long*)mmasks, indiv_ct);
            for (nn = 0; nn < jj; nn += MULTIPLEX_DIST_EXP / 2) {
	      glptr = (unsigned long*)ped_geno;
	      glptr2 = masks;
	      glptr3 = mmasks;
              giptr3 = indiv_missing;
	      for (oo = 0; oo < unfiltered_indiv_ct; oo++) {
		if (!excluded(indiv_exclude, oo)) {
		  kk = (oo % 4) * 2;
		  ulii = 0;
		  ulkk = 0;
		  gptr = &(pedbuf[oo / 4 + nn * unfiltered_indiv_ct4]);
		  for (mm = 0; mm < MULTIPLEX_DIST_EXP / 2; mm++) {
		    uljj = (gptr[mm * unfiltered_indiv_ct4] >> kk) & 3;
		    ulii |= uljj << (mm * 2);
		    if (uljj == 1) {
		      ulkk |= 1LLU << mm;
                      *giptr3 += wtbuf[mm + nn];
		    }
		  }
#if __LP64__
		  ulii ^= 0x5555555555555555LLU;
		  *glptr++ = ulii;
		  ulii = (ulii | (ulii >> 1)) & 0x5555555555555555LLU;
#else
                  ulii ^= 0x05555555;
		  *glptr++ = ulii;
		  ulii = (ulii | (ulii >> 1)) & 0x05555555;
#endif
		  *glptr2++ = ulii * 3;
		  *glptr3++ |= ulkk << nn;
                  giptr3++;
		}
	      }
	      fill_weights(weights, &(maf_buf[nn]), exponent);
	      for (ulii = 1; ulii < thread_ct; ulii++) {
		if (pthread_create(&(threads[ulii - 1]), NULL, &calc_dist_thread, (void*)ulii)) {
		  printf("Error: Could not create thread.\n");
		  retval = RET_THREAD_CREATE_FAIL;
		  goto wdist_ret_2;
		}
	      }
	      incr_dists(dists, (unsigned long*)ped_geno, 0);
	      for (oo = 0; oo < thread_ct - 1; oo++) {
		pthread_join(threads[oo], NULL);
	      }
	    }
	    weights_i = wtbuf;
	    for (ulii = 1; ulii < thread_ct; ulii++) {
	      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_distm_thread, (void*)ulii)) {
		printf("Error: Could not create thread.\n");
		retval = RET_THREAD_CREATE_FAIL;
		goto wdist_ret_2;
	      }
	    }
	    decr_dist_missing(missing_tot_weights, 0);
	    for (oo = 0; oo < thread_ct - 1; oo++) {
	      pthread_join(threads[oo], NULL);
	    }
          }
	  printf("\r%d markers complete.", pp);
	  fflush(stdout);
	}
      } else {
	printf("indiv-major distance calculation not done.\n");
        retval = RET_CALC_NOT_YET_SUPPORTED;
	goto wdist_ret_2;
      }
    } else {
      printf("text distance calculation not done (use --make-bed).\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto wdist_ret_2;
    }
    ulii = indiv_ct;
    ulii = ulii * (ulii - 1) / 2;
    giptr = missing_tot_weights;
    dptr2 = dists;
    if (exp0) {
      iptr = idists;
      for (ii = 1; ii < indiv_ct; ii++) {
        giptr2 = indiv_missing;
        uii = giptr2[ii];
        for (jj = 0; jj < ii; jj++) {
          *dptr2 = (4294967295.0 / ((4294967295U - uii - (*giptr2++)) + (*giptr++))) * (*iptr++);
          dptr2++;
        }
      }
    } else {
      for (ii = 1; ii < indiv_ct; ii++) {
        giptr2 = indiv_missing;
        uii = giptr2[ii];
        for (jj = 0; jj < ii; jj++) {
	  *dptr2 *= (4294967295.0 / ((4294967295U - uii - (*giptr2++)) + (*giptr++)));
	  dptr2++;
        }
      }
    }
    printf("\rDistance matrix calculation complete.\n");
  }

  if (calculation_type & CALC_DISTANCE_MASK) {
    if (calculation_type & CALC_DISTANCE_SQ0) {
      cptr2 = (char*)(&ulii);
      for (ii = 0; ii < sizeof(long); ii += 2) {
	cptr2[ii] = '\t';
	cptr2[ii + 1] = '0';
      }
      ii = (indiv_ct * 2 + sizeof(long) - 2) / sizeof(long);
      glptr2 = (unsigned long*)ped_geno;
      for (jj = 0; jj < ii; jj++) {
	*glptr2++ = ulii;
      }
    }
    kk = 1;
    if ((calculation_type & CALC_DISTANCE_GZBIN_MASK) == CALC_DISTANCE_GZ) {
      strcpy(outname_end, ".dist.gz");
      gz_outfile = gzopen(outname, "wb");
      if (!gz_outfile) {
	printf("Error: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (calculation_type & CALC_DISTANCE_SQ0) {
	gzwrite(gz_outfile, &(ped_geno[1]), indiv_ct * 2 - 1);
	gzprintf(gz_outfile, "\n");
      } else if (calculation_type & CALC_DISTANCE_SQ) {
        gzprintf(gz_outfile, "0");
        for (jj = 1; jj < indiv_ct; jj++) {
          gzprintf(gz_outfile, "\t%g", dists[(jj * (jj - 1)) / 2]);
        }
        gzprintf(gz_outfile, "\n");
      }
      dist_ptr = dists;
      for (ii = 1; ii < indiv_ct; ii++) {
	gzprintf(gz_outfile, "%g", *dist_ptr++);
	for (jj = 1; jj < ii; jj++) {
	  gzprintf(gz_outfile, "\t%g", *dist_ptr++);
	}
	if (calculation_type & CALC_DISTANCE_SQ0) {
	  gzwrite(gz_outfile, ped_geno, (indiv_ct - jj) * 2);
          if (ii * 100 >= (kk * indiv_ct)) {
            kk = (ii * 100) / indiv_ct;
            printf("\rWriting... %d%%", kk++);
            fflush(stdout);
          }
	} else {
          if (calculation_type & CALC_DISTANCE_SQ) {
            gzprintf(gz_outfile, "\t0");
            for (jj = ii + 1; jj < indiv_ct; jj++) {
              gzprintf(gz_outfile, "\t%g", dists[((jj * (jj - 1)) / 2) + ii]);
	    }
          }
          if ((long long)ii * (ii + 1) * 100 >= (long long)indiv_ct * (indiv_ct - 1) * kk) {
            kk = ((long long)ii * (ii + 1) * 100) / ((long long)indiv_ct * (indiv_ct - 1));
            printf("\rWriting... %d%%", kk++);
            fflush(stdout);
          }
        }
	gzprintf(gz_outfile, "\n");
      }
      gzclose(gz_outfile);
      gz_outfile = NULL;
    } else if ((calculation_type & CALC_DISTANCE_GZBIN_MASK) == CALC_DISTANCE_BIN) {
      strcpy(outname_end, ".dist.bin");
      outfile = fopen(outname, "wb");
      if (!outfile) {
	printf("Error: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      // ulii = indiv_ct;
      // ulii = ulii * ulii * sizeof(double);
      // uljj = fseeko(outfile, ulii - 1, SEEK_SET);
      // if (uljj == -1) {
      //   printf("Error: Failed to write final byte of %s.\n", outname);
      //   goto wdist_ret_2;
      // }
      // if (fwrite("\0", 1, 1, outfile) != 1) {
      //   printf("Error: Failed to write final byte of %s.\n", outname);
      // }
      // fseeko(outfile, 0, SEEK_SET);
      if (calculation_type & CALC_DISTANCE_SQ0) {
        fill_double_zero((double*)ped_geno, indiv_ct);
      } else {
        dxx = 0.0;
      }
      dist_ptr = dists;
      for (ii = 0; ii < indiv_ct; ii++) {
	fwrite(dist_ptr, 1, ii * sizeof(double), outfile);
	dist_ptr = &(dist_ptr[ii]);
	if (calculation_type & CALC_DISTANCE_SQ0) {
	  fwrite(ped_geno, 1, (indiv_ct - ii) * sizeof(double), outfile);
	} else {
	  fwrite(&dxx, 1, sizeof(double), outfile);
	  for (jj = ii + 1; jj < indiv_ct; jj++) {
	    fwrite(&(dists[(jj * (jj - 1)) / 2 + ii]), 1, sizeof(double), outfile);
	  }
	}
	if (ii * 100 >= (kk * indiv_ct)) {
	  kk = (ii * 100) / indiv_ct;
	  printf("\rWriting... %d%%", kk++);
	  fflush(stdout);
	}
      }
      fclose(outfile);
      outfile = NULL;
    } else {
      strcpy(outname_end, ".dist");
      outfile = fopen(outname, "w");
      if (!outfile) {
	printf("Error: Failed to open %s.\n", outname);
	goto wdist_ret_2;
      }
      if (calculation_type & CALC_DISTANCE_SQ0) {
	fwrite(&(ped_geno[1]), 1, indiv_ct * 2 - 1, outfile);
	fprintf(outfile, "\n");
      } else if (calculation_type & CALC_DISTANCE_SQ) {
        fprintf(outfile, "0");
        for (jj = 1; jj < indiv_ct; jj++) {
          fprintf(outfile, "\t%g", dists[(jj * (jj - 1)) / 2]);
        }
        fprintf(outfile, "\n");
      }
      dist_ptr = dists;
      for (ii = 1; ii < indiv_ct; ii++) {
	fprintf(outfile, "%g", *dist_ptr++);
	for (jj = 1; jj < ii; jj++) {
	  fprintf(outfile, "\t%g", *dist_ptr++);
	}
	if (calculation_type & CALC_DISTANCE_SQ0) {
	  fwrite(ped_geno, 1, (indiv_ct - jj) * 2, outfile);
	  if (ii * 100 >= (kk * indiv_ct)) {
	    kk = (ii * 100) / indiv_ct;
	    printf("\rWriting... %d%%", kk++);
	    fflush(stdout);
	  }
	} else {
          if (calculation_type & CALC_DISTANCE_SQ) {
            fprintf(outfile, "\t0");
            for (jj = ii + 1; jj < indiv_ct; jj++) {
              fprintf(outfile, "\t%g", dists[((jj * (jj - 1)) / 2) + ii]);
            }
          }
          if ((long long)ii * (ii + 1) * 100 >= (long long)indiv_ct * (indiv_ct - 1) * kk) {
            kk = ((long long)ii * (ii + 1) * 100) / ((long long)indiv_ct * (indiv_ct - 1));
            printf("\rWriting... %d%%", kk++);
            fflush(stdout);
          }
        }
	fprintf(outfile, "\n");
      }
      fclose(outfile);
      outfile = NULL;
    }
    printf("\rDistances written to %s.\n", outname);
    strcpy(outname_end, ".dist.id");
    outfile = fopen(outname, "w");
    if (!outfile) {
      printf("Error: Failed to open %s.\n", outname);
      goto wdist_ret_2;
    }
    for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
      if (!excluded(indiv_exclude, ii)) {
        fprintf(outfile, "%s\n", &(person_id[ii * max_person_id_len]));
      }
    }
    fclose(outfile);
    outfile = NULL;
  }
  wkspace_reset(wkspace_mark);

  if (calculation_type & CALC_GROUPDIST) {
    collapse_phenoc(pheno_c, indiv_exclude, unfiltered_indiv_ct);
    low_ct = 0;
    high_ct = 0;
    for (ii = 0; ii < indiv_ct; ii++) {
      if (pheno_c[ii] == 1) {
	high_ct++;
      } else if (pheno_c[ii] == 0) {
	low_ct++;
      }
    }
    ll_size = (low_ct * (low_ct - 1)) / 2;
    lh_size = low_ct * high_ct;
    hh_size = (high_ct * (high_ct - 1)) / 2;
    reg_tot_y = 0.0;
    reg_tot_xy = 0.0;
    reg_tot_x = 0.0;
    dhh_ssq = 0.0;
    dhl_ssq = 0.0;
    dll_ssq = 0.0;
    dist_ptr = dists;
    ll_pool = (double*)wkspace_alloc(ll_size * sizeof(double));
    lh_pool = (double*)wkspace_alloc(lh_size * sizeof(double));
    hh_pool = (double*)wkspace_alloc(hh_size * sizeof(double));
    ll_poolp = ll_pool;
    lh_poolp = lh_pool;
    hh_poolp = hh_pool;
    ped_geno = wkspace_base; // make this an explicit allocation
    for (ii = 1; ii < indiv_ct; ii++) {
      cptr = pheno_c;
      cptr2 = &(pheno_c[ii]);
      if (*cptr2 == 1) {
	while (cptr < cptr2) {
	  if (*cptr == 1) {
            dxx = *dist_ptr;
	    *hh_poolp++ = dxx;
	    reg_tot_x += dxx;
            dhh_ssq += dxx * dxx;
	  } else if (!(*cptr)) {
            dxx = *dist_ptr;
	    *lh_poolp++ = dxx;
	    reg_tot_xy += dxx;
            dhl_ssq += dxx * dxx;
	  }
	  cptr++;
	  dist_ptr++;
	}
      } else if (*cptr2 == 0) {
	while (cptr < cptr2) {
	  if (*cptr == 1) {
            dxx = *dist_ptr;
	    *lh_poolp++ = dxx;
	    reg_tot_xy += dxx;
            dhl_ssq += dxx * dxx;
	  } else if (!(*cptr)) {
            dxx = *dist_ptr;
	    *ll_poolp++ = dxx;
	    reg_tot_y += dxx;
            dll_ssq += dxx * dxx;
	  }
	  cptr++;
	  dist_ptr++;
	}
      } else {
	dist_ptr += ii;
      }
    }
    qsort(ll_pool, ll_size, sizeof(double), double_cmp);
    qsort(lh_pool, lh_size, sizeof(double), double_cmp);
    qsort(hh_pool, hh_size, sizeof(double), double_cmp);
    ll_med = get_dmedian(ll_pool, ll_size);
    lh_med = get_dmedian(lh_pool, lh_size);
    hh_med = get_dmedian(hh_pool, hh_size);
    printf("Group distance analysis (%d affected, %d unaffected):\n", high_ct, low_ct);
    if (high_ct < 2) {
      dxx = 0.0;
      dhh_sd = 0.0;
    } else {
      dww = (double)((high_ct * (high_ct - 1)) / 2);
      dxx = reg_tot_x / dww;
      dhh_sd = sqrt((dhh_ssq / dww - dxx * dxx) / (dww - 1.0));
    }
    if (!(high_ct * low_ct)) {
      dyy = 0.0;
      dhl_sd = 0.0;
    } else {
      dww = (double)(high_ct * low_ct);
      dyy = reg_tot_xy / dww;
      dhl_sd = sqrt((dhl_ssq / dww - dyy * dyy) / (dww - 1.0));
    }
    if (low_ct < 2) {
      dzz = 0.0;
      dll_sd = 0.0;
    } else {
      dww = (double)((low_ct * (low_ct - 1)) / 2);
      dzz = reg_tot_y / dww;
      dll_sd = sqrt((dll_ssq / dww - dzz * dzz) / (dww - 1.0));
    }
    printf("  Mean (sd), median dists between 2x affected     : %g (%g), %g\n", dxx, dhh_sd, hh_med);
    printf("  Mean (sd), median dists between aff. and unaff. : %g (%g), %g\n", dyy, dhl_sd, lh_med);
    printf("  Mean (sd), median dists between 2x unaffected   : %g (%g), %g\n\n", dzz, dll_sd, ll_med);
    if (2 * jackknife_d >= high_ct + low_ct) {
      printf("Delete-d jackknife skipped because d is too large.\n");
    } else {
      jackknife_iters = (groupdist_iters + thread_ct - 1) / thread_ct;

      // IN PROGRESS
      if (groupdist_d) {
	jackknife_d = groupdist_d;
      } else {
	jackknife_d = (int)pow((double)(high_ct + low_ct), 0.600000000001);
	printf("Setting d=%d for jackknife.\n", jackknife_d);
      }

      for (ulii = 1; ulii < thread_ct; ulii++) {
	if (pthread_create(&(threads[ulii - 1]), NULL, &groupdist_jack_thread, (void*)ulii)) {
	  printf("Error: Could not create thread.\n");
	  retval = RET_THREAD_CREATE_FAIL;
	  goto wdist_ret_2;
	}
      }
      ulii = 0;
      groupdist_jack_thread((void*)ulii);
      for (ii = 1; ii < thread_ct; ii++) {
        pthread_join(threads[ii - 1], NULL);
	calc_result[0] += calc_result[ii];
	calc_result2[0] += calc_result2[ii];
	calc_result3[0] += calc_result3[ii];
	calc_result4[0] += calc_result4[ii];
	calc_result5[0] += calc_result5[ii];
	calc_result6[0] += calc_result6[ii];
	calc_result7[0] += calc_result7[ii];
	calc_result8[0] += calc_result8[ii];
	calc_result9[0] += calc_result9[ii];
      }
      dxx = 1.0 / thread_ct;
      calc_result[0] *= dxx;
      calc_result2[0] *= dxx;
      calc_result3[0] *= dxx;
      dxx /= (jackknife_iters - 1) * thread_ct;
      calc_result4[0] *= dxx;
      calc_result5[0] *= dxx;
      calc_result6[0] *= dxx;
      calc_result7[0] *= dxx;
      calc_result8[0] *= dxx;
      calc_result9[0] *= dxx;
      printf("\r  AA - AU avg distance (avg s.e.): %g (%g)\n", calc_result[0] - calc_result2[0], sqrt(((high_ct + low_ct) / jackknife_d) * (calc_result4[0] + calc_result5[0] - 2 * calc_result7[0])));
      printf("  AA - UU avg distance (avg s.e.): %g (%g)\n", calc_result[0] - calc_result3[0], sqrt(((high_ct + low_ct) / jackknife_d) * (calc_result4[0] + calc_result6[0] - 2 * calc_result8[0])));
      printf("  AU - UU avg distance (avg s.e.): %g (%g)\n", calc_result2[0] - calc_result3[0], sqrt(((high_ct + low_ct) / jackknife_d) * (calc_result5[0] + calc_result6[0] - 2 * calc_result9[0])));
    }
  }
  if (calculation_type & CALC_REGRESS_DISTANCE) {
    // beta = (mean(xy) - mean(x)*mean(y)) / (mean(x^2) - mean(x)^2)
    if (unfiltered_indiv_ct != indiv_ct) {
      collapse_phenod(pheno_d, indiv_exclude, unfiltered_indiv_ct);
    }
    reg_tot_x = 0.0;
    reg_tot_xx = 0.0;
    for (ii = 0; ii < indiv_ct; ii++) {
      dxx = pheno_d[ii];
      reg_tot_x += dxx;
      reg_tot_xx += dxx * dxx;
    }
    printf("Phenotype stdev: %g\n", sqrt((reg_tot_xx - reg_tot_x * reg_tot_x / indiv_ct) / (indiv_ct - 1)));
    ulii = indiv_ct;
    ulii = ulii * (ulii - 1) / 2;
    reg_tot_xy = 0.0;
    reg_tot_x = 0.0;
    reg_tot_y = 0.0;
    reg_tot_xx = 0.0;
    dptr4 = dists;
    dptr3 = &(pheno_d[indiv_ct]);
    dist_ptr = pheno_d;
    while (++dist_ptr < dptr3) {
      dzz = *dist_ptr;
      dptr2 = pheno_d;
      while (dptr2 < dist_ptr) {
	dxx = (dzz + *dptr2++) * 0.5;
	dyy = (*dptr4++);
	reg_tot_xy += dxx * dyy;
	reg_tot_x += dxx;
	reg_tot_y += dyy;
	reg_tot_xx += dxx * dxx;
      }
    }

    dxx = ulii;
    printf("Regression slope (y = genetic distance, x = phenotype): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y / dxx) / (reg_tot_xx - reg_tot_x * reg_tot_x / dxx));

    jackknife_iters = (regress_iters + thread_ct - 1) / thread_ct;
    if (regress_d) {
      jackknife_d = regress_d;
    } else {
      jackknife_d = (int)pow((double)indiv_ct, 0.600000000001);
      printf("Setting d=%d for jackknife.\n", jackknife_d);
    }
    ped_geno = wkspace_base; // make this an explicit allocation
    for (ulii = 1; ulii < thread_ct; ulii++) {
      if (pthread_create(&(threads[ulii - 1]), NULL, &regress_jack_thread, (void*)ulii)) {
	printf("Error: Could not create thread.\n");
	retval = RET_THREAD_CREATE_FAIL;
	goto wdist_ret_2;
      }
    }
    ulii = 0;
    regress_jack_thread((void*)ulii);
    dyy = calc_result[0]; // sum
    dzz = calc_result2[0]; // sum of squares
    for (ii = 0; ii < thread_ct - 1; ii++) {
      pthread_join(threads[ii], NULL);
      dyy += calc_result[ii + 1];
      dzz += calc_result2[ii + 1];
    }
    regress_iters = jackknife_iters * thread_ct;
    printf("\rJackknife s.e.: %g\n", sqrt((indiv_ct / jackknife_d) * (dzz - dyy * dyy / regress_iters) / (regress_iters - 1)));
  }
  retval = RET_SUCCESS;
 wdist_ret_2:
  if (gz_outfile) {
    gzclose(gz_outfile);
  }
  if (person_missing_cts) {
    free(person_missing_cts);
  }
  if (marker_alleles_tmp) {
    free(marker_alleles_tmp);
  }
  if (bedtmpfile) {
    fclose(bedtmpfile);
  }
  if (bimtmpfile) {
    fclose(bimtmpfile);
  }
  if (famtmpfile) {
    fclose(famtmpfile);
  }
 wdist_ret_1:
  if (het_probs) {
    free(het_probs);
  }
  if (person_id) {
    free(person_id);
  }
  if (indiv_exclude) {
    free(indiv_exclude);
  }
  if (pheno_d) {
    free(pheno_d);
  }
  if (pheno_c) {
    free(pheno_c);
  }
  if (phenor_d) {
    free(phenor_d);
  }
  if (phenor_c) {
    free(phenor_c);
  }
  if (marker_allele_cts) {
    free(marker_allele_cts);
  }
  if (marker_alleles) {
    free(marker_alleles);
  }
  if (line_locs) {
    free(line_locs);
  }
  if (id_buf) {
    free(id_buf);
  }
  if (id_list) {
    free(id_list);
  }
  if (pid_list) {
    free(pid_list);
  }
  if (marker_exclude) {
    free(marker_exclude);
  }
  if (pedbuf) {
    free(pedbuf);
  }
  // if (marker_pos) {
  //   free(marker_pos);
  // }
  if (mafs) {
    free(mafs);
  }
  if (marker_id) {
    free(marker_id);
  }
  if (marker_chrom) {
    free(marker_chrom);
  }
 wdist_ret_0:
  if (outfile) {
    fclose(outfile);
  }
  if (freqfile) {
    fclose(freqfile);
  }
  if (filterfile) {
    fclose(filterfile);
  }
  if (phenofile) {
    fclose(phenofile);
  }
  if (famfile) {
    fclose(famfile);
  }
  if (mapfile) {
    fclose(mapfile);
  }
  if (pedfile) {
    fclose(pedfile);
  }
  return retval;
}

int param_count(int argc, char** argv, int flag_idx) {
  // Counts the number of optional parameters given to the flag at position
  // flag_idx, treating any parameter not beginning with "--" as optional.
  int opt_params = 0;
  int cur_idx = flag_idx + 1;
  while (cur_idx < argc) {
    if (argv[cur_idx][0] == '-') {
      if (argv[cur_idx][1] == '-') {
        break;
      }
    }
    opt_params++;
    cur_idx++;
  }
  return opt_params;
}

int main(int argc, char** argv) {
  unsigned char* wkspace_ua;
  char mapname[FNAMESIZE];
  char pedname[FNAMESIZE];
  char famname[FNAMESIZE];
  char outname[FNAMESIZE];
  char phenoname[FNAMESIZE];
  char extractname[FNAMESIZE];
  char excludename[FNAMESIZE];
  char filtername[FNAMESIZE];
  char freqname[FNAMESIZE];
  char* makepheno_str = NULL;
  char* filterval = NULL;
  char* argptr;
  char* sptr;
  int retval;
  int load_params = 0; // describes what file parameters have been provided
  int make_bed = 1; // STOPGAP
  int ped_col_1 = 1;
  int ped_col_34 = 1;
  int ped_col_5 = 1;
  int ped_col_6 = 1;
  int ped_col_7 = 0;
  int mpheno_col = 0;
  char* phenoname_str = NULL;
  int affection_01 = 0;
  double exponent = 0.0;
  double min_maf = 0.0;
  double geno_thresh = 1.0;
  double mind_thresh = 1.0;
  double hwe_thresh = 0.0;
  double grm_cutoff = 0.025;
  int cur_arg = 1;
  int calculation_type = 0;
  char* bubble;
  int filter_type = 0;
  int mfilter_col = 0;
  int tail_pheno = 0;
  int prune = 0;
  int missing_pheno = -9;
  unsigned char missing_geno = '0';
  double tail_bottom;
  double tail_top;
  int groupdist_iters = ITERS_DEFAULT;
  int groupdist_d = 0;
  int regress_iters = ITERS_DEFAULT;
  int regress_d = 0;
  double unrelated_herit_tol = 0.0000001;
  double unrelated_herit_covg = 0.45;
  double unrelated_herit_cove = 0.55;
  int ibc_type = 0;
  int ii;
  int jj;
  int kk;
  unsigned int chrom_mask = DEFAULT_CHROM_MASK;
  unsigned long int rseed = 0;
  FILE* scriptfile;
  int num_params;
  int in_param;
  thread_ct = sysconf(_SC_NPROCESSORS_ONLN);
  if (thread_ct == -1) {
    thread_ct = 1;
  } else if (thread_ct > MAX_THREADS) {
    thread_ct = MAX_THREADS;
  }
  strcpy(mapname, "wdist.map");
  strcpy(pedname, "wdist.ped");
  famname[0] = '\0';
  phenoname[0] = '\0';
  extractname[0] = '\0';
  excludename[0] = '\0';
  filtername[0] = '\0';
  freqname[0] = '\0';
  strcpy(outname, "wdist");
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (!strcmp(argptr, "--script")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --script parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (subst_argv) {
        printf("Error: Duplicate --script flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      scriptfile = fopen(argv[cur_arg + 1], "rb");
      if (!scriptfile) {
        printf("Error: Failed to open %s.\n", argv[cur_arg + 1]);
        return dispmsg(RET_OPENFAIL);
      }
      ii = fread(tbuf, 1, MAXLINELEN, scriptfile);
      fclose(scriptfile);
      if (ii == MAXLINELEN) {
        printf("Error: Script file too long (max %d bytes).\n", MAXLINELEN - 1);
        return dispmsg(RET_INVALID_FORMAT);
      }
      num_params = 0;
      in_param = 0;
      for (jj = 0; jj < ii; jj++) {
        if (is_space_or_eoln(tbuf[jj])) {
          if (in_param) {
            in_param = 0;
          }
        } else if (!in_param) {
          num_params += 1;
          in_param = 1;
        }
      }
      subst_argv = (char**)malloc((num_params + argc - cur_arg - 2) * sizeof(char*));
      num_params = 0;
      in_param = 0;
      for (jj = 0; jj < ii; jj++) {
        if (is_space_or_eoln(tbuf[jj])) {
          if (in_param) {
            tbuf[jj] = '\0';
            in_param = 0;
          }
        } else if (!in_param) {
          subst_argv[num_params++] = &(tbuf[jj]);
          in_param = 1;
        }
      }
      for (ii = cur_arg + 2; ii < argc; ii++) {
        subst_argv[num_params++] = argv[ii];
      }
      argc = num_params;
      cur_arg = 0;
      argv = subst_argv;
    } else if (!strcmp(argptr, "--file")) {
      if (load_params & 121) {
        if (load_params & 1) {
          printf("Error: Duplicate --file flag.\n");
        } else {
          printf("Error: --file flag cannot coexist with binary file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 1;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --file parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 5)) {
        printf("Error: --file parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      if (!(load_params & 2)) {
	strcpy(pedname, argv[cur_arg + 1]);
	strcat(pedname, ".ped");
      }
      if (!(load_params & 4)) {
	strcpy(mapname, argv[cur_arg + 1]);
	strcat(mapname, ".map");
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--help")) {
      return dispmsg(RET_MOREHELP);
    } else if (!strcmp(argptr, "--ped")) {
      if (load_params & 122) {
        if (load_params & 2) {
          printf("Error: Duplicate --ped flag.\n");
        } else {
          printf("Error: --ped flag cannot coexist with binary file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 2;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --ped parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --ped parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(pedname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--map")) {
      if (load_params & 124) {
        if (load_params & 4) {
          printf("Error: Duplicate --map flag.\n");
        } else {
          printf("Error: --map flag cannot coexist with binary file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 4;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --map parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --map parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(mapname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--make-bed")) {
      if (load_params & 120) {
        printf("Error: --make-bed cannot coexist with binary file flags.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      make_bed = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--bfile")) {
      if (load_params & 15) {
        if (load_params & 8) {
          printf("Error: Duplicate --bfile flag.\n");
        } else {
          printf("Error: --bfile flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 8;
      ii = param_count(argc, argv, cur_arg);
      if (ii > 1) {
        printf("Error: Too many --bfile parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
        sptr = argv[cur_arg + 1];
        if (strlen(sptr) > (FNAMESIZE - 5)) {
          printf("Error: --bfile parameter too long.\n");
          return dispmsg(RET_OPENFAIL);
        }
      } else {
        sptr = "wdist";
      }
      make_bed = 0;
      if (!(load_params & 16)) {
	strcpy(pedname, sptr);
	strcat(pedname, ".bed");
      }
      if (!(load_params & 32)) {
	strcpy(mapname, sptr);
	strcat(mapname, ".bim");
      }
      if (!(load_params & 64)) {
        strcpy(famname, sptr);
        strcat(famname, ".fam");
      }
      cur_arg += ii + 1;
    } else if (!strcmp(argptr, "--bed")) {
      if (load_params & 23) {
        if (load_params & 16) {
          printf("Error: Duplicate --bed flag.\n");
        } else {
          printf("Error: --bed flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 16;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --bed parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --bed parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      make_bed = 0;
      strcpy(pedname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--bim")) {
      if (load_params & 39) {
        if (load_params & 32) {
          printf("Error: Duplicate --bim flag.\n");
        } else {
          printf("Error: --bim flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 32;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --bim parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --bim parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      make_bed = 0;
      strcpy(mapname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--fam")) {
      if (load_params & 71) {
        if (load_params & 64) {
          printf("Error: Duplicate --fam flag.\n");
        } else {
          printf("Error: --fam flag cannot coexist with text input file flags.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      load_params |= 64;
      if (cur_arg == argc - 1) {
        printf("Error: Missing --fam parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --fam parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      make_bed = 0;
      strcpy(famname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--no-fid")) {
      ped_col_1 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-parents")) {
      ped_col_34 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-sex")) {
      ped_col_5 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--no-pheno")) {
      ped_col_6 = 0;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--liability")) {
      ped_col_7 = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--missing-genotype")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --missing-genotype parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (missing_geno != '0') {
        printf("Error: Duplicate --missing-genotype flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      missing_geno = (unsigned char)argv[cur_arg + 1][0];
      if ((strlen(argv[cur_arg + 1]) > 1) || (missing_geno <= ' ')) {
        printf("Error: Invalid --missing-genotype parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--missing-phenotype")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --missing-phenotype parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (missing_pheno != -9) {
        printf("Error: Duplicate --missing-phenotype flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      missing_pheno = atoi(argv[cur_arg + 1]);
      if ((missing_pheno == 0) || (missing_pheno == 1)) {
        printf("Error: Invalid --missing-phenotype parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--pheno")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --pheno parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname[0]) {
        if (makepheno_str) {
          printf("Error: --pheno and --make-pheno flags cannot coexist.\n");
        } else {
          printf("Error: Duplicate --pheno flag.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --pheno parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(phenoname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--make-pheno")) {
      if (cur_arg > argc - 2) {
        printf("Error: Not enough --make-pheno parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname[0]) {
        if (makepheno_str) {
          printf("Error: Duplicate --make-pheno flag.\n");
        } else {
          printf("Error: --pheno and --make-pheno flags cannot coexist.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (tail_pheno) {
        printf("Error: --make-pheno and --tail-pheno flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 1)) {
        printf("Error: --make-pheno parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(phenoname, argv[cur_arg + 1]);
      makepheno_str = argv[cur_arg + 2];
      cur_arg += 3;
    } else if (!strcmp(argptr, "--mpheno")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --mpheno parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (mpheno_col != 0) {
        printf("Error: Duplicate --mpheno flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname_str) {
        printf("Error: --mpheno and --pheno-name flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      mpheno_col = atoi(argv[cur_arg + 1]);
      if (mpheno_col < 1) {
        printf("Error: Invalid --mpheno parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--pheno-name")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --pheno-name parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (mpheno_col != 0) {
        printf("Error: --mpheno and --pheno-name flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (phenoname_str) {
        printf("Error: Duplicate --pheno-name flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      phenoname_str = argv[cur_arg + 1];
      cur_arg += 2;
    } else if (!strcmp(argptr, "--prune")) {
      prune = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--1")) {
      affection_01 = 1;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--tail-pheno")) {
      ii = param_count(argc, argv, cur_arg);
      if (ii == 0) {
        printf("Error: Not enough --tail-pheno parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 2) {
        printf("Error: Too many --tail-pheno paramaters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (tail_pheno) {
        printf("Error: Duplicate --tail-pheno flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (makepheno_str) {
        printf("Error: --make-pheno and --tail-pheno flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lg", &tail_bottom) != 1) {
        printf("Error: Invalid --tail-pheno parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii == 1) {
        tail_top = tail_bottom;
      } else {
	if (sscanf(argv[cur_arg + 2], "%lg", &tail_top) != 1) {
	  printf("Error: Invalid --tail-pheno parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
      }
      if (tail_bottom > tail_top) {
        printf("Error: Ltop cannot be larger than Hbottom for --tail-pheno.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      tail_pheno = 1;
      cur_arg += 3;
    } else if (!strcmp(argptr, "--extract")) {
      if (extractname[0]) {
        printf("Error: Duplicate --extract flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii == 0) {
        printf("Error: Missing --extract parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 1) {
        printf("Error: Too many --extract parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--extract filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(extractname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--exclude")) {
      if (excludename[0]) {
        printf("Error: Duplicate --exclude flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii == 0) {
        printf("Error: Missing --exclude parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 1) {
        printf("Error: Too many --exclude parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--exclude filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(excludename, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--write-snplist")) {
      if (calculation_type & CALC_WRITE_SNPLIST) {
        printf("Error: Duplicate --write-snplist flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      calculation_type |= CALC_WRITE_SNPLIST;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--keep")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --keep parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (filter_type) {
        if (filter_type == 1) {
          printf("Error: Duplicate --keep flag.\n");
        } else {
          printf("Error: --keep + --remove/--filter not supported.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = FILTER_KEEP;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--keep filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--remove")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --remove parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (filter_type) {
        if (filter_type == 2) {
          printf("Error: Duplicate --remove flag.\n");
        } else {
          printf("Error: --remove + --keep/--filter not supported.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = FILTER_REMOVE;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--remove filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--filter")) {
      if (cur_arg > argc - 2) {
        printf("Error: Not enough --filter parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (filter_type) {
        if (filter_type == 3) {
          printf("Error: Duplicate --filter flag.\n");
        } else {
          printf("Error: --filter + --keep/--remove not supported.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      filter_type = FILTER_CUSTOM;
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("Error: --filter filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(filtername, argv[cur_arg + 1]);
      filterval = argv[cur_arg + 2];
      cur_arg += 3;
    } else if (!strcmp(argptr, "--mfilter")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --mfilter parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (mfilter_col) {
        printf("Error: Duplicate --mfilter flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      mfilter_col = atoi(argv[cur_arg + 1]);
      if (mfilter_col < 1) {
        printf("Error: Invalid --mfilter parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--memory")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --memory parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      malloc_size_mb = atoi(argv[cur_arg + 1]);
      if (malloc_size_mb < 64) {
        printf("Error: Invalid --memory parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--threads")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --threads parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) {
        printf("Error: Invalid --threads parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > MAX_THREADS) {
        printf("Error: --threads parameter too large (max %d).\n", MAX_THREADS);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      thread_ct = ii;
      cur_arg += 2;
    } else if (!strcmp(argptr, "--exponent")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --exponent parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (sscanf(argv[cur_arg + 1], "%lg", &exponent) != 1) {
        printf("Error: Invalid --exponent parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--autosome")) {
      if (chrom_mask != DEFAULT_CHROM_MASK) {
        if (chrom_mask == 0x007ffffe) {
          printf("Error: Duplicate --autosome flag.\n");
        } else {
          printf("Error: --chr and --autosome flags cannot coexist.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      chrom_mask = 0x007ffffe;
      cur_arg += 1;
    } else if (!strcmp(argptr, "--chr")) {
      ii = param_count(argc, argv, cur_arg);
      if (chrom_mask != DEFAULT_CHROM_MASK) {
        if (chrom_mask == 0x007ffffe) {
          printf("Error: --chr and --autosome flags cannot coexist.\n");
        } else {
          printf("Error: Duplicate --chr flag.\n");
        }
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (!ii) {
        printf("Error: Missing --chr parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      chrom_mask = 0;
      for (jj = 0; jj < ii; jj++) {
        kk = marker_code(argv[cur_arg + 1 + jj]);
        // will allow some non-autosomal in future
        if ((kk == -1) || ((kk > 22) && (kk != 25))) {
          printf("Error: Invalid --chr parameter.\n");
          return dispmsg(RET_INVALID_CMDLINE);
        }
        chrom_mask |= 1 << kk;
      }
      cur_arg += 1 + ii;
    } else if (!strcmp(argptr, "--maf")) {
      ii = param_count(argc, argv, cur_arg);
      if (ii > 1) {
        printf("Error: Too many --maf parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
        if (sscanf(argv[cur_arg + 1], "%lg", &min_maf) != 1) {
          printf("Error: Invalid --maf parameter.\n");
          return dispmsg(RET_INVALID_CMDLINE);
        }
	if ((min_maf <= 0.0) || (min_maf > 0.5)) {
	  printf("Error: Invalid --maf parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
      } else {
        min_maf = 0.01;
      }
      cur_arg += ii + 1;
    } else if (!strcmp(argptr, "--geno")) {
      ii = param_count(argc, argv, cur_arg);
      if (ii > 1) {
        printf("Error: Too many --geno parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
        if (sscanf(argv[cur_arg + 1], "%lg", &geno_thresh) != 1) {
          printf("Error: Invalid --geno parameter.\n");
          return dispmsg(RET_INVALID_CMDLINE);
        }
	if ((geno_thresh < 0.0) || (geno_thresh > 1.0)) {
	  printf("Error: Invalid --geno parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
      } else {
        geno_thresh = 0.1;
      }
      cur_arg += ii + 1;
    } else if (!strcmp(argptr, "--mind")) {
      ii = param_count(argc, argv, cur_arg);
      if (ii > 1) {
        printf("Error: Too many --mind parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
	if (sscanf(argv[cur_arg + 1], "%lg", &mind_thresh) != 1) {
	  printf("Error: Invalid --mind parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
	if ((mind_thresh < 0.0) || (mind_thresh > 1.0)) {
	  printf("Error: Invalid --mind parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
      } else {
        mind_thresh = 0.1;
      }
      cur_arg += ii + 1;
    } else if (!strcmp(argptr, "--hwe")) {
      ii = param_count(argc, argv, cur_arg);
      if (ii > 1) {
        printf("Error: Too many --hwe parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
	if (sscanf(argv[cur_arg + 1], "%lg", &hwe_thresh) != 1) {
	  printf("Error: Invalid --hwe parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
	if ((hwe_thresh < 0.0) || (hwe_thresh >= 1.0)) {
	  printf("Error: Invalid --hwe parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
      } else {
        hwe_thresh = 0.001;
      }
      cur_arg += ii + 1;
    } else if ((!strcmp(argptr, "--rel-cutoff")) || (!strcmp(argptr, "--grm-cutoff"))) {
      ii = param_count(argc, argv, cur_arg);
      if (calculation_type & CALC_REL_CUTOFF) {
        printf("Error: Duplicate --rel-cutoff flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii > 1) {
        printf("Error: Too many --rel-cutoff parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
        if (sscanf(argv[cur_arg + 1], "%lg", &grm_cutoff) != 1) {
          printf("Error: Invalid --rel-cutoff parameter.\n");
          return dispmsg(RET_INVALID_CMDLINE);
        }
        if ((grm_cutoff <= 0.0) || (grm_cutoff >= 1.0)) {
          printf("Error: Invalid --rel-cutoff parameter.\n");
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      calculation_type |= CALC_REL_CUTOFF;
      cur_arg += ii + 1;
    } else if (!strcmp(argptr, "--rseed")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --rseed parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (rseed != 0) {
        printf("Error: Duplicate --rseed flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      rseed = (unsigned long int)atoi(argv[cur_arg + 1]);
      if (rseed == 0) {
        printf("Error: Invalid --rseed parameter.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 2;
    } else if (!strcmp(argptr, "--out")) {
      if (cur_arg == argc - 1) {
        printf("Error: Missing --out parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - 10)) {
        printf("Error: --out parameter too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(outname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--silent")) {
      freopen("/dev/null", "w", stdout);
      cur_arg += 1;
    } else if (!strcmp(argptr, "--make-rel")) {
      if (calculation_type & CALC_RELATIONSHIP_MASK) {
        printf("Error: --make-rel cannot coexist with another relationship matrix file\ncreation flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 4) {
        printf("Error: --make-rel accepts at most 4 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      jj = 0;
      for (kk = 1; kk <= ii; kk++) {
        if (!strcmp(argv[cur_arg + kk], "no-var-std")) {
          if (calculation_type & CALC_IBC) {
            printf("Error: --make-rel 'no-var-std' modifier cannot coexist with --ibc flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (calculation_type & CALC_UNRELATED_HERITABILITY) {
            printf("Error: --make-rel 'no-var-std' modifier cannot coexist with\n--unrelated-heritability flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-rel 'no-var-std' modifier cannot coexist with an IBC modifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_NO_VAR_STD;
        } else if (!strcmp(argv[cur_arg + kk], "gz")) {
          if (calculation_type & CALC_RELATIONSHIP_BIN) {
            printf("Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_GZ;
        } else if (!strcmp(argv[cur_arg + kk], "bin")) {
          if (calculation_type & CALC_RELATIONSHIP_GZ) {
            printf("Error: --make-rel 'gz' and 'bin' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_BIN;
        } else if (!strcmp(argv[cur_arg + kk], "square")) {
          if (calculation_type & CALC_RELATIONSHIP_SQ0) {
            printf("Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_SQ;
        } else if (!strcmp(argv[cur_arg + kk], "square0")) {
          if (calculation_type & CALC_RELATIONSHIP_SQ) {
            printf("Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_SQ0;
        } else if ((!strcmp(argv[cur_arg + kk], "ibc1")) || (!strcmp(argv[cur_arg + kk], "ibc2"))) {
          if (calculation_type & CALC_RELATIONSHIP_NO_VAR_STD) {
            printf("Error: --make-rel 'no-var-std' modifier cannot coexist with an IBC modifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-rel '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + kk], errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          ibc_type = argv[cur_arg + kk][3] - '0';
        } else {
          printf("Error: Invalid --make-rel parameter.%s", errstr_append);
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      cur_arg += ii + 1;
      if (!(calculation_type & CALC_RELATIONSHIP_MASK)) {
        calculation_type |= CALC_RELATIONSHIP_NULL;
      }
    } else if (!strcmp(argptr, "--make-grm")) {
      if (calculation_type & CALC_RELATIONSHIP_MASK) {
        printf("Error: --make-grm cannot coexist with another relationship matrix file\ncreation flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 3) {
        printf("Error: --make-grm accepts at most 3 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      calculation_type |= CALC_RELATIONSHIP_GZ;
      for (kk = 1; kk <= ii; kk++) {
        if (!strcmp(argv[cur_arg + kk], "no-var-std")) {
          if (calculation_type & CALC_IBC) {
            printf("Error: --make-grm 'no-var-std' modifier cannot coexist with --ibc flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (calculation_type & CALC_UNRELATED_HERITABILITY) {
            printf("Error: --make-grm 'no-var-std' modifier cannot coexist with\n--unrelated-heritability flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-grm 'no-var-std' modifier cannot coexist with an IBC modifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_NO_VAR_STD;
        } else if (!strcmp(argv[cur_arg + kk], "no-gz")) {
          calculation_type &= ~CALC_RELATIONSHIP_GZ;
        } else if ((!strcmp(argv[cur_arg + kk], "ibc1")) || (!strcmp(argv[cur_arg + kk], "ibc2"))) {
          if (calculation_type & CALC_RELATIONSHIP_NO_VAR_STD) {
            printf("Error: --make-grm 'no-var-std' modifier cannot coexist with an IBC modifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-grm '%s' modifier cannot coexist with another IBC modifier.%s", argv[cur_arg + kk], errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          ibc_type = argv[cur_arg + kk][3] - '0';
        } else {
          printf("Error: Invalid --make-grm parameter.%s", errstr_append);
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_RELATIONSHIP_GRM;
    } else if (!strcmp(argptr, "--make-grm-bin")) {
      if (calculation_type & CALC_RELATIONSHIP_MASK) {
        printf("Error: --make-grm-bin cannot coexist with another relationship matrix file\ncreation flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 3) {
        printf("Error: --make-grm-bin accepts at most 3 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      for (kk = 1; kk <= ii; kk++) {
        if (!strcmp(argv[cur_arg + kk], "no-var-std")) {
          if (calculation_type & CALC_IBC) {
            printf("Error: --make-grm-bin 'no-var-std' modifier cannot coexist with --ibc flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (calculation_type & CALC_UNRELATED_HERITABILITY) {
            printf("Error: --make-grm-bin 'no-var-std' modifier cannot coexist with\n--unrelated-heritability flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-grm-bin 'no-var-std' modifier cannot coexist with an IBC\nmodifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_NO_VAR_STD;
        } else if (!strcmp(argv[cur_arg + kk], "square")) {
          if (calculation_type & CALC_RELATIONSHIP_SQ0) {
            printf("Error: --make-grm-bin 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_SQ;
        } else if (!strcmp(argv[cur_arg + kk], "square0")) {
          if (calculation_type & CALC_RELATIONSHIP_SQ) {
            printf("Error: --make-grm-bin 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_SQ0;
        } else if ((!strcmp(argv[cur_arg + kk], "ibc1")) || (!strcmp(argv[cur_arg + kk], "ibc2"))) {
          if (calculation_type & CALC_RELATIONSHIP_NO_VAR_STD) {
            printf("Error: --make-grm-bin 'no-var-std' modifier cannot coexist with an IBC\nmodifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-grm-bin '%s' modifier cannot coexist with another IBC\nmodifier.%s", argv[cur_arg + kk], errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          ibc_type = argv[cur_arg + kk][3] - '0';
        } else {
          printf("Error: Invalid --make-grm-bin parameter.%s", errstr_append);
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_RELATIONSHIP_GRM | CALC_RELATIONSHIP_BIN;
    } else if (!strcmp(argptr, "--ibc")) {
      if (calculation_type & CALC_RELATIONSHIP_NO_VAR_STD) {
        printf("Error: --ibc flag cannot coexist with a 'no-var-std' relationship matrix\ncalculation.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 1;
      calculation_type |= CALC_IBC;
    } else if (!strcmp(argptr, "--distance")) {
      if (calculation_type & CALC_DISTANCE_MASK) {
        printf("Error: --distance cannot coexist with another distance matrix file creation\nflag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 2) {
        printf("Error: --distance accepts at most 2 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      jj = 0;
      for (kk = 1; kk <= ii; kk++) {
        if (!strcmp(argv[cur_arg + kk], "square")) {
          if (calculation_type & CALC_DISTANCE_SQ0) {
            printf("Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_DISTANCE_SQ;
        } else if (!strcmp(argv[cur_arg + kk], "square0")) {
          if (calculation_type & CALC_DISTANCE_SQ) {
            printf("Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_DISTANCE_SQ0;
        } else if (!strcmp(argv[cur_arg + kk], "gz")) {
          if (calculation_type & CALC_DISTANCE_BIN) {
            printf("Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_DISTANCE_GZ;
        } else if (!strcmp(argv[cur_arg + kk], "bin")) {
          if (calculation_type & CALC_DISTANCE_GZ) {
            printf("Error: --distance 'gz' and 'bin' flags cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_DISTANCE_BIN;
        } else {
          printf("Error: Invalid --distance parameter.%s", errstr_append);
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      cur_arg += ii + 1;
      if (!(calculation_type & CALC_DISTANCE_MASK)) {
        calculation_type |= CALC_DISTANCE_NULL;
      }
    } else if (!strcmp(argptr, "--groupdist")) {
      if (calculation_type & CALC_GROUPDIST) {
        printf("Error: Duplicate --groupdist flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 2) {
        printf("Error: --groupdist accepts at most 2 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
	groupdist_iters = atoi(argv[cur_arg + 1]);
	if (groupdist_iters < 2) {
	  printf("Error: Invalid --groupdist jackknife iteration count.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
        if (ii == 2) {
	  groupdist_d = atoi(argv[cur_arg + 2]);
	  if (groupdist_d <= 0) {
	    printf("Error: Invalid --groupdist jackknife delete parameter.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_GROUPDIST;
    } else if (!strcmp(argptr, "--regress")) {
      printf("Error: --regress flag has been renamed to --regress-distance.%s", errstr_append);
      return dispmsg(RET_INVALID_CMDLINE);
    } else if (!strcmp(argptr, "--regress-distance")) {
      if (calculation_type & CALC_REGRESS_DISTANCE) {
        printf("Error: Duplicate --regress-distance flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 2) {
        printf("Error: --regress-distance accepts at most 2 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
	regress_iters = atoi(argv[cur_arg + 1]);
	if (regress_iters < 2) {
	  printf("Error: Invalid --regress-distance jackknife iteration count.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
        if (ii == 2) {
	  regress_d = atoi(argv[cur_arg + 2]);
	  if (regress_d <= 0) {
	    printf("Error: Invalid --regress-distance jackknife delete parameter.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_REGRESS_DISTANCE;
    } else if (!strcmp(argptr, "--unrelated-heritability")) {
      if (calculation_type & CALC_RELATIONSHIP_NO_VAR_STD) {
        printf("Error: --unrelated-heritability flag cannot coexist with a 'no-var-std'\nrelationship matrix calculation.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (calculation_type & CALC_UNRELATED_HERITABILITY) {
        printf("Error: Duplicate --unrelated-heritability flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 4) {
        printf("Error: --unrelated-heritability accepts at most four parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
        if (!strcmp(argv[cur_arg + 1], "strict")) {
          calculation_type |= CALC_UNRELATED_HERITABILITY_STRICT;
          jj = 2;
        } else {
          jj = 1;
        }
        if (ii >= jj) {
	  if (sscanf(argv[cur_arg + jj], "%lg", &unrelated_herit_tol) != 1) {
	    printf("Error: Invalid --unrelated-heritability EM tolerance parameter.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
	  if (unrelated_herit_tol <= 0.0) {
	    printf("Error: Invalid --unrelated-heritability EM tolerance parameter.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
	  if (ii > jj) {
	    if (sscanf(argv[cur_arg + jj + 1], "%lg", &unrelated_herit_covg) != 1) {
	      printf("Error: Invalid --unrelated-heritability genetic covariance prior.\n");
	      return dispmsg(RET_INVALID_CMDLINE);
	    }
	    if ((unrelated_herit_covg <= 0.0) || (unrelated_herit_covg > 1.0)) {
	      printf("Error: Invalid --unrelated-heritability genetic covariance prior.\n");
	      return dispmsg(RET_INVALID_CMDLINE);
	    }
	    if (ii == jj + 2) {
	      if (sscanf(argv[cur_arg + jj + 2], "%lg", &unrelated_herit_cove) != 1) {
		printf("Error: Invalid --unrelated-heritability residual covariance prior.\n");
		return dispmsg(RET_INVALID_CMDLINE);
	      }
	      if ((unrelated_herit_cove <= 0.0) || (unrelated_herit_cove > 1.0)) {
		printf("Error: Invalid --unrelated-heritability residual covariance prior.\n");
		return dispmsg(RET_INVALID_CMDLINE);
	      }
	    } else {
	      unrelated_herit_cove = 1.0 - unrelated_herit_covg;
	    }
          }
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_UNRELATED_HERITABILITY;
    } else if (!strcmp(argptr, "--freq")) {
      if (freqname[0]) {
        printf("Error: --freq and --update-freq flags cannot coexist.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (calculation_type & CALC_FREQ) {
        printf("Error: Duplicate --freq flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 1) {
        printf("Error: --freq accepts at most one parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
        if (!strcmp(argv[cur_arg + 1], "gcta")) {
          calculation_type |= CALC_FREQ_GCTA;
        } else {
          printf("Error: Invalid --freq parameter.\n");
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_FREQ;
    } else if (!strcmp(argptr, "--update-freq")) {
      if (calculation_type & CALC_FREQ) {
        printf("Error: --freq and --update-freq flags cannot coexist.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (freqname[0]) {
        printf("Error: Duplicate --update-freq flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 1) {
        printf("Error: --update-freq accepts at most one parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("Error: --update-freq filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(freqname, argv[cur_arg + 1]);
      cur_arg += 1 + ii;
    } else if (!strcmp(argptr, "--map3")) {
      printf("Note: --map3 flag unnecessary (.map file format is autodetected).\n");
      cur_arg += 1;
    } else if (!strcmp(argptr, "--compound-genotypes")) {
      printf("Note: --compound-genotypes flag unnecessary (spaces between alleles in .ped\nare optional).\n");
      cur_arg += 1;
    } else {
      printf("Error: Invalid argument (%s).%s", argv[cur_arg], errstr_append);
      return dispmsg(RET_INVALID_CMDLINE);
    }
  }
  if ((calculation_type & ~CALC_REL_CUTOFF) == 0) {
    return dispmsg(RET_HELP);
  }
  if (prune && (!phenoname[0]) && (!ped_col_6)) {
    printf("Error: --prune and --no-pheno cannot coexist without an alternate phenotype\nfile.\n");
    return dispmsg(RET_INVALID_CMDLINE);
  }
  if (subst_argv) {
    free(subst_argv);
  }
#ifndef __LP64__
  popcount[0] = 0;
  for (ii = 0; ii < 65536; ii++) {
    popcount[ii] = (ii & 1) + popcount[ii / 2];
  }
#endif

  bubble = (char*)malloc(67108864 * sizeof(char));
  if (!bubble) {
    return dispmsg(RET_NOMEM);
  }
  wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  if ((malloc_size_mb > MALLOC_DEFAULT_MB) && !wkspace_ua) {
    printf("%lld MB memory allocation failed.  Using default allocation behavior.\n", malloc_size_mb);
    malloc_size_mb = MALLOC_DEFAULT_MB;
  }
  while (!wkspace_ua) {
    if (malloc_size_mb > 128) {
      malloc_size_mb -= 64;
    } else {
      malloc_size_mb = 64;
    }
    wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (wkspace_ua) {
      printf("Allocated %lld MB successfully.\n", malloc_size_mb);
    }
  }
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((unsigned long)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = malloc_size_mb * 1048576 - (unsigned long)(wkspace - wkspace_ua);
  free(bubble);

  if (!rseed) {
    rseed = (unsigned long int)time(NULL);
  }
  init_genrand(rseed);

  // famname[0] indicates binary vs. text
  // filtername[0] indicates existence of filter
  // freqname[0] signals --update-freq
  retval = wdist(pedname, mapname, famname, phenoname, extractname, excludename, filtername, freqname, makepheno_str, filter_type, filterval, mfilter_col, make_bed, ped_col_1, ped_col_34, ped_col_5, ped_col_6, ped_col_7, (char)missing_geno, missing_pheno, mpheno_col, phenoname_str, prune, affection_01, chrom_mask, exponent, min_maf, geno_thresh, mind_thresh, hwe_thresh, grm_cutoff, tail_pheno, tail_bottom, tail_top, outname, calculation_type, groupdist_iters, groupdist_d, regress_iters, regress_d, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_cove, ibc_type);
  free(wkspace_ua);
  return dispmsg(retval);
}
