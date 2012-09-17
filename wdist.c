// WDIST weighted genomic distance calculator
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


// Uncomment this to build without CBLAS/CLAPACK:
// #define NOLAPACK


// The key ideas behind this calculator's design are:
//
// 1. Incremental processing of SNPs.  Each element A_{jk} of a distance or
// relationship matrix is a sum of M terms, one for each SNP, multiplied by a
// missingness correction at the end.  We can walk through the SNPs
// sequentially without keeping much in memory beyond partial sums;
// conveniently, this plays well with the decision made by the PLINK team a few
// years ago to switch to SNP-major binary files.
//
// 2. Multiplexing  of markers using bitwise operations.  For instance, there
// are only seven possible ways SNP i can affect the relationship matrix entry
// between individuals j and k:
//    a. j and k are both homozygous for the rare allele
//    b. one is homozygous rare, one is heterozygous
//    c. one is homozygous rare, one is homozygous common
//    d. both are heterozygous
//    e. one is heterozygous, one is homozygous common
//    f. both are homozygous common
//    g. one or both has missing genotype data
// Seven cases can be distinguished by three bits, so one can compose a
// sequence of bitwise operations that maps a pair of padded 2-bit PLINK
// genotypes to seven different 3-bit values according to this breakdown.
// On 64-bit machines, this allows integer operations to act on 20 markers
// simultaneously.  (There's space for a 21st, but we currently choose not to
// use it.)
//
// 3. Lookup tables describing the effect of 5-7 markers at a time on a
// distance or relationship, rather than just one.  For relationship matrices,
// idea #2 allows computation of a single 64-bit integer where bits 0-2
// describe the relationship on marker #0, bits 3-5 describe the relationship
// on marker #1, ..., all the way up to bits 57-59 describing the relationship
// on marker #19.  We then want to perform the update
//    A_{jk} := A_{jk} + f_0(x_0) + f_1(x_1) + ... + f_19(x_19)
// where the x_i's are bit trios, and the f_i's map them to floating point
// terms.  We could do this with 20 table lookups and floating point additions.
// Or, we could structure this as
//    A_{jk} := A_{jk} + f_{0-4}(x_{0-4}) + ... + f_{15-19}(x_{15-19})
// where x_{0-4} denotes the lowest order *15* bits, and f_{0-4} maps them
// directly to f_0(x_0) + f_1(x_1) + f_2(x_2) + f_3(x_3) + f_4(x_4); similarly
// for f_{5-9}, f_{10-14}, and f_{15-19}.  This requires precomputation of four
// lookup tables of size 2^15, total size 1 MB (which fits comfortably in a
// typical L2 cache these days), and lets you get away with four table lookups
// and adds instead of twenty.  These updates are performed N(N-1)/2 times
// before the tables need to be recalculated, so the precomputation cost is
// often recouped hundreds or thousands of times over.
//
// 4. Splitting the distance/relationship matrix into pieces of roughly equal
// size and assigning one thread to each piece.  This is an "embarrassingly
// parallel" problem with no need for careful synchronization.  Cluster
// computing is supported in essentially the same manner.
//
// We also note that IBS/zero exponent distance matrices are special cases,
// reducing to Hamming weight computations plus a bit of dressing to handle
// missing markers.  Hamming weight computation on x86 CPUs has been studied
// extensively by others; a good reference as of this writing is
// http://www.dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html .
// We use a variant of the Kim Walisch/Cedric Lauradoux bitslicing algorithm
// discussed there (with most 64-bit integer operations replaced by SSE2
// instructions), which runs several times faster than our corresponding
// nonzero exponent distance matrix computation.
//
//
// In the end, we can't get around the O(MN^2) nature of these calculations,
// but we believe we have beaten down the leading constant by a large enough
// factor to matter.

#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
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
// floating point comparison-to-nonzero tolerance, currently 2^{-30}
#define EPSILON 0.000000000931322574615478515625
#define CACHELINE 64 // assumed number of bytes per cache line, for alignment
#define CACHELINE_DBL (CACHELINE / sizeof(double))

#define RET_SUCCESS 0
#define RET_NOMEM 1
#define RET_OPENFAIL 2
#define RET_INVALID_FORMAT 3
#define RET_CALC_NOT_YET_SUPPORTED 4
#define RET_INVALID_CMDLINE 5
#define RET_WRITE_FAIL 6
#define RET_READ_FAIL 7
#define RET_MOREHELP 8
#define RET_THREAD_CREATE_FAIL 9

#define FILTER_NONE 0
#define FILTER_KEEP 1
#define FILTER_REMOVE 2

#define CALC_RELATIONSHIP_COV 1
#define CALC_RELATIONSHIP_SQ 2
#define CALC_RELATIONSHIP_SQ0 4
#define CALC_RELATIONSHIP_TRI 6
#define CALC_RELATIONSHIP_SHAPEMASK 6
#define CALC_RELATIONSHIP_GZ 8
#define CALC_RELATIONSHIP_BIN 0x10
#define CALC_RELATIONSHIP_GRM 0x20
#define CALC_RELATIONSHIP_MASK 0x3f
#define CALC_IBC 0x40
#define CALC_DISTANCE_SQ 0x80
#define CALC_DISTANCE_SQ0 0x100
#define CALC_DISTANCE_TRI 0x180
#define CALC_DISTANCE_SHAPEMASK 0x180
#define CALC_DISTANCE_GZ 0x200
#define CALC_DISTANCE_BIN 0x400
#define CALC_DISTANCE_IBS 0x800
#define CALC_DISTANCE_1_MINUS_IBS 0x1000
#define CALC_DISTANCE_SNPS 0x2000
#define CALC_DISTANCE_FORMATMASK 0x3800
#define CALC_DISTANCE_MASK 0x3f80
#define CALC_PLINK_DISTANCE_MATRIX 0x4000
#define CALC_PLINK_IBS_MATRIX 0x8000
#define CALC_GDISTANCE_MASK 0xff80
#define CALC_LOAD_DISTANCES 0x10000
#define CALC_GROUPDIST 0x20000
#define CALC_REGRESS_DISTANCE 0x40000
#define CALC_UNRELATED_HERITABILITY 0x80000
#define CALC_UNRELATED_HERITABILITY_STRICT 0x100000
#define CALC_FREQ 0x200000
#define CALC_REL_CUTOFF 0x400000
#define CALC_WRITE_SNPLIST 0x800000
#define CALC_GENOME 0x1000000
#define CALC_REGRESS_REL 0x2000000

#define _FILE_OFFSET_BITS 64
#define MAX_THREADS 63
#define MAX_THREADS_P1 64
#define PEDBUFBASE 256
#define FNAMESIZE 2048
#define MALLOC_DEFAULT_MB 2176
// size of generic text line load buffer.  .ped lines can of course be longer
#define MAXLINELEN 131072

// number of different types of jackknife values to precompute (x^2, x, y, xy)
#define JACKKNIFE_VALS_DIST 5
#define JACKKNIFE_VALS_REL 5

// allow .mdist.bin.xxxxxxxxxx extension
#define MAX_POST_EXT 22

#define MAX_EM_ACCEL 100.0

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^31 - 1
#define PARALLEL_MAX 32768

// default jackknife iterations
#define ITERS_DEFAULT 100000

// 0-22, plus XY
#define DEFAULT_CHROM_MASK 0x027fffff

#define DEFAULT_PPC_GAP 500000

// Number of snp-major .bed lines to read at once for distance calc if exponent
// is zero.  Currently assumed to be a multiple of 192, and no larger than
// 1920, by the popcount_..._multiword functions.  (The optimal value depends
// on both system-specific properties such as cache sizes, as well as the
// number of individuals in the current calculation, so in principle it's best
// to select this value at runtime.  But 960 usually works well in practice in
// my experience.)
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)

// Multiple of 192, no larger than 1728.  (Can be increased to 1920 with a
// minor dot product implementation change.)
#define MULTIPLEX_COV 768
#define MULTIPLEX_2COV (MULTIPLEX_COV * 2)

// Must be multiple of 384, no larger than 3840.
#define GENOME_MULTIPLEX 1152
#define GENOME_MULTIPLEX2 (GENOME_MULTIPLEX * 2)

#if __LP64__
#define FIVEMASK 0x5555555555555555LU
#define AAAAMASK 0xaaaaaaaaaaaaaaaaLU
#define BITCT 64
// number of snp-major .bed lines to read at once for distance calc if exponent
// is nonzero.
#define MULTIPLEX_DIST_EXP 64
// number of snp-major .bed lines to read at once for relationship calc
#define MULTIPLEX_REL 60
#else
// N.B. 32-bit version not as carefully tested or optimized, but I'll try to
// make sure it works properly
#define FIVEMASK 0x55555555
#define AAAAMASK 0xaaaaaaaa
#define BITCT 32
#define MULTIPLEX_DIST_EXP 28
#define MULTIPLEX_REL 30
#endif

#define BITCT2 (BITCT / 2)

#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - 1)))
#define CACHEALIGN_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1)))

#define MIN(aa, bb) ((aa) > (bb))? (bb) : (aa)

const char ver_str[] =
#ifdef NOLAPACK
  "WDIST genomic distance calculator, v0.9.5 "
#else
  "WDIST genomic distance calculator, v0.9.5L "
#endif
#if __LP64__
  "64-bit"
#else
  "32-bit"
#endif
  " (15 September 2012)\n"
  "(C) 2012 Christopher Chang, BGI Cognitive Genomics Lab    GNU GPL, version 3\n";
const char errstr_append[] = "\nRun 'wdist --help | more' for more information.\n";
const char errstr_fopen[] = "Error: Failed to open %s.\n";
const char errstr_map_format[] = "Error: Improperly formatted .map file.\n";
const char errstr_fam_format[] = "Error: Improperly formatted .fam file.\n";
const char errstr_ped_format[] = "Error: Improperly formatted .ped file.\n";
const char errstr_phenotype_format[] = "Error: Improperly formatted phenotype file.\n";
const char errstr_filter_format[] = "Error: Improperly formatted filter file.\n";
const char errstr_freq_format[] = "Error: Improperly formatted frequency file.\n";
// fit 4 pathologically long IDs plus other stuff
char tbuf[MAXLINELEN * 4 + 256];

int fopen_checked(FILE** target_ptr, char* fname, char* mode) {
  *target_ptr = fopen(fname, mode);
  if (!(*target_ptr)) {
    printf(errstr_fopen, fname);
    return -1;
  }
  return 0;
}

int gzopen_checked(gzFile* target_ptr, char* fname, char* mode) {
  *target_ptr = gzopen(fname, mode);
  if (!(*target_ptr)) {
    printf(errstr_fopen, fname);
    return -1;
  }
  return 0;
}

inline int fwrite_checked(void* buf, size_t len, FILE* outfile) {
  if ((!len) || fwrite(buf, len, 1, outfile)) {
    return 0;
  }
  return -1;
}

inline int gzwrite_checked(gzFile gz_outfile, void* buf, size_t len) {
  if ((!len) || gzwrite(gz_outfile, buf, len)) {
    return 0;
  }
  return -1;
}

inline int fclose_null(FILE** fptr_ptr) {
  int ii;
  ii = fclose(*fptr_ptr);
  *fptr_ptr = NULL;
  return ii;
}

inline void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

inline void gzclose_cond(gzFile gz_outfile) {
  if (gz_outfile) {
    gzclose(gz_outfile);
  }
}

inline void free_cond(void* memptr) {
  if (memptr) {
    free(memptr);
  }
}

// manually managed, very large stack
unsigned char* wkspace;
unsigned char* wkspace_base;
long long malloc_size_mb = MALLOC_DEFAULT_MB;
long long wkspace_left;

unsigned char* wkspace_alloc(unsigned long size) {
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

inline int wkspace_alloc_uc_checked(unsigned char** ucp_ptr, unsigned long size) {
  *ucp_ptr = wkspace_alloc(size);
  if (!(*ucp_ptr)) {
    return 1;
  }
  return 0;
}

inline int wkspace_alloc_i_checked(int** ip_ptr, unsigned long size) {
  *ip_ptr = (int*)wkspace_alloc(size);
  if (!(*ip_ptr)) {
    return 1;
  }
  return 0;
}

inline int wkspace_alloc_ui_checked(unsigned int** uip_ptr, unsigned long size) {
  *uip_ptr = (unsigned int*)wkspace_alloc(size);
  if (!(*uip_ptr)) {
    return 1;
  }
  return 0;
}

inline int wkspace_alloc_ul_checked(unsigned long** ulp_ptr, unsigned long size) {
  *ulp_ptr = (unsigned long*)wkspace_alloc(size);
  if (!(*ulp_ptr)) {
    return 1;
  }
  return 0;
}

inline int wkspace_alloc_d_checked(double** dp_ptr, unsigned long size) {
  *dp_ptr = (double*)wkspace_alloc(size);
  if (!(*dp_ptr)) {
    return 1;
  }
  return 0;
}


void wkspace_reset(void* new_base) {
  unsigned long freed_bytes = wkspace_base - (unsigned char*)new_base;
  wkspace_base = (unsigned char*)new_base;
  wkspace_left += freed_bytes;
}

char** subst_argv = NULL;

int dispmsg(int retval) {
  switch(retval) {
  case RET_NOMEM:
    printf("Error: Out of memory.  Try the --memory and/or --parallel flags.\n");
    break;
  case RET_WRITE_FAIL:
    printf("\nError: File write failure.\n");
    break;
  case RET_READ_FAIL:
    printf("\nError: File read failure.\n");
    break;
  case RET_MOREHELP:
    printf(
"\nwdist [flags...]\n\n"
"In the command line flag definitions that follow,\n"
"  * [square brackets] denote a required parameter, where the text between the\n"
"    brackets describes its nature.\n"
"  * <angle brackets> denote an optional modifier (or if '|' is present, a set\n"
"    of mutually exclusive optional modifiers).  Use the EXACT text in the\n"
"    definition, e.g. '--distance square0'.\n"
"  * {curly braces} denote an optional parameter, where the text between the\n"
"    braces describes its nature.\n\n"
"Each run must invoke at least one of the following calculations:\n\n"
"  --distance <square | square0 | triangle> <gz | bin> <ibs> <1-ibs> <snps>\n"
"    Writes a lower-triangular tab-delimited table of (weighted) genomic\n"
"    distances in SNP units to {output prefix}.dist, and a list of the\n"
"    corresponding family/individual IDs to {output prefix}.dist.id.  The first\n"
"    row of the .dist file contains a single {genotype 1-genotype 2} distance,\n"
"    the second row has the {genotype 1-genotype 3} and {genotype 2-genotype 3}\n"
"    distances in that order, etc.\n"
"    * If the 'square' or 'square0' modifier is present, a square matrix is\n"
"      written instead; 'square0' fills the upper right triangle with zeroes.\n"
"    * If the 'gz' modifier is present, a compressed .dist.gz file is written\n"
"      instead of a plain text file.\n"
"    * If the 'bin' modifier is present, a binary (square) matrix of\n"
"      double-precision floating point values, suitable for loading from R, is\n"
"      instead written to {output prefix}.dist.bin.  This can be combined with\n"
"      'square0' if you still want the upper right zeroed out, or 'triangle' if\n"
"      you don't want to pad the upper right at all.\n"
"    * If the 'ibs' modifier is present, an identity-by-state matrix is written\n"
"      to {output prefix}.mibs.  '1-ibs' causes distances expressed as genomic\n"
"      proportions (i.e. 1 - IBS) to be written to {output prefix}.mdist.\n"
"      Combine with 'snps' if you want to generate the usual .dist file as well.\n"
"  --matrix\n"
"  --distance-matrix\n"
"    These generate space-delimited text matrices, and are included for\n"
"    backwards compatibility with old scripts relying on the corresponding PLINK\n"
"    flags.  New scripts should migrate to '--distance ibs' and '--distance\n"
"    1-ibs', which support output formats better suited to parallel computation\n"
"    and have more accurate handling of missing markers.\n\n"
"  --genome <gz> <full> <unbounded>\n"
"    Identity-by-descent analysis.  This yields the same output as PLINK's\n"
"    --genome or --Z-genome, and the 'full' and 'unbounded' modifiers have the\n"
"    same effect as PLINK's --full-genome and --unbounded flags.\n\n"
"  --ibc\n"
"    Calculates inbreeding coefficients in three different ways.\n"
"    * For more details, see Yang J, Lee SH, Goddard ME and Visscher PM.  GCTA:\n"
"    a tool for Genome-wide Complex Trait Analysis.  Am J Hum Genet. 2011 Jan\n"
"    88(1): 76-82.  This paper also describes the relationship matrix\n"
"    computation we implement.\n\n"
"  --make-rel <square | square0 | triangle> <gz | bin>\n"
"             <cov | ibc1 | ibc2>\n"
"    Writes a lower-triangular variance-standardized relationship (coancestry)\n"
"    matrix to {output prefix}.rel, and corresponding IDs to\n"
"    {output prefix}.rel.id.\n"
"    * 'square', 'square0', 'triangle', 'gz', and 'bin' act as they do on\n"
"      --distance.\n"
"    * The 'cov' modifier removes the variance standardization step, causing a\n"
"      covariance matrix to be calculated instead.\n"
"    * By default, the diagonal elements in the relationship matrix are based on\n"
"      --ibc's Fhat3; use the 'ibc1' or 'ibc2' modifiers to base them on Fhat1\n"
"      or Fhat2 instead.\n"
"  --make-grm <no-gz> <cov | ibc1 | ibc2>\n"
"    Writes the relationships in a gzipped list format, describing one pair per\n"
"    line.  Note that this file explicitly stores the number of valid\n"
"    observations (where neither individual has a missing call) for each pair,\n"
"    which is useful input for some scripts.\n\n"
"  --groupdist {iters} {d}\n"
"    Considers three subsets of the distance matrix: pairs of affected\n"
"    individuals, affected-unaffected pairs, and pairs of unaffected\n"
"    individuals.  Each of these subsets has an average pairwise genomic\n"
"    distance; --groupdist computes the differences between those three\n"
"    averages, and estimates standard errors via delete-d jackknife.  Binary\n"
"    phenotype data is required.\n"
"    * If only one parameter is provided, d defaults to {number of people}^0.6\n"
"      rounded down.  With no parameters, 100k iterations are run.\n\n"
"  --regress-distance {iters} {d}\n"
"    Linear regression of pairwise genomic distances on pairwise average\n"
"    phenotypes and vice versa, using delete-d jackknife for standard errors.\n"
"    Scalar phenotype data is required.  Defaults for iters and d are the same\n"
"    as for --groupdist.\n\n"
"  --regress-rel {iters} {d}\n"
"    Linear regression of pairwise genomic relationships on pairwise average\n"
"    phenotypes, and vice versa.\n\n"
#ifndef NOLAPACK
"  --unrelated-heritability <strict> {tol} {initial covg} {initial covr}\n"
"    REML estimate of additive heritability, iterating with an accelerated\n"
"    variant of the EM algorithm until the rate of change of the log likelihood\n"
"    function is less than tol.  Scalar phenotype data is required.\n"
"    * The 'strict' modifier forces regular EM to be used.  tol defaults to\n"
"      10^{-7}, genomic covariance prior defaults to 0.45, and residual\n"
"      covariance prior defaults to (1 - covg).\n"
"    * For more details, see Vattikuti S, Guo J, Chow CC (2012) Heritability and\n"
"      Genetic Correlations Explained by Common SNPs for Metabolic Syndrome\n"
"      Traits.  PLoS Genet 8(3): e1002637.  doi:10.1371/journal.pgen.1002637\n\n"
#endif
"  --freq\n"
"    Generates an allele frequency report (identical to PLINK --freq).\n\n"
"  --write-snplist\n"
"    Write a .snplist file listing the names of all SNPs that pass the filters\n"
"    you've specified.\n\n"
"The following other flags are supported.  (Order of operations conforms to the\n"
"PLINK specification at http://pngu.mgh.harvard.edu/~purcell/plink/flow.shtml.\n"
"Or, to be more precise, it does not yet conform, but I am currently fixing\n"
"that.)\n"
"  --script [fname] : Include command-line options from file.\n"
"  --file [prefix]  : Specify prefix for .ped and .map files.  (When this flag\n"
"                     isn't present, the prefix is assumed to be 'wdist'.)\n"
"  --ped [filename] : Specify full name of .ped file.\n"
"  --map [filename] : Specify full name of .map file.\n"
"  --no-fid         : .ped file does not contain column 1 (family ID).\n"
"  --no-parents     : .ped file does not contain columns 3-4 (parents).\n"
"  --no-sex         : .ped file does not contain column 5 (sex).\n"
"  --no-pheno       : .ped file does not contain column 6 (phenotype).\n"
"  --liability      : .ped file does contain liability (column 7).\n"
"  --bfile {prefix} : Specify .bed/.bim/.fam prefix (default 'wdist').\n"
"  --bed [filename] : Specify full name of .bed file.\n"
"  --bim [filename] : Specify full name of .bim file.\n"
"  --fam [filename] : Specify full name of .fam file.\n"
"  --load-dists [f] : Load a binary TRIANGULAR distance matrix for --groupdist\n"
"                     or --regress-distance analysis, instead of recalculating\n"
"                     it from scratch.\n"
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
"                     e.g. '--chr 1 3 xy'.  The X, Y, and MT chromosome values\n"
"                     are currently unsupported, but this will change in the\n"
"                     future.\n"
"  --maf {val}      : Minor allele frequency minimum threshold (default 0.01).\n"
"                     Note that the default threshold is only applied if --maf\n"
"                     is used without an accompanying value; if you do not\n"
"                     invoke --maf, WDIST does not perform any MAF-based SNP\n"
"                     filtering at all.  Other QC flags work the same way.\n"
"  --geno {val}     : Maximum per-SNP missing (default 0.1).\n"
"  --mind {val}     : Maximum per-person missing (default 0.1).\n"
"  --hwe {val}      : Minimum Hardy-Weinberg disequilibrium p-value (exact),\n"
"                     default 0.001.\n"
"  --nonfounders    : Include all individuals in MAF/HWE calculations.\n"
"  --rel-cutoff {v} : Exclude individuals until no remaining pairs have\n"
"  --grm-cutoff {v}   relatedness greater than the given cutoff value (default\n"
"                     0.025).  Note that maximizing the remaining sample size is\n"
"                     equivalent to the NP-hard maximum independent set problem,\n"
"                     so we use a greedy algorithm instead of guaranteeing\n"
"                     optimality.  (Use the --make-rel and --keep/--remove flags\n"
"                     if you want to try to do better.)\n"
"  --ppc-gap [val]  : Minimum number of base pairs, in thousands, between\n"
"                     informative pairs of SNPs used in --genome PPC test.  500\n"
"                     if unspecified.\n"
"  --rseed [val]    : Set random number seed (relevant for jackknife standard\n"
"                     error estimation).\n"
"  --memory [val]   : Size, in MB, of initial malloc attempt.\n"
"  --threads [val]  : Maximum number of concurrent threads.\n"
"  --extract [file] : Only include SNPs in the given list.\n"
"  --exclude [file] : Exclude all SNPs in the given list.\n"
"  --keep [fname]   : Only include individuals in the given list.\n"
"  --remove [fname] : Exclude all individuals in the given list.\n"
"  --exponent [val] : When computing genomic distances, each marker has a weight\n"
"                     of (2q(1-q))^{-val}, where q is the observed MAF.  (Use\n"
"                     --read-freq if you want to explicitly specify some or all\n"
"                     of the MAFs.)\n"
"  --read-freq [filename]    : Loads MAFs from the given PLINK-style frequency\n"
"  --update-freq [filename]    file, instead of just setting them to frequencies\n"
"                              observed in the .ped/.bed file.  (This can be\n"
"                              important when performing multiple rounds of\n"
"                              filtering, because MAFs affect the behavior of\n"
"                              some filters, and observed MAFs change when you\n"
"                              filter people out.)\n"
"  --parallel [k] [n]        : Divide the output matrix into n pieces, and only\n"
"                              compute the kth piece.  The primary output file\n"
"                              will have the piece number included in its name,\n"
"                              e.g. wdist.rel.13 or wdist.rel.13.gz if k is 13.\n"
"                              Concatenating these files in order will yield the\n"
"                              full matrix of interest.  (Yes, this can be done\n"
"                              before unzipping.)\n"
"                              N.B. This cannot be used to directly write a\n"
"                              symmetric square matrix.  Choose the square0 or\n"
"                              triangle format instead, and postprocess as\n"
"                              necessary.\n"
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
"One final note.  If you provide input in a different format, WDIST will always\n"
"convert to PLINK SNP-major binary files before proceeding with its calculations\n"
"(i.e. --make-bed is effectively always on).  These files will be saved to\n"
"{output prefix}.bed, .bim, and .fam, and you are encouraged to use them\n"
"directly in future analyses.\n"
         );
    break;
  }
  free_cond(subst_argv);
  return retval;
}

// (copied from PLINK helper.cpp, modified to perform threshold test and avoid
// repeated allocation of the same memory.  Substantial further optimization
// should be possible, but it's unimportant for now.)
//
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

  // start at midpoint
  int mid = (int)((rare_copies_ll * (2 * genotypes_ll - rare_copies_ll)) / (2 * genotypes_ll));
  
  // check to ensure that midpoint and rare alleles have same parity
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
    // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
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
  // p-value calculation for p_hwe
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

int indiv_major_to_snp_major(char* indiv_major_fname, char* outname, FILE** outfile_ptr, int unfiltered_marker_ct) {
  int in_fd = open(indiv_major_fname, O_RDONLY);
  unsigned char* in_contents = MAP_FAILED;
  int unfiltered_marker_ct4 = (unfiltered_marker_ct + 3) / 4;
  struct stat sb;
  unsigned char* icoff;
  int retval;
  int indiv_ct;
  int indiv_ct4;
  int indiv_ct4l;
  long max_4blocks_in_mem;
  int superblock_offset;
  int block_last_marker;
  int ii;
  int add_val;
  int rshift_val;
  int jj;
  unsigned char* write_ptr;
  *outfile_ptr = NULL;
  if (in_fd == -1) {
    printf(errstr_fopen, indiv_major_fname);
    return RET_OPENFAIL;
  }
  // obtain file size, see example in OS X mmap() documentation
  if (fstat(in_fd, &sb) == -1) {
    goto indiv_major_to_snp_major_ret_READ_FAIL;
  }
  in_contents = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, in_fd, 0);
  if (in_contents == MAP_FAILED) {
    goto indiv_major_to_snp_major_ret_READ_FAIL;
  }
  if (fopen_checked(outfile_ptr, outname, "wb")) {
    goto indiv_major_to_snp_major_ret_OPENFAIL;
  }
  if ((in_contents[0] == 'l') && (in_contents[1] == '\x1b') && (!in_contents[2])) {
    ii = 3;
  } else if ((!in_contents[0]) && (!((sb.st_size - 1) % unfiltered_marker_ct4))) {
    ii = 1;
  } else {
    ii = 0;
  }
  icoff = &(in_contents[ii]);
  if ((sb.st_size - ii) % unfiltered_marker_ct4) {
    goto indiv_major_to_snp_major_ret_INVALID_FORMAT;
  }
  if (fwrite_checked("l\x1b\x01", 3, *outfile_ptr)) {
    goto indiv_major_to_snp_major_ret_WRITE_FAIL;
  }
  indiv_ct = sb.st_size / unfiltered_marker_ct4;
  indiv_ct4l = indiv_ct / 4;
  indiv_ct4 = (indiv_ct + 3) / 4;
  // 4 * indiv_ct4 bytes needed per 4-marker block
  max_4blocks_in_mem = wkspace_left / (4 * indiv_ct4);
  superblock_offset = 0;
  while (superblock_offset < unfiltered_marker_ct4) {
    block_last_marker = unfiltered_marker_ct - (superblock_offset * 4);
    if (block_last_marker > (max_4blocks_in_mem * 4)) {
      block_last_marker = max_4blocks_in_mem * 4;
    }
    write_ptr = wkspace_base;
    for (ii = 0; ii < block_last_marker; ii++) {
      rshift_val = (ii % 4) * 2;
      add_val = ii / 4;
      for (jj = 0; jj < indiv_ct4l; jj++) {
        *write_ptr++ = ((icoff[4 * jj * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) + (((icoff[(4 * jj + 1) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 2) + (((icoff[(4 * jj + 2) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 4) + (((icoff[(4 * jj + 3) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << 6);
      }
      if (indiv_ct % 4) {
	*write_ptr = 0;
	for (jj = 0; jj < (indiv_ct % 4); jj++) {
	  *write_ptr |= ((icoff[(jj + indiv_ct) * unfiltered_marker_ct4 + add_val] >> rshift_val) & 3) << (jj * 2);
	}
	write_ptr++;
      }
    }
    if (fwrite_checked(wkspace_base, ((long long)block_last_marker) * indiv_ct4, *outfile_ptr)) {
      goto indiv_major_to_snp_major_ret_WRITE_FAIL;
    }
    superblock_offset += max_4blocks_in_mem;
  }
  retval = 0;
  while (0) {
  indiv_major_to_snp_major_ret_INVALID_FORMAT:
    printf("Error: %s's file size is inconsistent with the marker count.\n", indiv_major_fname);
    retval = RET_INVALID_FORMAT;
    break;
  indiv_major_to_snp_major_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  indiv_major_to_snp_major_ret_OPENFAIL:
    retval = RET_OPENFAIL;
    break;
  indiv_major_to_snp_major_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  }
  if (*outfile_ptr) {
    fclose(*outfile_ptr);
    *outfile_ptr = NULL;
  }
  if (in_contents != MAP_FAILED) {
    munmap(in_contents, sb.st_size);
  }
  close(in_fd);
  return retval;
}

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
  // main_arr = packed array of equal-length items to sort
  // arr_length = number of items
  // item_length = byte count of each main_arr item
  // comparator_deref = returns positive if first > second, 0 if equal,
  //                    negative if first < second
  // secondary_arr = packed array of fixed-length records associated with the
  //                 main_arr items, to be resorted in the same way.  (e.g.
  //                 if one is building an index, this could start as a sorted
  //                 0..(n-1) sequence of integers; then, post-sort, this would
  //                 be a lookup table for the original position of each
  //                 main_arr item.)
  // secondary_item_len = byte count of each secondary_arr item
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

inline char* read_next_unsafe(char* target, char* source) {
  // assumes space- or tab-termination
  while ((*source != ' ') && (*source != '\t')) {
    *target++ = *source++;
  }
  return target;
}

inline void read_next_terminate(char* target, char* source) {
  while (!is_space_or_eoln(*source)) {
    *target++ = *source++;
  }
  *target = '\0';
}

inline char* read_next_unsafe_upd(char* target, char** source_ptr) {
  // assumes space- or tab-termination
  while ((**source_ptr != ' ') && (**source_ptr != '\t')) {
    *target++ = **source_ptr;
    *source_ptr += 1;
  }
  return target;
}

inline char* read_next_upd(char* target, char** source_ptr) {
  while (!is_space_or_eoln(**source_ptr)) {
    *target++ = **source_ptr;
    *source_ptr += 1;
  }
  return target;
}

int bsearch_fam_indiv(char* id_buf, char* lptr, int max_id_len, int filter_line_ct, char* fam_id, char* indiv_id) {
  // id_buf = workspace
  // lptr = packed, sorted list of ID strings to search over
  // fam_id and indiv_id are considered terminated by any space/eoln character
  int ii;
  int jj;
  if (!filter_line_ct) {
    return -1;
  }
  ii = strlen_se(fam_id);
  jj = strlen_se(indiv_id);
  if (ii + jj + 2 > max_id_len) {
    return -1;
  }
  memcpy(id_buf, fam_id, ii);
  id_buf[ii] = '\t';
  memcpy(&(id_buf[ii + 1]), indiv_id, jj);
  id_buf[ii + jj + 1] = '\0';
  return bsearch_str(id_buf, lptr, max_id_len, 0, filter_line_ct - 1);
}

inline int is_contained(char* id_buf, char* lptr, int max_id_len, int filter_line_ct, char* fam_id, char* indiv_id) {
  return (bsearch_fam_indiv(id_buf, lptr, max_id_len, filter_line_ct, fam_id, indiv_id) != -1);
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

inline int no_more_items(char* sptr) {
  return ((!sptr) || (*sptr == '\n') || (*sptr == '\0'));
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

int determine_max_id_len(FILE* filterfile, char* filterval, int mfilter_col, int* filter_line_ct_ptr) {
  int cur_max = 4;
  char* bufptr;
  int ii;
  int jj;

  while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
    if (*tbuf == '\n') {
      continue;
    }
    bufptr = tbuf;
    ii = 2 + strlen_se(tbuf);
    bufptr = next_item(tbuf);
    if (no_more_items(bufptr)) {
      printf(errstr_filter_format);
      return -1;
    }
    ii += strlen_se(bufptr);
    if (ii > cur_max) {
      cur_max = ii;
    }
    if (filterval) {
      for (jj = 0; jj < mfilter_col; jj++) {
	bufptr = next_item(bufptr);
      }
      if (no_more_items(bufptr)) {
        printf(errstr_filter_format);
	return -1;
      }
      if (!strncmp(filterval, bufptr, jj)) {
	*filter_line_ct_ptr += 1;
      }
    } else {
      *filter_line_ct_ptr += 1;
    }
    if (!tbuf[MAXLINELEN - 1]) {
      printf("Error: Excessively long line in filter file (max %d chars).\n", MAXLINELEN - 3);
      return -1;
    }
  }
  return cur_max;
}

char* resize_id_buf(char* id_buf, int max_id_len, int max_pid_len) {
  if (max_pid_len) {
    if (max_id_len > max_pid_len) {
      free(id_buf);
    } else {
      return id_buf;
    }
  }
  return (char*)malloc(max_id_len * sizeof(char));
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
    *larr++ = ~0LU;
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
#else
  int oo;
#endif
  if (((unsigned long)wtarr) & 15) {
    // force 16-byte alignment; can't do this at compile-time since stack
    // pointer has no 16-byte align guarantee
    wtarr++;
  }
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if (((mafs[ii] != 0.0) && (mafs[ii] < (1.0 - EPSILON))) || (!var_std)) {
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
      if (mafs[ii] == 0.0) {
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
            for (oo = 0; oo < 8; oo++) {
              *weights++ = twt4 + wtptr[oo];
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
  if ((!lhi) && ((!lli) || (!hhi))) {
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
unsigned int* missing_dbl_excluded = NULL;
int* idists;
double* dists = NULL;
char* pheno_c = NULL;
double* pheno_d = NULL;
unsigned char* ped_geno = NULL;
unsigned long* geno = NULL;
unsigned long* glptr;
double weights[2048 * BITCT];
unsigned int* weights_i = (unsigned int*)weights;
int thread_start[MAX_THREADS_P1];
int thread_start0[MAX_THREADS_P1];
double reg_tot_xy;
double reg_tot_x;
double reg_tot_y;
double reg_tot_xx;
double reg_tot_yy;
int low_ct;
int high_ct;
int jackknife_iters;
int jackknife_d;
double calc_result[9][MAX_THREADS];
unsigned long* masks;
unsigned long* mmasks;
double* marker_weights;
unsigned int* marker_weights_i = NULL;
unsigned int* missing_tot_weights;
unsigned int* indiv_missing;
unsigned int* indiv_missing_unwt = NULL;
double* jackknife_precomp = NULL;
unsigned int* genome_main;
unsigned int* marker_pos = NULL;
unsigned long marker_window[GENOME_MULTIPLEX * 2];
double* pheno_packed;

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
    if ((mafs[ii] != 0.0) && (mafs[ii] < (1.0 - EPSILON))) {
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
          if (mafs[ii] == 0.0) {
            wtarr[ii * 8] = 0;
            wtarr[ii * 8 + 3] = INFINITY;
          } else {
            wtarr[ii * 8] = INFINITY;
            wtarr[ii * 8 + 3] = 0;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 2] = -INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      } else {
        if (mafs[ii] == 0.0) {
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 3] = INFINITY;
        } else {
          wtarr[ii * 8] = INFINITY;
          wtarr[ii * 8 + 3] = 1;
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

void exclude(unsigned char* exclude_arr, int loc, int* exclude_ct_ptr) {
  int maj = loc / 8;
  int min = 1 << (loc % 8);
  if (!(exclude_arr[maj] & min)) {
    exclude_arr[maj] |= min;
    *exclude_ct_ptr += 1;
  }
}

inline int excluded(unsigned char* exclude_arr, int loc) {
  return exclude_arr[loc / 8] & (1 << (loc % 8));
}

// note that this is unsafe if you don't know there's another included marker
// or person remaining
inline int next_non_excluded(unsigned char* exclude_arr, int loc) {
  while (excluded(exclude_arr, loc)) {
    loc++;
  }
  return loc;
}

int write_ids(char* outname, int unfiltered_indiv_ct, unsigned char* indiv_exclude, char* person_id, int max_person_id_len) {
  FILE* outfile;
  int ii;
  if (fopen_checked(&outfile, outname, "w")) {
    return RET_OPENFAIL;
  }
  for (ii = 0; ii < unfiltered_indiv_ct; ii++) {
    if (!excluded(indiv_exclude, ii) && (fprintf(outfile, "%s\n", &(person_id[ii * max_person_id_len])) < 0)) {
      return RET_WRITE_FAIL;
    }
  }
  if (fclose(outfile)) {
    return RET_WRITE_FAIL;
  }
  return 0;
}

void exclude_multi(unsigned char* exclude_arr, int* new_excl, int indiv_ct, int* exclude_ct_ptr) {
  int ii;
  int true_loc = 0;
  for (ii = 0; ii < indiv_ct; ii++) {
    true_loc = next_non_excluded(exclude_arr, true_loc);
    if (new_excl[ii] == -1) {
      exclude(exclude_arr, true_loc, exclude_ct_ptr);
    }
    true_loc++;
  }
}

int sort_marker_ids(char** sorted_ids_ptr, int** id_map_ptr, int unfiltered_marker_ct, unsigned char* marker_exclude, int marker_exclude_ct, char* marker_ids, unsigned long max_marker_id_len) {
  // note that this allocates from wkspace
  int marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  int ii;
  int jj;
  *sorted_ids_ptr = (char*)wkspace_alloc(marker_ct * max_marker_id_len);
  if (!(*sorted_ids_ptr)) {
    return RET_NOMEM;
  }
  *id_map_ptr = (int*)wkspace_alloc(marker_ct * sizeof(int));
  if (!(*id_map_ptr)) {
    return RET_NOMEM;
  }
  ii = 0;
  for (jj = 0; jj < marker_ct; jj++) {
    ii = next_non_excluded(marker_exclude, ii);
    memcpy(&((*sorted_ids_ptr)[jj * max_marker_id_len]), &(marker_ids[ii * max_marker_id_len]), max_marker_id_len);
    (*id_map_ptr)[jj] = ii++;
  }
  if (qsort_ext(*sorted_ids_ptr, marker_ct, max_marker_id_len, strcmp_deref, (char*)(*id_map_ptr), sizeof(int))) {
    return RET_NOMEM;
  }
  return 0;
}

int extract_exclude_markers(char* fname, char* sorted_ids, int sorted_ids_len, unsigned long max_marker_id_len, int* id_map, int unfiltered_marker_ct, unsigned char* marker_exclude, int* marker_exclude_ct_ptr, int do_exclude) {
  FILE* infile = NULL;
  unsigned char* wkspace_mark = wkspace_base;
  unsigned char* marker_exclude_new = NULL;
  int marker_include_ct = 0;
  int unfiltered_marker_ct8 = (unfiltered_marker_ct + 7) / 8;
  int ii;
  int jj;

  if (!do_exclude) {
    marker_exclude_new = wkspace_alloc(unfiltered_marker_ct8);
    if (!marker_exclude_new) {
      return RET_NOMEM;
    }
    memset(marker_exclude_new, 255, unfiltered_marker_ct8);
  }
  if (fopen_checked(&infile, fname, "r")) {
    return RET_OPENFAIL;
  }
  while (fgets(tbuf, MAXLINELEN, infile) != NULL) {
    read_next_terminate(tbuf, tbuf);
    ii = bsearch_str(tbuf, sorted_ids, max_marker_id_len, 0, sorted_ids_len - 1);
    if (ii != -1) {
      jj = id_map[ii];
      if (do_exclude) {
        exclude(marker_exclude, jj, marker_exclude_ct_ptr);
      } else if ((!excluded(marker_exclude, jj)) && excluded(marker_exclude_new, jj)) {
	marker_exclude_new[jj / 8] -= (1 << (jj % 8));
	marker_include_ct++;
      }
    }
  }
  if (!do_exclude) {
    memcpy(marker_exclude, marker_exclude_new, unfiltered_marker_ct8);
    *marker_exclude_ct_ptr = unfiltered_marker_ct - marker_include_ct;
    wkspace_reset(wkspace_mark);
  }
  fclose(infile);
  return 0;
}

int read_freq(char* freqname, FILE** freqfile_ptr, int unfiltered_marker_ct, unsigned char* marker_exclude, int marker_exclude_ct, char* marker_ids, unsigned long max_marker_id_len, int* marker_chroms, char* marker_alleles, double* mafs) {
  unsigned char* wkspace_mark;
  char* sorted_ids;
  int* id_map;
  int ii;
  int jj;
  char cc;
  char* bufptr;
  char* bufptr2;
  double maf;
  if (fopen_checked(freqfile_ptr, freqname, "r")) {
    return RET_OPENFAIL;
  }
  if (fgets(tbuf, MAXLINELEN, *freqfile_ptr) == NULL) {
    printf("Error: Empty --read-freq file.\n");
    return RET_INVALID_FORMAT;
  }
  wkspace_mark = wkspace_base;
  if (sort_marker_ids(&sorted_ids, &id_map, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len)) {
    return RET_NOMEM;
  }
  if (!memcmp(tbuf, " CHR  ", 6)) {
    while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL) {
      bufptr = next_item(tbuf); // just initial spaces
      jj = marker_code(bufptr);
      bufptr = next_item(bufptr); // now at beginning of SNP name
      bufptr2 = next_item(bufptr);
      if (!bufptr2) {
        goto read_freq_ret_INVALID_FORMAT;
      }
      read_next_terminate(bufptr, bufptr); // destructive read (\0 at end of item)
      ii = bsearch_str(bufptr, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        if ((jj == marker_chroms[ii]) || (!jj) || (!marker_chroms[ii])) {
          cc = *bufptr2;
          bufptr2 = next_item(bufptr2);
          bufptr = next_item(bufptr2);
          if (!bufptr) {
            goto read_freq_ret_INVALID_FORMAT;
          }
          if (sscanf(bufptr, "%lg", &maf) != 1) {
            goto read_freq_ret_INVALID_FORMAT;
          }
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
    // Also support simpler frequency files with each line in the following
    // format:
    // [marker ID]\t[reference allele]\t[frequency of reference allele]\n
    do {
      bufptr = next_item(tbuf);
      if (!bufptr) {
        goto read_freq_ret_INVALID_FORMAT;
      }
      read_next_terminate(tbuf, tbuf); // destructive read
      ii = bsearch_str(tbuf, sorted_ids, max_marker_id_len, 0, unfiltered_marker_ct - marker_exclude_ct - 1);
      if (ii != -1) {
        ii = id_map[ii];
        cc = *bufptr;
        bufptr = next_item(bufptr);
	if (!bufptr) {
          goto read_freq_ret_INVALID_FORMAT;
	}
        if (sscanf(bufptr, "%lg", &maf) != 1) {
          goto read_freq_ret_INVALID_FORMAT;
        }
        if (marker_alleles[ii * 2] == cc) {
          mafs[ii] = 1.0 - maf;
        } else if (marker_alleles[ii * 2 + 1] == cc) {
          mafs[ii] = maf;
        }
      }
    } while (fgets(tbuf, MAXLINELEN, *freqfile_ptr) != NULL);
  }
  fclose_null(freqfile_ptr);
  wkspace_reset(wkspace_mark);
  return 0;
 read_freq_ret_INVALID_FORMAT:
  printf(errstr_freq_format);
  return RET_INVALID_FORMAT;
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

void collapse_arr(char* item_arr, int fixed_item_len, unsigned char* exclude_arr, int exclude_ct) {
  // collapses array of fixed-length items
  int ii = 0;
  int jj;
  while ((ii < exclude_ct) && (!excluded(exclude_arr, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < exclude_ct) {
    if (!excluded(exclude_arr, ii)) {
      memcpy(&(item_arr[(jj++) * fixed_item_len]), &(item_arr[ii * fixed_item_len]), fixed_item_len);
    }
  }
}

void collapse_bitarr(unsigned char* bitarr, unsigned char* exclude_arr, int orig_ct) {
  int ii = 0;
  int jj;
  while ((ii < orig_ct) && (!excluded(exclude_arr, ii))) {
    ii++;
  }
  jj = ii;
  while (++ii < orig_ct) {
    if (!excluded(exclude_arr, ii)) {
      if (excluded(bitarr, ii)) {
        // set bit jj
        bitarr[jj / 8] |= (1 << (jj % 8));
      } else {
	bitarr[jj / 8] &= (~(1 << (jj % 8)));
      }
      jj++;
    }
  }
}

void collapse_copy_phenod(double *target, double* pheno_d, unsigned char* indiv_exclude, int unfiltered_indiv_ct) {
  int ii = 0;
  double* target_end = &(target[indiv_ct]);
  while (target < target_end) {
    ii = next_non_excluded(indiv_exclude, ii);
    *target++ = pheno_d[ii++];
  }
}

int collapse_marker_pos(unsigned int* marker_pos, unsigned char* marker_exclude, int unfiltered_marker_ct, int marker_exclude_ct, int marker_pos_start) {
  unsigned int remap_old_idx = 0;
  unsigned int remap_new_idx = 0;
  unsigned int* mpptr = marker_pos;
  int ii;
  unsigned int uii;
  for (ii = 0; ii < unfiltered_marker_ct; ii++) {
    if (!excluded(marker_exclude, ii)) {
      uii = marker_pos[ii];
      while (remap_old_idx < uii) {
	if (marker_pos_start == remap_old_idx) {
	  marker_pos_start = remap_new_idx;
	}
	if (!excluded(marker_exclude, remap_old_idx)) {
	  remap_new_idx++;
	  if (remap_new_idx == unfiltered_marker_ct - marker_exclude_ct) {
	    remap_new_idx = remap_new_idx + GENOME_MULTIPLEX - (remap_new_idx % GENOME_MULTIPLEX);
	  }
	}
        remap_old_idx++;
      }
      *mpptr++ = remap_new_idx;
    }
  }
  return marker_pos_start;
}

void print_pheno_stdev(double* pheno_d) {
  double reg_tot_x = 0.0;
  double reg_tot_xx = 0.0;
  double dxx;
  int ii;
  for (ii = 0; ii < indiv_ct; ii++) {
    dxx = pheno_d[ii];
    reg_tot_x += dxx;
    reg_tot_xx += dxx * dxx;
  }
  printf("Phenotype stdev: %g\n", sqrt((reg_tot_xx - reg_tot_x * reg_tot_x / indiv_ct) / (indiv_ct - 1)));
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
  pick_d(tmp_cbuf, ct, dd);
  for (ii = 0; ii < ct; ii++) {
    if (tmp_cbuf[ii]) {
      *ibuf++ = ii;
    }
  }
  *ibuf = ct;
}

static inline unsigned int popcount_long(unsigned long val) {
  // the simple version, and good enough for all non-time-critical stuff
  val = val - ((val >> 1) & FIVEMASK);
#if __LP64__
  val = (val & 0x3333333333333333LU) + ((val >> 2) & 0x3333333333333333LU);
  return (((val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fLU) * 0x0101010101010101) >> 56;
#else
  val = (val & 0x33333333) + ((val >> 2) & 0x33333333);
  return (((val + (val >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24;
#endif
}


#if __LP64__
typedef union {
  __m128i vi;
  unsigned long u8[2];
  unsigned int u4[4];
} __uni16;

// SSE2 implementations of Lauradoux/Walisch popcount, combined with xor to
// handle Hamming distance, and masking to handle missingness.
// Note that the size of the popcounted buffer is a hardcoded constant
// (specifically, (MULTIPLEX_DIST / BITCT) * 16 bytes).  The current code
// assumes (MULTIPLEX / BITCT) is a multiple of 3, and no greater than 30.
static inline unsigned int popcount_xor_1mask_multiword(__m128i** xor1p, __m128i* xor2, __m128i** maskp) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  do {
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
    // Four bits represent 0-15, so we can safely add four 0-3 partial
    // bitcounts together.
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    // Accumulator stores sixteen 0-255 counts in parallel.
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  } while (xor2 < xor2_end);
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  // can get away with this when MULTIPLEX_DIST <= 960, since the 8-bit counts
  // are guaranteed to be <= 120, thus adding two together does not overflow
  // 255.
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return (unsigned int)(acc.u8[0] + acc.u8[1]);
}

static inline unsigned int popcount_xor_2mask_multiword(__m128i** xor1p, __m128i* xor2, __m128i** mask1p, __m128i* mask2) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, half1, half2;
  __uni16 acc;
  __m128i* xor2_end = &(xor2[MULTIPLEX_2DIST / 128]);

  acc.vi = _mm_setzero_si128();
  do {
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
  } while (xor2 < xor2_end);
#if MULTIPLEX_DIST > 960
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return (unsigned int)(acc.u8[0] + acc.u8[1]);
}

static inline unsigned int vec_dot_product_biased(__m128i** vec1_ptr, __m128i* vec2) {
  // Assumes vec1_ptr and vec2_ptr point to vectors of packed 2-bit values in
  // {-1, 0, 1} which we wish to compute the dot product of.  These values are
  // assumed to be stored with a bias of 2, i.e. -1 -> 01, 0 -> 10, 1 -> 11.
  //
  // MULTIPLEX_COV pairs of values are handled per function call; if fewer
  // values are present, pad the end with 10s.  Note that no missingness mask
  // is needed, because zero times anything is zero.
  //
  // Return value is ((MULTIPLEX_COV * 4) / 3) + the true dot product.
  const __m128i mx0a = {AAAAMASK, AAAAMASK};
  const __m128i mx04 = {0x4444444444444444LU, 0x4444444444444444LU};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i count1, count2, loader1, loader2;
  __uni16 acc;
  __m128i* vec2_end = &(vec2[MULTIPLEX_2COV / 128]);
  acc.vi = _mm_setzero_si128();
  do {
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;

    // This is the key two-bit product operation: (x & y) | ((x + y) & 10).
    //
    //      01  10  11
    //      ----------
    // 01 | 11  10  01
    // 10 | 10  10  10
    // 11 | 01  10  11
    //
    // Unless I made a mistake in my inspection of the possibilities, there is
    // no other combination of four x86_64 integer operations which acts as a
    // multiplexed ternary product.  (Note that, while + is NOT an admissible
    // operation by itself because it corrupts an adjacent value, masking the
    // low-order bits out immediately afterward eliminates the corruption.)  So
    // we build our heavy-duty dot product computations around the excess-2
    // encoding this implies.
    //
    // A slight improvement may be possible using newer vector instructions
    // that enable multiplexed usage of lookup tables.  Such instructions would
    // have a greater impact on the Hamming weight calculation, though (since
    // they would eliminate the need for separate mask arrays there; I have not
    // found a way to implement any permutation of
    //
    //      00  01  10  11
    //      --------------
    // 00 | 00  00  00  00
    // 01 | 00  00  10  01
    // 10 | 00  10  00  01
    // 11 | 00  01  01  00
    //
    // that is as fast as precomputing missingness masks).
    //
    // There may be other sequences of operations which directly yield partial
    // sums of two or more dot product terms, but it seems unlikely to me that
    // any would be faster than the current implementation.
    count1 = _mm_or_si128(_mm_and_si128(loader1, loader2), _mm_and_si128(_mm_add_epi64(loader1, loader2), mx0a));

    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    count2 = _mm_or_si128(_mm_and_si128(loader1, loader2), _mm_and_si128(_mm_add_epi64(loader1, loader2), mx0a));
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    count1 = _mm_add_epi64(_mm_and_si128(count1, m2), _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
    count1 = _mm_add_epi64(count1, _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    count2 = _mm_or_si128(_mm_and_si128(loader1, loader2), _mm_and_si128(_mm_add_epi64(loader1, loader2), mx0a));
    // count1 currently has 4-bit values in the 4..12 range.  To avoid
    // overflows we must subtract from these values before adding more partial
    // sums...
    count1 = _mm_add_epi64(_mm_sub_epi64(count1, mx04), _mm_add_epi64(_mm_and_si128(count2, m2), _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
    // and now count1 now has 32 four-bit 2..14 partial sums (each equal to
    // [six dot product terms] + 8)
    acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(count1, m4), _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
  } while (vec2 < vec2_end);
#if MULTIPLEX_COV > 768
  // 20..140 + 20..140 does not fit in 8 bits
  acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
  // 16..112 + 16..112 fits in 8 bits
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
  acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 16)), m16);
  acc.vi = _mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 32));
  return (unsigned int)(acc.u8[0] + acc.u8[1]);
}
#else
static inline unsigned int popcount_xor_1mask_multiword(unsigned long** xor1p, unsigned long* xor2, unsigned long** maskp) {
  // The humble 16-bit lookup table actually beats
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
  // on my development machine by a hair.
  // However, if we take the hint from Walisch-Lauradoux and postpone the
  // multiply and right shift, this is no longer true.  Ah well.
  unsigned long* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  unsigned int bit_count = 0;
  unsigned int tmp_stor;
  unsigned long loader;
  unsigned long ulii;
  unsigned long uljj;
  do {
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & (*((*maskp)++));
    ulii += (loader >> 1) & FIVEMASK;
    uljj += loader & FIVEMASK;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    // Each 8 bit slot stores a number in 0..48.  Multiplying by 0x01010101 is
    // equivalent to the left-shifts and adds we need to sum those four 8-bit
    // numbers in the high-order slot.
    bit_count += (tmp_stor * 0x01010101) >> 24;
  } while (xor2 < xor2_end);
  return bit_count;
}

static inline unsigned int popcount_xor_2mask_multiword(unsigned long** xor1p, unsigned long* xor2, unsigned long** mask1p, unsigned long* mask2) {
  unsigned long* xor2_end = &(xor2[MULTIPLEX_DIST / 16]);
  unsigned int bit_count = 0;
  unsigned int tmp_stor;
  unsigned long loader;
  unsigned long ulii;
  unsigned long uljj;
  do {
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    loader -= ((loader >> 1) & FIVEMASK);
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    ulii += (loader & 0x33333333) + ((loader >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    ulii = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    uljj = loader - ((loader >> 1) & FIVEMASK);
    loader = ((*((*xor1p)++)) ^ *xor2++) & ((*((*mask1p)++)) & *mask2++);
    loader -= ((loader >> 1) & FIVEMASK);
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    ulii += (loader & 0x33333333) + ((loader >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

    bit_count += (tmp_stor * 0x01010101) >> 24;
  } while (xor2 < xor2_end);
  return bit_count;
}

static inline unsigned int vec_dot_product_biased(unsigned long** vec1_ptr, unsigned long* vec2) {
  unsigned long* vec2_end = &(vec2[MULTIPLEX_COV / 16]);
  unsigned int biased_dot_product = 0;
  unsigned long loader1;
  unsigned long loader2;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long tmp_stor;
  do {
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    printf("%lu %lu\n", loader2 >> 16, (loader1 + loader2) & 65535);
    ulii = (loader1 & loader2) | ((loader1 + loader2) & AAAAMASK);
    printf("%lu %lu\n", ulii >> 16, ulii & 65535);
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    uljj = (loader1 & loader2) | ((loader1 + loader2) & AAAAMASK);
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    uljj = (loader1 & loader2) | ((loader1 + loader2) & AAAAMASK);
    ulii = (ulii - 0x44444444) + (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    ulii = (loader1 & loader2) | ((loader1 + loader2) & AAAAMASK);
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    uljj = (loader1 & loader2) | ((loader1 + loader2) & AAAAMASK);
    loader1 = *((*vec1_ptr)++);
    loader2 = *vec2++;
    ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
    ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    uljj = (loader1 & loader2) | ((loader1 + loader2) & AAAAMASK);
    ulii = (ulii - 0x44444444) + (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
    tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);
    biased_dot_product += (tmp_stor * 0x01010101) >> 24;
  } while (vec2 < vec2_end);
  return biased_dot_product;
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
	*idists += popcount_xor_2mask_multiword(&glptr, glptr2, &mptr, mcptr_start);
	idists++;
      }
    } else {
      while (glptr < glptr2) {
	*idists += popcount_xor_1mask_multiword(&glptr, glptr2, &mptr);
	idists++;
      }
    }
  }
}

void* calc_idist_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  int jj = thread_start[0];
  incr_dists_i(&(idists[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (unsigned long*)ped_geno, (int)tidx);
  return NULL;
}

void incr_cov(int* idists, unsigned long* geno, int tidx) {
#if __LP64__
  __m128i* glptr;
  __m128i* glptr2;
#else
  unsigned long* glptr;
  unsigned long* glptr2;
#endif
  int ii;
  int jj;
  for (ii = thread_start0[tidx]; ii < thread_start0[tidx + 1]; ii++) {
#if __LP64__
    jj = ii * (MULTIPLEX_2COV / BITCT);
    glptr = (__m128i*)geno;
    glptr2 = (__m128i*)(&(geno[jj]));
#else
    jj = ii * (MULTIPLEX_2COV / BITCT);
    glptr = geno;
    glptr2 = &(geno[jj]);
#endif
    do {
      *idists += vec_dot_product_biased(&glptr, glptr2);
      idists++;
    } while (glptr <= glptr2);
  }
}

void* calc_cov_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start0[tidx];
  int jj = thread_start0[0];
  incr_cov(&(idists[((long long)ii * (ii + 1) - (long long)jj * (jj + 1)) / 2]), geno, (int)tidx);
  return NULL;
}

void incr_genome(unsigned int* genome_main, unsigned long* geno, int tidx) {
#if __LP64__
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LU, 0x3333333333333333LU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLU, 0x0f0f0f0f0f0f0f0fLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLU, 0x00ff00ff00ff00ffLU};
  const __m128i m16 = {0x0000ffff0000ffffLU, 0x0000ffff0000ffffLU};
  __m128i xor_buf[GENOME_MULTIPLEX / BITCT];
  __m128i* xor_buf_end = &(xor_buf[GENOME_MULTIPLEX / BITCT]);
  __m128i* maskptr;
  __m128i* maskptr_fixed;
  __m128i* maskptr_fixed_tmp;
  __m128i* xor_ptr;
  __m128i* glptr_fixed_tmp;
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i count_ibs1;
  __m128i count_ibs0;
  __m128i count2_ibs1;
  __m128i count2_ibs0;
  __uni16 acc_ibs1;
  __uni16 acc_ibs0;
  unsigned long* lptr;
  __m128i* glptr;
  __m128i* glptr_fixed;
  __m128i* glptr_end;
#else
  unsigned long* glptr;
  unsigned long* glptr_fixed;
  unsigned long* glptr_end;
  unsigned long* maskptr;
  unsigned long* maskptr_fixed;
  unsigned long* maskptr_fixed_tmp;
  unsigned long* glptr_fixed_tmp;
  unsigned long xor_buf[GENOME_MULTIPLEX2 / BITCT];
  unsigned long* xor_buf_end = &(xor_buf[GENOME_MULTIPLEX2 / BITCT]);
  unsigned long* xor_ptr;
  unsigned int bit_count_ibs1 = 0;
  unsigned int bit_count_ibs0 = 0;
  unsigned long bitfield_ibs1;
  unsigned long bitfield_ibs0;
  unsigned long loader;
  unsigned long loader2;
  unsigned long tmp_stor_ibs1;
  unsigned long tmp_stor_ibs0;
#endif
  unsigned long* glptr_back;
  unsigned long ibs_incr;
  int ii;
  int jj;
  int offset;
  unsigned long uland;
  unsigned long ulval;
  unsigned long next_ppc_marker_hybrid;
  unsigned long mask_fixed_test;
  unsigned long* marker_window_ptr;
  int lowct2 = low_ct * 2;
  int highct2 = high_ct * 2;
#if __LP64__
  glptr_end = (__m128i*)(&(geno[indiv_ct * (GENOME_MULTIPLEX2 / BITCT)]));
#else
  glptr_end = &(geno[indiv_ct * (GENOME_MULTIPLEX2 / BITCT)]);
#endif
  for (ii = thread_start[tidx]; ii < thread_start[tidx + 1]; ii++) {
    jj = ii * (GENOME_MULTIPLEX2 / BITCT);
#if __LP64__
    glptr_fixed = (__m128i*)(&(geno[jj]));
    glptr = (__m128i*)(&(geno[jj + (GENOME_MULTIPLEX2 / BITCT)]));
    lptr = &(masks[jj]);
    maskptr = (__m128i*)(&(masks[jj + (GENOME_MULTIPLEX2 / BITCT)]));
    maskptr_fixed = (__m128i*)lptr;
    mask_fixed_test = *lptr++;
    for (jj = 0; jj < GENOME_MULTIPLEX2 / BITCT - 1; jj++) {
      mask_fixed_test &= *lptr++;
    }
#else
    glptr_fixed = &(geno[jj]);
    glptr = &(geno[jj + (GENOME_MULTIPLEX2 / BITCT)]);
    maskptr_fixed = &(masks[jj]);
    maskptr = maskptr_fixed;
    mask_fixed_test = *maskptr++;
    for (jj = 0; jj < GENOME_MULTIPLEX2 / BITCT - 1; jj++) {
      mask_fixed_test &= *maskptr++;
    }
#endif
    if (~mask_fixed_test) {
      while (glptr < glptr_end) {
	xor_ptr = xor_buf;
	glptr_back = (unsigned long*)glptr;
	glptr_fixed_tmp = glptr_fixed;
	maskptr_fixed_tmp = maskptr_fixed;
#if __LP64__
	acc_ibs1.vi = _mm_setzero_si128();
	acc_ibs0.vi = _mm_setzero_si128();
	do {
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
	  count_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count_ibs0;

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
	  count2_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count2_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count2_ibs0;

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), _mm_and_si128(*maskptr++, *maskptr_fixed_tmp++));
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

          count_ibs1 = _mm_add_epi64(_mm_and_si128(count_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count_ibs1, 2), m2));
          count_ibs0 = _mm_add_epi64(_mm_and_si128(count_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count_ibs0, 2), m2));
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_add_epi64(_mm_and_si128(count2_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs1, 2), m2)));
          count_ibs0 = _mm_add_epi64(count_ibs0, _mm_add_epi64(_mm_and_si128(count2_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs0, 2), m2)));
          acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_add_epi64(_mm_and_si128(count_ibs1, m4), _mm_and_si128(_mm_srli_epi64(count_ibs1, 4), m4)));
          acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_add_epi64(_mm_and_si128(count_ibs0, m4), _mm_and_si128(_mm_srli_epi64(count_ibs0, 4), m4)));
	} while (xor_ptr < xor_buf_end);
#if GENOME_MULTIPLEX > 1920
	acc_ibs1.vi = _mm_add_epi64(_mm_and_si128(acc_ibs1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs1.vi, 8), m8));
	acc_ibs0.vi = _mm_add_epi64(_mm_and_si128(acc_ibs0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs0.vi, 8), m8));
#else
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 8)), m8);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 8)), m8);
#endif
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 16)), m16);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 16)), m16);
	acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 32));
	acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 32));
	*genome_main += (unsigned int)(acc_ibs1.u8[0] + acc_ibs1.u8[1]);
	genome_main++;
        *genome_main += (unsigned int)(acc_ibs0.u8[0] + acc_ibs0.u8[1]);
#else
        bit_count_ibs1 = 0;
	bit_count_ibs0 = 0;
	do {
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 = (loader ^ loader2) & FIVEMASK;
	  bitfield_ibs0 = (loader & loader2) & FIVEMASK;
	  *xor_ptr++ = bitfield_ibs0;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & ((*maskptr++) & (*maskptr_fixed_tmp++));
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
          bitfield_ibs1 = (bitfield_ibs1 & 0x33333333) + ((bitfield_ibs1 >> 2) & 0x33333333);
	  bitfield_ibs0 = (bitfield_ibs0 & 0x33333333) + ((bitfield_ibs0 >> 2) & 0x33333333);
	  tmp_stor_ibs1 = (bitfield_ibs1 + (bitfield_ibs1 >> 4)) & 0x0f0f0f0f;
	  tmp_stor_ibs0 = (bitfield_ibs0 + (bitfield_ibs0 >> 4)) & 0x0f0f0f0f;

          bit_count_ibs1 += (tmp_stor_ibs1 * 0x01010101) >> 24;
	  bit_count_ibs0 += (tmp_stor_ibs0 * 0x01010101) >> 24;
	} while (xor_ptr < xor_buf_end);
	*genome_main += bit_count_ibs1;
	genome_main++;
	*genome_main += bit_count_ibs0;
#endif
	genome_main++;
	next_ppc_marker_hybrid = *genome_main - lowct2;
	if (next_ppc_marker_hybrid < GENOME_MULTIPLEX2) {
	  ibs_incr = 0; // hethet low-order, ibs0 high-order
	  do {
	    offset = next_ppc_marker_hybrid / BITCT;
	    marker_window_ptr = &(marker_window[offset * BITCT]);
	    next_ppc_marker_hybrid = ~0LU << (next_ppc_marker_hybrid & (BITCT - 1));
	  incr_genome_2mask_loop:
	    uland = glptr_back[offset] & (((unsigned long*)glptr_fixed)[offset]);
	    // het is represented as 11, so
	    //   (uland & (uland << 1)) & 0xaaaaaaaaaaaaaaaa
	    // stores whether a particular marker is a hethet hit in the
	    // corresponding odd bit.
	    //
	    // homozygotes are represented as 01 and 10, so
	    //   (ulxor & (ulxor >> 1)) & 0x5555555555555555
	    // stores whether a particular marker is a hom1-hom2 hit in the
	    // corresponding even bit.  (het-missing pairs also set that bit,
	    // but the masking filters that out.)
	    //
	    // ~0LU << xx masks out the bottom xx bits.
	    ulval = (((uland & (uland << 1)) & AAAAMASK) | (((unsigned long*)xor_buf)[offset]));
	    do {
	      ulval &= next_ppc_marker_hybrid;
	      if (ulval) {
		jj = __builtin_ctzl(ulval);
		next_ppc_marker_hybrid = marker_window_ptr[jj];
		ibs_incr += (1LU << ((jj & 1) * BITCT2));
	      } else if (offset < ((GENOME_MULTIPLEX2 - BITCT) / BITCT)) {
		offset++;
		next_ppc_marker_hybrid = ~0LU;
		marker_window_ptr = &(marker_window_ptr[BITCT]);
		goto incr_genome_2mask_loop;
	      } else {
		*genome_main = highct2;
		goto incr_genome_2mask_exit;
	      }
	    } while (next_ppc_marker_hybrid & (1LU << (BITCT - 1)));
	  } while (next_ppc_marker_hybrid < GENOME_MULTIPLEX2);
	  *genome_main = next_ppc_marker_hybrid + lowct2;
        incr_genome_2mask_exit:
	  genome_main++;
          *genome_main += ibs_incr & ((~0LU) >> BITCT2);
	  genome_main++;
	  *genome_main += ibs_incr >> BITCT2;
	  genome_main++;
	} else {
	  genome_main = &(genome_main[3]);
	}
      }
    } else {
      while (glptr < glptr_end) {
	xor_ptr = xor_buf;
	glptr_back = (unsigned long*)glptr;
	glptr_fixed_tmp = glptr_fixed;
#if __LP64__
	acc_ibs1.vi = _mm_setzero_si128();
	acc_ibs0.vi = _mm_setzero_si128();
	do {
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
	  count_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count_ibs0;
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count_ibs0 = _mm_add_epi64(count_ibs0, loader3);

	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
	  count2_ibs1 = _mm_and_si128(_mm_xor_si128(loader, loader2), m1);
	  count2_ibs0 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = count2_ibs0;
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);
	  loader = _mm_and_si128(_mm_xor_si128(*glptr++, *glptr_fixed_tmp++), *maskptr++);
          loader2 = _mm_srli_epi64(loader, 1);
          count2_ibs1 = _mm_add_epi64(count2_ibs1, _mm_and_si128(_mm_xor_si128(loader, loader2), m1));
	  loader3 = _mm_and_si128(_mm_and_si128(loader, loader2), m1);
	  *xor_ptr++ = loader3;
	  count2_ibs0 = _mm_add_epi64(count2_ibs0, loader3);

          count_ibs1 = _mm_add_epi64(_mm_and_si128(count_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count_ibs1, 2), m2));
          count_ibs0 = _mm_add_epi64(_mm_and_si128(count_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count_ibs0, 2), m2));
          count_ibs1 = _mm_add_epi64(count_ibs1, _mm_add_epi64(_mm_and_si128(count2_ibs1, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs1, 2), m2)));
          count_ibs0 = _mm_add_epi64(count_ibs0, _mm_add_epi64(_mm_and_si128(count2_ibs0, m2), _mm_and_si128(_mm_srli_epi64(count2_ibs0, 2), m2)));
          acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_add_epi64(_mm_and_si128(count_ibs1, m4), _mm_and_si128(_mm_srli_epi64(count_ibs1, 4), m4)));
          acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_add_epi64(_mm_and_si128(count_ibs0, m4), _mm_and_si128(_mm_srli_epi64(count_ibs0, 4), m4)));
	} while (xor_ptr < xor_buf_end);
#if GENOME_MULTIPLEX > 1920
	acc_ibs1.vi = _mm_add_epi64(_mm_and_si128(acc_ibs1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs1.vi, 8), m8));
	acc_ibs0.vi = _mm_add_epi64(_mm_and_si128(acc_ibs0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_ibs0.vi, 8), m8));
#else
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 8)), m8);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 8)), m8);
#endif
	acc_ibs1.vi = _mm_and_si128(_mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 16)), m16);
	acc_ibs0.vi = _mm_and_si128(_mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 16)), m16);
	acc_ibs1.vi = _mm_add_epi64(acc_ibs1.vi, _mm_srli_epi64(acc_ibs1.vi, 32));
	acc_ibs0.vi = _mm_add_epi64(acc_ibs0.vi, _mm_srli_epi64(acc_ibs0.vi, 32));
	*genome_main += (unsigned int)(acc_ibs1.u8[0] + acc_ibs1.u8[1]);
	genome_main++;
        *genome_main += (unsigned int)(acc_ibs0.u8[0] + acc_ibs0.u8[1]);
#else
        bit_count_ibs1 = 0;
	bit_count_ibs0 = 0;
	do {
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 = (loader ^ loader2) & FIVEMASK;
	  bitfield_ibs0 = (loader & loader2) & FIVEMASK;
	  *xor_ptr++ = bitfield_ibs0;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
	  loader = ((*glptr++) ^ (*glptr_fixed_tmp++)) & (*maskptr++);
	  loader2 = loader >> 1;
	  bitfield_ibs1 += (loader ^ loader2) & FIVEMASK;
	  loader2 = (loader & loader2) & FIVEMASK;
	  bitfield_ibs0 += loader2;
	  *xor_ptr++ = loader2;
          bitfield_ibs1 = (bitfield_ibs1 & 0x33333333) + ((bitfield_ibs1 >> 2) & 0x33333333);
	  bitfield_ibs0 = (bitfield_ibs0 & 0x33333333) + ((bitfield_ibs0 >> 2) & 0x33333333);
	  tmp_stor_ibs1 = (bitfield_ibs1 + (bitfield_ibs1 >> 4)) & 0x0f0f0f0f;
	  tmp_stor_ibs0 = (bitfield_ibs0 + (bitfield_ibs0 >> 4)) & 0x0f0f0f0f;

          bit_count_ibs1 += (tmp_stor_ibs1 * 0x01010101) >> 24;
	  bit_count_ibs0 += (tmp_stor_ibs0 * 0x01010101) >> 24;
	} while (xor_ptr < xor_buf_end);
	*genome_main += bit_count_ibs1;
	genome_main++;
	*genome_main += bit_count_ibs0;
#endif
	genome_main++;
	next_ppc_marker_hybrid = *genome_main - lowct2;
	if (next_ppc_marker_hybrid < GENOME_MULTIPLEX2) {
	  ibs_incr = 0; // hethet low-order, ibs0 high-order
	  do {
	    offset = next_ppc_marker_hybrid / BITCT;
	    marker_window_ptr = &(marker_window[offset * BITCT]);
	    next_ppc_marker_hybrid = ~0LU << (next_ppc_marker_hybrid & (BITCT - 1));
	  incr_genome_1mask_loop:
	    uland = glptr_back[offset] & (((unsigned long*)glptr_fixed)[offset]);
	    ulval = ((uland & (uland << 1)) & AAAAMASK) | (((unsigned long*)xor_buf)[offset]);
	    do {
	      ulval &= next_ppc_marker_hybrid;
	      if (ulval) {
		jj = __builtin_ctzl(ulval);
		next_ppc_marker_hybrid = marker_window_ptr[jj];
		ibs_incr += (1LU << ((jj & 1) * BITCT2));
	      } else if (offset < ((GENOME_MULTIPLEX2 - BITCT) / BITCT)) {
		offset++;
		next_ppc_marker_hybrid = ~0LU;
		marker_window_ptr = &(marker_window_ptr[BITCT]);
		goto incr_genome_1mask_loop;
	      } else {
		*genome_main = highct2;
		goto incr_genome_1mask_exit;
	      }
	    } while (next_ppc_marker_hybrid & (1LU << (BITCT - 1)));
	  } while (next_ppc_marker_hybrid < GENOME_MULTIPLEX2);
	  *genome_main = next_ppc_marker_hybrid + lowct2;
	incr_genome_1mask_exit:
	  genome_main++;
	  *genome_main += ibs_incr & ((~0LU) >> BITCT2);
	  genome_main++;
	  *genome_main += ibs_incr >> BITCT2;
	  genome_main++;
	} else {
	  genome_main = &(genome_main[3]);
	}
      }
    }
  }
}

void* calc_genome_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  int jj = thread_start[0];
  incr_genome(&(genome_main[((long long)indiv_ct * (ii - jj) - ((long long)ii * (ii + 1) - (long long)jj * (jj + 1)) / 2) * 5]), (unsigned long*)ped_geno, (int)tidx);
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
    if (mask_fixed == ~0LU) {
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
  int jj = thread_start[0];
  incr_dists(&(dists[(ii * (ii - 1) - jj * (jj - 1)) / 2]), (unsigned long*)ped_geno, (int)tidx);
  return NULL;
}

void incr_wt_dist_missing(unsigned int* mtw, int tidx) {
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
  int jj = thread_start[0];
  incr_wt_dist_missing(&(missing_tot_weights[(ii * (ii - 1) - jj * (jj - 1)) / 2]), (int)tidx);
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
  int jj = thread_start[0];
  incr_dists_r(&(rel_dists[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (unsigned long*)ped_geno, masks, (int)tidx, weights);
  return NULL;
}

void incr_dists_rm(unsigned int* idists, int tidx, int* thread_start) {
  // count missing intersection, optimized for sparsity
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
          idists[jj] += 1;
          uljj &= uljj - 1;
        }
      }
    }
    idists = &(idists[ii]);
  }
}

void* calc_missing_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  int jj = thread_start[0];
  incr_dists_rm(&(missing_dbl_excluded[((long long)ii * (ii - 1) - (long long)jj * (jj - 1)) / 2]), (int)tidx, thread_start);
  return NULL;
}

int triangle_divide(long long cur_prod, int modif) {
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod, and neither term in the product is negative.  (Note the
  // lack of a divide by two; cur_prod should also be double its "true" value
  // as a result.)
  long long vv;
  if (cur_prod == 0) {
    if (modif < 0) {
      return -modif;
    } else {
      return 0;
    }
  }
  vv = (long long)sqrt((double)cur_prod);
  while ((vv - 1) * (vv + modif - 1) >= cur_prod) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod) {
    vv++;
  }
  return vv;
}

void parallel_bounds(int ct, int start, int parallel_idx, int parallel_tot, int* bound_start_ptr, int* bound_end_ptr) {
  int modif = 1 - start * 2;
  long long ct_tot = (long long)ct * (ct + modif);
  *bound_start_ptr = triangle_divide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = triangle_divide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

// set align to 1 for no alignment
void triangle_fill(int* target_arr, int ct, int pieces, int parallel_idx, int parallel_tot, int start, int align) {
  long long ct_tr;
  long long cur_prod;
  int modif = 1 - start * 2;
  int cur_piece = 1;
  int lbound;
  int ubound;
  int ii;
  int align_m1;
  parallel_bounds(ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[pieces] = ubound;
  cur_prod = (long long)lbound * (lbound + modif);
  ct_tr = ((long long)ubound * (ubound + modif) - cur_prod) / pieces;
  while (cur_piece < pieces) {
    cur_prod += ct_tr;
    lbound = triangle_divide(cur_prod, modif);
    ii = (lbound - start) & align_m1;
    if ((ii) && (ii != align_m1)) {
      lbound = start + ((lbound - start) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (lbound > ct) {
      lbound = ct;
    }
    target_arr[cur_piece++] = lbound;
  }
}

void incr_cov_consts(double* centered_freq_buf, unsigned long* geno, double* cov_linear_terms) {
  // Note that this subtracts, not adds, so during the final matrix write step
  // additions should be formed.
  int ii;
  int jj;
  int kk;
  int mm;
  double* centered_freq_ptr;
  double* wt_ptr = &(weights[64]);
  double twt[3];
  double partial_sum = 0.0;
  unsigned long geno_val;
  for (ii = 0; ii < BITCT / 4; ii++) {
    // 192 lookup table slots required per set of 4 markers
    centered_freq_ptr = &(centered_freq_buf[ii * 4]);
    twt[0] = centered_freq_ptr[0];
    for (jj = 1; jj < 4; jj++) {
      twt[1] = twt[0] + centered_freq_ptr[1];
      wt_ptr = &(wt_ptr[16]);
      for (kk = 1; kk < 4; kk++) {
        twt[2] = twt[1] + centered_freq_ptr[2];
	wt_ptr = &(wt_ptr[4]);
        for (mm = 1; mm < 4; mm++) {
          wt_ptr++;
	  *wt_ptr++ = twt[2] + centered_freq_ptr[3];
	  *wt_ptr++ = twt[2];
	  *wt_ptr++ = twt[2] - centered_freq_ptr[3];
          twt[2] -= centered_freq_ptr[2];
	}
	twt[1] -= centered_freq_ptr[1];
      }
      twt[0] -= centered_freq_ptr[0];
    }
  }
  for (ii = 0; ii < indiv_ct; ii++) {
    wt_ptr = weights;
    partial_sum = 0.0;
    for (jj = 0; jj < 2; jj++) {
      geno_val = *geno;
      geno = &(geno[MULTIPLEX_COV / BITCT2]);
#if __LP64__
      partial_sum += wt_ptr[1344 + (geno_val >> 56)] + wt_ptr[1152 + ((geno_val >> 48) & 255)] + wt_ptr[960 + ((geno_val >> 40) & 255)] + wt_ptr[768 + ((geno_val >> 32) & 255)] + wt_ptr[576 + ((geno_val >> 24) & 255)] + wt_ptr[384 + ((geno_val >> 16) & 255)] + wt_ptr[192 + ((geno_val >> 8) & 255)] + wt_ptr[geno_val & 255];
#else
      partial_sum += wt_ptr[576 + ((geno_val >> 24) & 255)] + wt_ptr[384 + ((geno_val >> 16) & 255)] + wt_ptr[192 + ((geno_val >> 8) & 255)] + wt_ptr[geno_val & 255];
#endif
      wt_ptr = &(wt_ptr[24 * BITCT]);
    }

    *cov_linear_terms += partial_sum;
    cov_linear_terms++;
  }
}

inline int flexwrite_checked(FILE* outfile, gzFile gz_outfile, char* contents, unsigned long len) {
  if (outfile) {
    return fwrite_checked(contents, len, outfile);
  } else {
    return gzwrite_checked(gz_outfile, contents, len);
  }
}

inline int flexclose_null(FILE** outfile_ptr, gzFile* gz_outfile_ptr) {
  int ii;
  if (*outfile_ptr) {
    return fclose_null(outfile_ptr);
  } else {
    ii = gzclose(*gz_outfile_ptr);
    *gz_outfile_ptr = NULL;
    return (ii != Z_OK);
  }
}

int calc_cov(FILE* bedfile, int bed_offset, unsigned int unfiltered_marker_ct, unsigned char* marker_exclude, int marker_ct, double* mafs, int unfiltered_indiv_ct, unsigned char* indiv_exclude, int parallel_idx, int parallel_tot, char* outname, char* outname_end, char* person_id, int max_person_id_len, int calculation_type, pthread_t* threads) {
  unsigned char* wkspace_mark = wkspace_base;
  int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  int indiv_ct4long = (indiv_ct + 4 * sizeof(long) - 1) / (4 * sizeof(long));
  unsigned char* loadbuf = NULL;
  double centered_freq_buf[MULTIPLEX_COV];
  int marker_unfiltered_offset = 0;
  int markers_processed = 0;
  int retval;
  int marker_batch_size;
  int marker_buf_offset;
  int marker_buf_idx;
  int indiv_unfiltered_idx;
  int indiv_idx;
  int ii;
  int jj;
  unsigned long new_geno;
  unsigned long new_missing;
  unsigned long ulii;
  unsigned long uljj;
  unsigned char* loadbuf_ref;
  int loadbuf_shift;
  unsigned long* geno_writer;
  unsigned long* missing_writer;
  unsigned int* missing_cts;
  unsigned int* missing_ct_ref;
  double* cov_linear_terms;
  double row_and_fixed_adj;
  double centered_freq_sq_sum = 0.0;
  double dxx;
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  int write_shape = calculation_type & CALC_RELATIONSHIP_SHAPEMASK;
  int write_start_iidx;
  int write_end_iidx;
  long long write_midx_offset;
  long long write_midx_ct;
  int row_num;
  int col_num;
  int pct;
  long long pct_thresh;
  long long pct_cur_compare;
  int write_gz = calculation_type & CALC_RELATIONSHIP_GZ;
  int fixed_missing_ct;
  int* biased_dot_prod_reader;
  unsigned int* missing_dbl_exclude_reader;
  double* cov_linear_reader;
  int slen;
  char* cptr;

  fseeko(bedfile, bed_offset, SEEK_SET);
  triangle_fill(thread_start, indiv_ct, thread_ct, parallel_idx, parallel_tot, 1, 1);
  triangle_fill(thread_start0, indiv_ct, thread_ct, parallel_idx, parallel_tot, 0, 1);
  write_start_iidx = thread_start0[0];
  write_end_iidx = thread_start0[thread_ct];
  write_midx_offset = ((long long)write_start_iidx * (write_start_iidx + 1)) / 2;
  write_midx_ct = (((long long)write_end_iidx * (write_end_iidx + 1)) / 2) - write_midx_offset;
  ulii = ((unsigned long)indiv_ct * (indiv_ct - 1)) / 2;
  uljj = ulii + indiv_ct;
  if (wkspace_alloc_i_checked(&idists, uljj * sizeof(int))) {
    goto calc_cov_ret_NOMEM;
  }
  fill_int_zero(idists, uljj);
  if (wkspace_alloc_ui_checked(&missing_dbl_excluded, ulii * sizeof(int))) {
    goto calc_cov_ret_NOMEM;
  }
  fill_int_zero((int*)missing_dbl_excluded, ulii);
  if (wkspace_alloc_d_checked(&cov_linear_terms, indiv_ct * sizeof(double))) {
    goto calc_cov_ret_NOMEM;
  }
  fill_double_zero(cov_linear_terms, indiv_ct);
  if (wkspace_alloc_ul_checked(&geno, MULTIPLEX_COV * indiv_ct4long * sizeof(long))) {
    goto calc_cov_ret_NOMEM; 
  }
  if (wkspace_alloc_ul_checked(&mmasks, indiv_ct * sizeof(long))) {
    goto calc_cov_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&missing_cts, indiv_ct * sizeof(int))) {
    goto calc_cov_ret_NOMEM;
  }
  fill_int_zero((int*)missing_cts, indiv_ct);
  if (wkspace_alloc_uc_checked(&loadbuf, MULTIPLEX_COV * unfiltered_indiv_ct4 * sizeof(char))) {
    goto calc_cov_ret_NOMEM;
  }
  do {
    marker_batch_size = marker_ct - markers_processed;
    if (marker_batch_size >= MULTIPLEX_COV) {
      marker_batch_size = MULTIPLEX_COV;
    } else {
      // 0xaa = quadruple zero in our excess-2 ternary representation
      // does not increase missing counts
      memset(&(loadbuf[marker_batch_size * unfiltered_indiv_ct4]), 0xaa, (MULTIPLEX_COV - marker_batch_size) * unfiltered_indiv_ct4);
      fill_double_zero(&(centered_freq_buf[marker_batch_size]), MULTIPLEX_COV - marker_batch_size);
      fill_long_zero((long*)mmasks, indiv_ct);
    }

    for (marker_buf_idx = 0; marker_buf_idx < marker_batch_size; marker_buf_idx++) {
      if (excluded(marker_exclude, marker_unfiltered_offset)) {
	ii = 0;
	do {
	  marker_unfiltered_offset++;
	  ii++;
	} while (excluded(marker_exclude, marker_unfiltered_offset));
	fseeko(bedfile, ii * unfiltered_indiv_ct4, SEEK_CUR);
      }
      if (fread(&(loadbuf[marker_buf_idx * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	goto calc_cov_ret_READ_FAIL;
      }
      // mafs[] is actually currently the frequency of the set allele
      dxx = 2 * mafs[marker_unfiltered_offset] - 1;
      centered_freq_sq_sum += dxx * dxx;
      centered_freq_buf[marker_buf_idx] = dxx;
      marker_unfiltered_offset++;
    }

    for (marker_buf_offset = 0; marker_buf_offset < marker_batch_size; marker_buf_offset += BITCT) {
      geno_writer = &(geno[marker_buf_offset / BITCT2]);
      missing_writer = mmasks;
      missing_ct_ref = missing_cts;
      indiv_idx = 0;
      for (indiv_unfiltered_idx = 0; indiv_idx < indiv_ct; indiv_unfiltered_idx++) {
        indiv_unfiltered_idx = next_non_excluded(indiv_exclude, indiv_unfiltered_idx);
        loadbuf_ref = &(loadbuf[indiv_unfiltered_idx / 4 + marker_buf_offset * unfiltered_indiv_ct4]);
        loadbuf_shift = (indiv_unfiltered_idx % 4) * 2;
	new_geno = 0;
	new_missing = 0;
        for (marker_buf_idx = 0; marker_buf_idx < BITCT / 2; marker_buf_idx++) {
          ulii = (loadbuf_ref[marker_buf_idx * unfiltered_indiv_ct4] >> loadbuf_shift) & 3;
	  uljj = 1 ^ (ulii >> 1);
	  // ulii & uljj == 1 iff missing
	  new_missing |= (ulii & uljj) << marker_buf_idx;
          // plink 00 -> -1, 01/10 -> 0, 11 -> +1
	  new_geno |= (ulii + uljj) << (marker_buf_idx * 2);
	}
        *geno_writer++ = new_geno;
	loadbuf_ref = &(loadbuf_ref[unfiltered_indiv_ct4 * (BITCT / 2)]);
	for (marker_buf_idx = 0; marker_buf_idx < BITCT / 2; marker_buf_idx++) {
          ulii = (loadbuf_ref[marker_buf_idx * unfiltered_indiv_ct4] >> loadbuf_shift) & 3;
	  uljj = 1 ^ (ulii >> 1);
	  new_missing |= (ulii & uljj) << (marker_buf_idx + (BITCT / 2));
	  new_geno |= (ulii + uljj) << (marker_buf_idx * 2);
	}
	*geno_writer = new_geno;
	*missing_writer++ = new_missing;
	*missing_ct_ref += popcount_long(new_missing);
	missing_ct_ref++;
	geno_writer = &(geno_writer[(MULTIPLEX_COV / BITCT2) - 1]);
	indiv_idx++;
      }

      for (ulii = 1; ulii < thread_ct; ulii++) {
	if (pthread_create(&(threads[ulii - 1]), NULL, &calc_missing_thread, (void*)ulii)) {
	  goto calc_cov_ret_THREAD_CREATE_FAIL;
	}
      }
      incr_dists_rm(missing_dbl_excluded, 0, thread_start);
      for (ii = 0; ii < thread_ct - 1; ii++) {
	pthread_join(threads[ii], NULL);
      }

      incr_cov_consts(&(centered_freq_buf[marker_buf_offset]), &(geno[marker_buf_offset / BITCT2]), cov_linear_terms);
    }

    for (ulii = 1; ulii < thread_ct; ulii++) {
      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_cov_thread, (void*)ulii)) {
	goto calc_cov_ret_THREAD_CREATE_FAIL;
      }
    }
    incr_cov(idists, geno, 0);

    for (ii = 0; ii < thread_ct - 1; ii++) {
      pthread_join(threads[ii], NULL);
    }

    markers_processed += marker_batch_size;
    printf("\r%d markers complete.", markers_processed);
    fflush(stdout);
  } while (markers_processed < marker_ct);
  printf("\r                            \r");

  // subtract bias from all terms
  centered_freq_sq_sum -= ((markers_processed - marker_batch_size + MULTIPLEX_COV) * 4) / 3;

  row_num = write_start_iidx;
  biased_dot_prod_reader = idists;
  missing_dbl_exclude_reader = missing_dbl_excluded;
  pct_cur_compare = 0;

  if (calculation_type & CALC_RELATIONSHIP_BIN) {
    if (write_shape == CALC_RELATIONSHIP_SQ0) {
      fill_double_zero((double*)loadbuf, indiv_ct - 1);
    }
    strcpy(outname_end, ".rel.bin");
    if (parallel_tot > 1) {
      sprintf(&(outname_end[8]), ".%d", parallel_idx + 1);
    }
    if (fopen_checked(&outfile, outname, "wb")) {
      goto calc_cov_ret_OPENFAIL;
    }
    for (pct = 1; pct <= 100; pct++) {
      if (write_shape == CALC_RELATIONSHIP_TRI) {
        pct_thresh = write_midx_ct * pct;
      } else {
        pct_thresh = ((long long)pct) * (write_end_iidx - write_start_iidx);
      }

      while (pct_cur_compare < pct_thresh) {
	row_and_fixed_adj = cov_linear_terms[row_num] + centered_freq_sq_sum;
	cov_linear_reader = cov_linear_terms;
	fixed_missing_ct = marker_ct - missing_cts[row_num];
	missing_ct_ref = missing_cts;
	for (col_num = 0; col_num <= row_num; col_num++) {
	  if (col_num < row_num) {
	    ii = fixed_missing_ct - (*missing_ct_ref++) + (*missing_dbl_exclude_reader++);
	  } else {
	    ii = fixed_missing_ct;
	  }
	  dxx = ((*biased_dot_prod_reader++) + row_and_fixed_adj + (*cov_linear_reader++)) / ((double)ii);
	  if (fwrite_checked(&dxx, sizeof(double), outfile)) {
	    goto calc_cov_ret_WRITE_FAIL;
	  }
	}
	if (write_shape == CALC_RELATIONSHIP_TRI) {
	  pct_cur_compare += 100 * (row_num + 1);
	} else {
	  if (write_shape == CALC_RELATIONSHIP_SQ0) {
	    if (fwrite_checked(loadbuf, (indiv_ct - row_num - 1) * sizeof(double), outfile)) {
	      goto calc_cov_ret_WRITE_FAIL;
	    }
	  } else {
	    ulii = ((unsigned long)row_num * (row_num + 3)) / 2;
            for (col_num = row_num + 1; col_num < indiv_ct; col_num++) {
	      ulii += col_num;
	      dxx = (idists[ulii] + row_and_fixed_adj + (*cov_linear_reader++)) / ((double)(fixed_missing_ct - (*missing_ct_ref++) + missing_dbl_excluded[ulii]));
              if (fwrite_checked(&dxx, sizeof(double), outfile)) {
		goto calc_cov_ret_WRITE_FAIL;
	      }
	    }
	  }
	  pct_cur_compare += 100;
	}
	row_num++;
      }
      printf("\rWriting... %d%%", pct);
      fflush(stdout);
    }
  } else if (calculation_type & CALC_RELATIONSHIP_GRM) {
    if (write_gz) {
      if (parallel_tot > 1) {
        sprintf(outname_end, ".grm.%d.gz", parallel_idx + 1);
      } else {
	strcpy(outname_end, ".grm.gz");
      }
      if (gzopen_checked(&gz_outfile, outname, "wb")) {
	goto calc_cov_ret_OPENFAIL;
      }
    } else {
      if (parallel_tot > 1) {
	sprintf(outname_end, ".grm.%d", parallel_idx + 1);
      } else {
	strcpy(outname_end, ".grm");
      }
      if (fopen_checked(&outfile, outname, "wb")) {
	goto calc_cov_ret_OPENFAIL;
      }
    }
    for (pct = 1; pct <= 100; pct++) {
      printf("\rWriting... %d%%", pct - 1);
      fflush(stdout);
      pct_thresh = write_midx_ct * pct;
      while (pct_cur_compare < pct_thresh) {
	row_and_fixed_adj = cov_linear_terms[row_num] + centered_freq_sq_sum;
	cov_linear_reader = cov_linear_terms;
	fixed_missing_ct = marker_ct - missing_cts[row_num];
	missing_ct_ref = missing_cts;
	for (col_num = 0; col_num <= row_num; col_num++) {
	  if (col_num < row_num) {
	    ii = fixed_missing_ct - (*missing_ct_ref++) + (*missing_dbl_exclude_reader++);
	  } else {
	    ii = fixed_missing_ct;
	  }
	  dxx = ((*biased_dot_prod_reader++) + row_and_fixed_adj + (*cov_linear_reader++)) / ((double)ii);
	  slen = sprintf(tbuf, "%d\t%d\t%d\t%g\n", row_num + 1, col_num + 1, ii, dxx);
	  if (flexwrite_checked(outfile, gz_outfile, tbuf, slen)) {
	    goto calc_cov_ret_WRITE_FAIL;
	  }
	}
	row_num++;
      }
    }
  } else {
    if (write_shape == CALC_RELATIONSHIP_SQ0) {
      cptr = (char*)(&ulii);
      for (ii = 0; ii < sizeof(long); ii += 2) {
	cptr[ii] = '\t';
	cptr[ii + 1] = '0';
      }
      ii = (indiv_ct * 2 + sizeof(long) - 4) / sizeof(long);
      geno_writer = geno;
      for (jj = 0; jj < ii; jj++) {
	*geno_writer++ = ulii;
      }
    }
    if (write_gz) {
      if (parallel_tot > 1) {
	sprintf(outname_end, ".rel.%d.gz", parallel_idx + 1);
      } else {
	strcpy(outname_end, ".rel.gz");
      }
      if (gzopen_checked(&gz_outfile, outname, "wb")) {
	goto calc_cov_ret_OPENFAIL;
      }
    } else {
      if (parallel_tot > 1) {
	sprintf(outname_end, ".rel.%d", parallel_idx + 1);
      } else {
	strcpy(outname_end, ".rel");
      }
      if (fopen_checked(&outfile, outname, "w")) {
	goto calc_cov_ret_OPENFAIL;
      }      
    }
    for (pct = 1; pct <= 100; pct++) {
      printf("\rWriting... %d%%", pct - 1);
      fflush(stdout);
      if (write_shape == CALC_RELATIONSHIP_TRI) {
        pct_thresh = write_midx_ct * pct;
      } else {
        pct_thresh = ((long long)pct) * (write_end_iidx - write_start_iidx);
      }
      while (pct_cur_compare < pct_thresh) {
	row_and_fixed_adj = cov_linear_terms[row_num] + centered_freq_sq_sum;
	cov_linear_reader = cov_linear_terms;
	fixed_missing_ct = marker_ct - missing_cts[row_num];
	missing_ct_ref = missing_cts;
        for (col_num = 0; col_num < row_num; col_num++) {
	  if (row_num < 3) {
	    printf("%d %g %g\n", (*biased_dot_prod_reader), row_and_fixed_adj, (*cov_linear_reader));
	  }
          slen = sprintf(tbuf, "%g\t", ((*biased_dot_prod_reader++) + row_and_fixed_adj + (*cov_linear_reader++)) / ((double)(fixed_missing_ct - (*missing_ct_ref++) + (*missing_dbl_exclude_reader++))));
	  if (flexwrite_checked(outfile, gz_outfile, tbuf, slen)) {
	    goto calc_cov_ret_WRITE_FAIL;
	  }
	}
	if (row_num < 3) {
	  printf("%d %g %g\n", (*biased_dot_prod_reader), row_and_fixed_adj, (*cov_linear_reader));
	}
	slen = sprintf(tbuf, "%g", ((*biased_dot_prod_reader++) + row_and_fixed_adj + (*cov_linear_reader++)) / ((double)fixed_missing_ct));
	if (flexwrite_checked(outfile, gz_outfile, tbuf, slen)) {
	  goto calc_cov_ret_WRITE_FAIL;
	}
	if (write_shape == CALC_RELATIONSHIP_TRI) {
	  pct_cur_compare += (row_num + 1) * 100;
	} else {
	  if (write_shape == CALC_RELATIONSHIP_SQ0) {
	    if (flexwrite_checked(outfile, gz_outfile, (char*)geno, (indiv_ct - row_num - 1) * 2)) {
	      goto calc_cov_ret_WRITE_FAIL;
	    }
	  } else {
	    ulii = ((unsigned long)row_num * (row_num + 3)) / 2;
	    for (col_num = row_num + 1; col_num < indiv_ct; col_num++) {
	      ulii += col_num;
	      slen = sprintf(tbuf, "\t%g", (idists[ulii] + row_and_fixed_adj + (*cov_linear_reader++)) / ((double)(fixed_missing_ct - (*(++missing_ct_ref)) + missing_dbl_excluded[ulii])));
	      if (flexwrite_checked(outfile, gz_outfile, tbuf, slen)) {
		goto calc_cov_ret_WRITE_FAIL;
	      }
	    }
	  }
	  pct_cur_compare += 100;
	}
	if (flexwrite_checked(outfile, gz_outfile, "\n", 1)) {
	  goto calc_cov_ret_WRITE_FAIL;
	}
        row_num++;
      }
    }
  }
  if (flexclose_null(&outfile, &gz_outfile)) {
    goto calc_cov_ret_WRITE_FAIL;
  }
  printf("\rCovariance matrix written to %s.\n", outname);
  if (!parallel_idx) {
  }
  strcpy(&(outname_end[4]), ".id");
  retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);

  while (0) {
  calc_cov_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_cov_ret_OPENFAIL:
    retval = RET_OPENFAIL;
    break;
  calc_cov_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  calc_cov_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_cov_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    // todo: join threads
    break;
  }
  gzclose_cond(gz_outfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return retval;
}

void incr_dists_rm_inv(unsigned int* idists, int tidx) {
  // inverted loops for --genome --parallel
  unsigned long* glptr;
  unsigned long ulii;
  unsigned long uljj;
  int kk = indiv_ct - 1;
  int ii;
  int jj;
  for (ii = thread_start[tidx]; ii < kk; ii++) {
    ulii = mmasks[ii];
    if (ulii) {
      glptr = &(mmasks[ii + 1]);
      // jj is deliberately biased down by 1
      for (jj = ii; jj < kk; jj++) {
        uljj = (*glptr++) & ulii;
        while (uljj) {
          idists[jj] += 1;
          uljj &= uljj - 1;
        }
      }
    }
    idists = &(idists[kk - ii - 1]);
  }
}

void* calc_genomem_thread(void* arg) {
  long tidx = (long)arg;
  int ii = thread_start[tidx];
  int jj = thread_start[0];
  // f(0) = 0
  // f(1) = ic - 2
  // f(2) = 2ic - 5
  // f(3) = 3ic - 9
  // ...
  // f(n) = nic - (n+1)(n+2)/2 + 1
  incr_dists_rm_inv(&(missing_dbl_excluded[(long long)indiv_ct * (ii - jj) - ((long long)(ii + 1) * (ii + 2) - (long long)(jj + 1) * (jj + 2)) / 2]), (int)tidx);
  return NULL;
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
      *ibuf++ = curpos;
    }
    curpos++;
  } while (ibuf < ibuf_end);
}

void* groupdist_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(ped_geno[tidx * CACHEALIGN(high_ct + low_ct + (jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(ped_geno[tidx * CACHEALIGN(high_ct + low_ct + (jackknife_d + 1) * sizeof(int)) + (jackknife_d + 1) * sizeof(int)]);
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
  for (ulii = 0; ulii < 9; ulii++) {
    calc_result[ulii][tidx] = results[ulii];
  }
  return NULL;
}

// double regress_jack(int* ibuf) {
double regress_jack(int* ibuf, double* ret2_ptr) {
  int* iptr = ibuf;
  int* jptr = &(ibuf[jackknife_d]);
  int ii;
  int jj;
  int kk;
  double* dptr;
  double* dptr2;
  double neg_tot_xy = 0.0;
  double neg_tot_x = 0.0;
  double neg_tot_y = 0.0;
  double neg_tot_xx = 0.0;
  double neg_tot_yy = 0.0;
  double dxx;
  double dxx1;
  double dyy;
  while (iptr < jptr) {
    dptr2 = &(jackknife_precomp[(*iptr++) * JACKKNIFE_VALS_DIST]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  iptr = ibuf;
  for (ii = 1; ii < jackknife_d; ii++) {
    jj = *(++iptr);
    dxx1 = pheno_d[jj];
    jptr = ibuf;
    dptr = &(dists[(jj * (jj - 1)) / 2]);
    while (jptr < iptr) {
      kk = *jptr++;
      dxx = (dxx1 + pheno_d[kk]) * 0.5;
      dyy = dptr[kk];
      neg_tot_xy -= dxx * dyy;
      neg_tot_x -= dxx;
      neg_tot_y -= dyy;
      neg_tot_xx -= dxx * dxx;
      neg_tot_yy -= dyy * dyy;
    }
  }
  dxx = reg_tot_y - neg_tot_y;
  dyy = indiv_ct - jackknife_d;
  dyy = dyy * (dyy - 1.0) * 0.5;
  *ret2_ptr = ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_x - neg_tot_x) / dyy) / ((reg_tot_yy - neg_tot_yy) - dxx * dxx / dyy);
  dxx = reg_tot_x - neg_tot_x;
  return ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_y - neg_tot_y) / dyy) / ((reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

void* regress_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(ped_geno[tidx * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(ped_geno[tidx * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int)) + (jackknife_d + 1) * sizeof(int)]);
  unsigned long long ulii;
  unsigned long long uljj = jackknife_iters / 100;
  double sum = 0.0;
  double sum_sq = 0.0;
  double sum2 = 0;
  double sum2_sq = 0.0;
  double dxx;
  double ret2;
  for (ulii = 0; ulii < jackknife_iters; ulii++) {
    pick_d_small(cbuf, ibuf, indiv_ct, jackknife_d);
    dxx = regress_jack(ibuf, &ret2);
    // dxx = regress_jack(ibuf);
    sum += dxx;
    sum_sq += dxx * dxx;
    sum2 += ret2;
    sum2_sq += ret2 * ret2;
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100) / jackknife_iters;
      printf("\r%lld%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1) * jackknife_iters) / 100;
    }
  }
  calc_result[0][tidx] = sum;
  calc_result[1][tidx] = sum_sq;
  calc_result[2][tidx] = sum2;
  calc_result[3][tidx] = sum2_sq;
  return NULL;
}

double regress_rel_jack(int* ibuf, double* ret2_ptr) {
  int* iptr = ibuf;
  int* jptr = &(ibuf[jackknife_d]);
  int ii;
  int jj;
  int kk;
  double* dptr;
  double* dptr2;
  double neg_tot_xy = 0.0;
  double neg_tot_x = 0.0;
  double neg_tot_y = 0.0;
  double neg_tot_xx = 0.0;
  double neg_tot_yy = 0.0;
  double dxx;
  double dxx1;
  double dyy;
  while (iptr < jptr) {
    dptr2 = &(jackknife_precomp[(*iptr++) * JACKKNIFE_VALS_REL]);
    neg_tot_xy += *dptr2++;
    neg_tot_x += *dptr2++;
    neg_tot_y += *dptr2++;
    neg_tot_xx += *dptr2++;
    neg_tot_yy += *dptr2++;
  }
  iptr = ibuf;
  for (ii = 1; ii < jackknife_d; ii++) {
    jj = *(++iptr);
    dxx1 = pheno_packed[jj];
    jptr = ibuf;
    dptr = &(rel_dists[(jj * (jj - 1)) / 2]);
    while (jptr < iptr) {
      kk = *jptr++;
      dxx = (dxx1 + pheno_packed[kk]) * 0.5;
      dyy = dptr[kk];
      neg_tot_xy -= dxx * dyy;
      neg_tot_x -= dxx;
      neg_tot_y -= dyy;
      neg_tot_xx -= dxx * dxx;
      neg_tot_yy -= dyy * dyy;
    }
  }
  dxx = reg_tot_y - neg_tot_y;
  dyy = indiv_ct - jackknife_d;
  dyy = dyy * (dyy - 1.0) * 0.5;
  *ret2_ptr = ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_x - neg_tot_x) / dyy) / ((reg_tot_yy - neg_tot_yy) - dxx * dxx / dyy);
  dxx = reg_tot_x - neg_tot_x;
  return ((reg_tot_xy - neg_tot_xy) - dxx * (reg_tot_y - neg_tot_y) / dyy) / ((reg_tot_xx - neg_tot_xx) - dxx * dxx / dyy);
}

void* regress_rel_jack_thread(void* arg) {
  long tidx = (long)arg;
  int* ibuf = (int*)(&(ped_geno[tidx * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int))]));
  unsigned char* cbuf = &(ped_geno[tidx * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int)) + (jackknife_d + 1) * sizeof(int)]);
  unsigned long long ulii;
  unsigned long long uljj = jackknife_iters / 100;
  double sum = 0.0;
  double sum_sq = 0.0;
  double sum2 = 0;
  double sum2_sq = 0.0;
  double dxx;
  double ret2;
  for (ulii = 0; ulii < jackknife_iters; ulii++) {
    pick_d_small(cbuf, ibuf, indiv_ct, jackknife_d);
    dxx = regress_rel_jack(ibuf, &ret2);
    sum += dxx;
    sum_sq += dxx * dxx;
    sum2 += ret2;
    sum2_sq += ret2 * ret2;
    if ((!tidx) && (ulii >= uljj)) {
      uljj = (ulii * 100) / jackknife_iters;
      printf("\r%lld%%", uljj);
      fflush(stdout);
      uljj = ((uljj + 1) * jackknife_iters) / 100;
    }
  }
  calc_result[0][tidx] = sum;
  calc_result[1][tidx] = sum_sq;
  calc_result[2][tidx] = sum2;
  calc_result[3][tidx] = sum2_sq;
  return NULL;
}


int set_default_jackknife_d(int ct) {
  int dd = (int)pow((double)ct, 0.600000000001);
  printf("Setting d=%d for jackknife.\n", dd);
  return dd;
}

int regress_rel_main(unsigned char* indiv_exclude, int unfiltered_indiv_ct, int regress_rel_iters, int regress_rel_d, pthread_t* threads) {
  double* rel_ptr;
  double* pheno_ptr;
  double* pheno_ptr2;
  double* jp_fixed_ptr;
  double* jp_moving_ptr;
  int ii;
  unsigned long ulii;
  unsigned long trimatrix_size;
  double trimatrix_size_recip;
  double half_avg_pheno;
  double dxx;
  double dyy;
  double dxxyy;
  double dxxsq;
  double dyysq;
  pheno_packed = (double*)wkspace_alloc(indiv_ct * sizeof(double));
  if (!pheno_packed) {
    return RET_NOMEM;
  }
  collapse_copy_phenod(pheno_packed, pheno_d, indiv_exclude, unfiltered_indiv_ct);
  print_pheno_stdev(pheno_packed);
  trimatrix_size = ((unsigned long)indiv_ct * (indiv_ct - 1)) / 2;
  reg_tot_xy = 0.0;
  reg_tot_x = 0.0;
  reg_tot_y = 0.0;
  reg_tot_xx = 0.0;
  reg_tot_yy = 0.0;
  rel_ptr = rel_dists;
  pheno_ptr = pheno_packed;
  jackknife_precomp = (double*)wkspace_alloc(indiv_ct * JACKKNIFE_VALS_REL * sizeof(double));
  if (!jackknife_precomp) {
    return RET_NOMEM;
  }
  fill_double_zero(jackknife_precomp, indiv_ct * JACKKNIFE_VALS_REL);
  for (ii = 1; ii < indiv_ct; ii++) {
    half_avg_pheno = *(++pheno_ptr);
    pheno_ptr2 = pheno_packed;
    jp_fixed_ptr = &(jackknife_precomp[ii * JACKKNIFE_VALS_REL]);
    jp_moving_ptr = jackknife_precomp;
    while (pheno_ptr2 < pheno_ptr) {
      dxx = (half_avg_pheno + (*pheno_ptr2++)) * 0.5;
      dyy = (*rel_ptr++);
      dxxyy = dxx * dyy;
      dxxsq = dxx * dxx;
      dyysq = dyy * dyy;
      reg_tot_xy += dxxyy;
      jp_fixed_ptr[0] += dxxyy;
      *jp_moving_ptr += dxxyy;
      jp_moving_ptr++;
      reg_tot_x += dxx;
      jp_fixed_ptr[1] += dxx;
      *jp_moving_ptr += dxx;
      jp_moving_ptr++;
      reg_tot_y += dyy;
      jp_fixed_ptr[2] += dyy;
      *jp_moving_ptr += dyy;
      jp_moving_ptr++;
      reg_tot_xx += dxxsq;
      jp_fixed_ptr[3] += dxxsq;
      *jp_moving_ptr += dxxsq;
      jp_moving_ptr++;
      reg_tot_yy += dyysq;
      jp_fixed_ptr[4] += dyysq;
      *jp_moving_ptr += dyysq;
      jp_moving_ptr++;
    }
  }
  trimatrix_size_recip = 1.0 / (double)trimatrix_size;
  printf("Regression slope (y = genomic relationship, x = avg phenotype): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y * trimatrix_size_recip) / (reg_tot_xx - reg_tot_x * reg_tot_x * trimatrix_size_recip));
  printf("                 (y = avg phenotype, x = genomic relationship): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y * trimatrix_size_recip) / (reg_tot_yy - reg_tot_y * reg_tot_y * trimatrix_size_recip));
  jackknife_iters = (regress_rel_iters + thread_ct - 1) / thread_ct;
  if (regress_rel_d) {
    jackknife_d = regress_rel_d;
  } else {
    jackknife_d = set_default_jackknife_d(indiv_ct);
  }
  ped_geno = wkspace_alloc(thread_ct * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int)));
  if (!ped_geno) {
    return RET_NOMEM;
  }
  for (ulii = 1; ulii < thread_ct; ulii++) {
    if (pthread_create(&(threads[ulii - 1]), NULL, &regress_rel_jack_thread, (void*)ulii)) {
      // todo: better error handling
      return RET_THREAD_CREATE_FAIL;
    }
  }
  ulii = 0;
  regress_rel_jack_thread((void*)ulii);
  dxx = calc_result[0][0]; // relationship on pheno
  dxxsq = calc_result[1][0];

  dyy = calc_result[2][0]; // pheno on relationship
  dyysq = calc_result[3][0];

  for (ii = 0; ii < thread_ct - 1; ii++) {
    pthread_join(threads[ii], NULL);
    dxx += calc_result[0][ii + 1];
    dxxsq += calc_result[1][ii + 1];
    dyy += calc_result[2][ii + 1];
    dyysq += calc_result[3][ii + 1];
  }
  ulii = jackknife_iters * thread_ct;
  printf("\rJackknife s.e. (y = genomic relationship): %g\n", sqrt((indiv_ct / (double)jackknife_d) * (dxxsq - dxx * dxx / (double)ulii) / ((double)ulii - 1)));
  printf("               (y = phenotype): %g\n", sqrt((indiv_ct / (double)jackknife_d) * (dyysq - dyy * dyy / (double)ulii) / ((double)ulii - 1)));
  return 0;
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

// sums[idx] = matrix[idx][1] + matrix[idx][2] + ....  Ideally, we can assume
// the matrix is symmetric and only reference the FORTRAN upper right, but for
// now we can't.
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
void reml_em_one_trait(double* wkbase, double* pheno, double* covg_ref, double* covr_ref, double tol, int strict) {
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
  double covr_cur_change = 1.0;
  double covg_last_change;
  double covr_last_change;
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
    matrix_const_mult_add(wkbase, *covg_ref, *covr_ref);
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
    covr_last_change = covr_cur_change;
    covg_cur_change = (*covg_ref) * (*covg_ref) * dlg * indiv_ct_d;
    covr_cur_change = (*covr_ref) * (*covr_ref) * dle * indiv_ct_d;
    if (strict) {
      max_jump = 1.0;
    } else {
      // acceleration factor:
      // min(half covg distance to 0 or 1, covr distance to 0 or 1, pi/4 divided
      // by last angular change, 1.0 / (1 - ratio of last two step lengths),
      // MAX_EM_ACCEL)
      dxx = atan2(covg_last_change, covr_last_change) - atan2(covg_cur_change, covr_cur_change);
      if (dxx < 0.0) {
	dxx = -dxx;
      }
      if (dxx > PI) {
	dxx = 2 * PI - dxx;
      }
      dyy = sqrt((covg_cur_change * covg_cur_change + covr_cur_change * covr_cur_change) / (covg_last_change * covg_last_change + covr_last_change * covr_last_change));
      if (covg_cur_change < 0.0) {
	max_jump = *covg_ref * (-0.5) / covg_cur_change;
      } else {
	max_jump = (1.0 - *covg_ref) * 0.5 / covg_cur_change;
      }
      if (covr_cur_change < 0.0) {
	dzz = *covr_ref * (-0.5) / covr_cur_change;
      } else {
	dzz = (1.0 - *covr_ref) * 0.5 / covr_cur_change;
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
    *covr_ref += covr_cur_change * max_jump;
    ll_change = (covg_cur_change * dlg) + (covr_cur_change * dle);
    printf("\b\b\b\b\b\b      \rcovg: %g  covr: %g  EM step log likelihood change: %g", *covg_ref, *covr_ref, ll_change);
    fflush(stdout);
  } while (ll_change > tol);
  printf("\n");
}
#endif // NOLAPACK

int load_map_or_bim(FILE** mapfile_ptr, char* mapname, int binary_files, int* map_cols_ptr, int* unfiltered_marker_ct_ptr, int* marker_exclude_ct_ptr, unsigned long* max_marker_id_len_ptr, int* plink_maxsnp_ptr, unsigned char** marker_exclude_ptr, double** mafs_ptr, int** marker_allele_cts_ptr, char** marker_alleles_ptr, char** marker_ids_ptr, int** marker_chroms_ptr, unsigned int chrom_mask, int* marker_pos_start_ptr, char extractname_0, char excludename_0, char freqname_0, int ppc_gap, int calculation_type) {
  int marker_ids_needed = (extractname_0 || excludename_0 || freqname_0 || (calculation_type & (CALC_FREQ | CALC_WRITE_SNPLIST)));
  int unfiltered_marker_ct = 0;
  unsigned long max_marker_id_len = 0;
  int plink_maxsnp = 4;
  int loaded_chrom_mask = 0;
  int last_chrom = -1;
  unsigned int* mpptr = NULL;
  unsigned int* mpptr2 = NULL;
  char* bufptr;
  char* bufptr2;
  unsigned long ulii;
  int ii;
  int jj;
  int kk;
  if (fopen_checked(mapfile_ptr, mapname, "r")) {
    return RET_OPENFAIL;
  }
  // first pass: count columns, determine raw marker count, determine maximum
  // marker ID length if necessary.
  while (fgets(tbuf, MAXLINELEN - 5, *mapfile_ptr) != NULL) {
    if (!tbuf[MAXLINELEN - 6]) {
      printf("Error: Excessively long line in .map/.bim file (max %d chars).\n", MAXLINELEN - 8);
      return RET_INVALID_FORMAT;
    }
    if (tbuf[0] > ' ') {
      if (marker_ids_needed || (!unfiltered_marker_ct)) {
	bufptr = next_item(tbuf);
	if (no_more_items(bufptr)) {
	  printf(errstr_map_format);
	  return RET_INVALID_FORMAT;
	}
	ulii = strlen_se(bufptr) + 1;
	if (ulii > max_marker_id_len) {
	  max_marker_id_len = ulii;
	}
	if (ulii > (plink_maxsnp + 1)) {
	  plink_maxsnp = ulii + 1;
	}
	if (!unfiltered_marker_ct) {
	  bufptr = next_item(next_item(bufptr));
	  if (binary_files) {
	    bufptr = next_item(next_item(bufptr));
	  }
	  if (!bufptr) {
	    printf(errstr_map_format);
	    return RET_INVALID_FORMAT;
	  }
	  if (*bufptr > ' ') {
	    *map_cols_ptr = 4;
	  }
	}
      }
      unfiltered_marker_ct += 1;
    }
  }
  if (!feof(*mapfile_ptr)) {
    return RET_READ_FAIL;
  }
  if (!unfiltered_marker_ct) {
    printf("Error: No markers in .map/.bim file.");
    return RET_INVALID_FORMAT;
  }
  *unfiltered_marker_ct_ptr = unfiltered_marker_ct;
  *max_marker_id_len_ptr = max_marker_id_len;
  *plink_maxsnp_ptr = plink_maxsnp;
  rewind(*mapfile_ptr);
  ii = (unfiltered_marker_ct + 7) / 8;

  // unfiltered_marker_ct can be very large, so use wkspace for all allocations
  // that are a multiple of it

  // permanent stack allocation #1: marker_exclude
  *marker_exclude_ptr = wkspace_alloc(ii * sizeof(char));
  if (!(*marker_exclude_ptr)) {
    return RET_NOMEM;
  }
  memset(*marker_exclude_ptr, 0, ii);
  // permanent stack allocation #2: mafs
  *mafs_ptr = (double*)wkspace_alloc(unfiltered_marker_ct * sizeof(double));
  if (!(*mafs_ptr)) {
    return RET_NOMEM;
  }
  for (ii = 0; ii < unfiltered_marker_ct; ii++) {
    (*mafs_ptr)[ii] = -1.0;
  }
  // permanent stack allocation #3, if needed: marker_pos
  if (calculation_type & CALC_GENOME) {
    marker_pos = (unsigned int*)wkspace_alloc(unfiltered_marker_ct * sizeof(int));
    if (!marker_pos) {
      return RET_NOMEM;
    }
    // mpptr loads the distance values in the .map/.bim file into the
    // marker_pos array.  mpptr2 makes a second pass through the array,
    // converting these raw positions into PPC next-marker-index values.
    mpptr = marker_pos;
    mpptr2 = mpptr;
    mpptr--;
  }
  // semi-permanent stack allocation (#4): marker_weights
  marker_weights = (double*)wkspace_alloc(unfiltered_marker_ct * sizeof(double));
  if (!marker_weights) {
    return RET_NOMEM;
  }
  if (binary_files) {
    *marker_allele_cts_ptr = (int*)wkspace_alloc(unfiltered_marker_ct * 2 * sizeof(int));
    if (!(*marker_allele_cts_ptr)) {
      return RET_NOMEM;
    }
    fill_int_zero(*marker_allele_cts_ptr, unfiltered_marker_ct * 2);
    if (freqname_0 || (calculation_type & CALC_FREQ)) {
      *marker_alleles_ptr = (char*)wkspace_alloc(unfiltered_marker_ct * 2 * sizeof(char));
      if (!(*marker_alleles_ptr)) {
	return RET_NOMEM;
      }
      memset(*marker_alleles_ptr, 0, unfiltered_marker_ct * 2);
    }
  }
  if (marker_ids_needed) {
    *marker_ids_ptr = (char*)wkspace_alloc(unfiltered_marker_ct * max_marker_id_len);
    if (!(*marker_ids_ptr)) {
      return RET_NOMEM;
    }
    if (freqname_0 || (calculation_type & CALC_FREQ)) {
      *marker_chroms_ptr = (int*)wkspace_alloc(unfiltered_marker_ct * sizeof(int));
      if (!(*marker_chroms_ptr)) {
	return RET_NOMEM;
      }
    }
  }

  // second pass: actually load stuff
  for (ii = 0; ii < unfiltered_marker_ct; ii++) {
    do {
      if (fgets(tbuf, MAXLINELEN, *mapfile_ptr) == NULL) {
        return RET_READ_FAIL;
      }
    } while (tbuf[0] <= ' ');
    jj = marker_code(tbuf);
    if (jj == -1) {
      printf("Error: Invalid chromosome index in .map/.bim file.\n");
      return RET_INVALID_FORMAT;
    }
    if (calculation_type & CALC_GENOME) {
      mpptr++;
      if (jj != last_chrom) {
        if (loaded_chrom_mask & (1 << jj)) {
	  printf("Error: .map/.bim is unsorted.  Use PLINK --make-bed to remedy this.\n");
	  return RET_INVALID_FORMAT;
	}
	loaded_chrom_mask &= 1 << jj;
	if (jj) {
	  if (*marker_pos_start_ptr == -1) {
	    *marker_pos_start_ptr = ii;
	  }
	  while (mpptr2 < mpptr) {
	    *mpptr2++ = ii;
	  }
	}
	last_chrom = jj;
      }
    }
    if (*marker_chroms_ptr) {
      (*marker_chroms_ptr)[ii] = jj;
    }
    if (!(chrom_mask & (1 << marker_code(tbuf)))) {
      exclude(*marker_exclude_ptr, ii, marker_exclude_ct_ptr);
    } else {
      bufptr = next_item(tbuf);
      if (*marker_ids_ptr) {
        if (no_more_items(bufptr)) {
	  printf(errstr_map_format);
          return RET_INVALID_FORMAT;
        }
        read_next_terminate(&((*marker_ids_ptr)[ii * max_marker_id_len]), bufptr);
      }
      bufptr = next_item(bufptr);
      if (*map_cols_ptr == 4) {
	bufptr = next_item(bufptr);
      }
      if (no_more_items(bufptr)) {
        printf(errstr_map_format);
        return RET_INVALID_FORMAT;
      }
      if (*bufptr == '-') {
	if (binary_files) {
	  printf("Error: Negative marker position in .bim file.\n");
	  return RET_INVALID_FORMAT;
	}
	exclude(*marker_exclude_ptr, ii, marker_exclude_ct_ptr);
      } else {
	if ((calculation_type & CALC_GENOME) && jj) {
	  kk = atoi(bufptr);
	  *mpptr = kk;
	  if (kk >= ppc_gap) {
	    kk -= ppc_gap;
	    while (*mpptr2 < kk) {
	      *mpptr2++ = ii;
	      if (mpptr2 >= mpptr) {
		break;
	      }
	    }
	  }
	}
        if (binary_files && (*marker_alleles_ptr)) {
	  bufptr = next_item(bufptr);
	  bufptr2 = next_item(bufptr);
	  if (no_more_items(bufptr2)) {
	    printf(errstr_map_format);
	    return RET_INVALID_FORMAT;
	  }
	  (*marker_alleles_ptr)[ii * 2] = *bufptr;
	  (*marker_alleles_ptr)[ii * 2 + 1] = *bufptr2;
	}
      }
    }
  }
  if (calculation_type & CALC_GENOME) {
    ii = unfiltered_marker_ct + GENOME_MULTIPLEX - (unfiltered_marker_ct % GENOME_MULTIPLEX);
    while (mpptr2 <= mpptr) {
      *mpptr2++ = ii;
    }
  }
  if (*marker_exclude_ct_ptr) {
    printf("%d markers loaded (after excluding %d).\n", unfiltered_marker_ct - *marker_exclude_ct_ptr, *marker_exclude_ct_ptr);
  } else {
    printf("%d markers loaded.\n", unfiltered_marker_ct);
  }
  return 0;
}

int distance_open(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, char* outname, char* outname_end, char* varsuffix, char* mode, int calculation_type, int parallel_idx, int parallel_tot) {
  if (calculation_type & CALC_DISTANCE_SNPS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".dist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".dist%s", varsuffix);
    }
    strcpy(tbuf, outname_end);
    if (fopen_checked(outfile_ptr, outname, mode)) {
      return 1;
    }
  }
  if (calculation_type & CALC_DISTANCE_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mibs%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mibs%s", varsuffix);
    }
    strcpy(&(tbuf[MAX_POST_EXT]), outname_end);
    if (fopen_checked(outfile2_ptr, outname, mode)) {
      return 1;
    }
  }
  if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mdist%s.%d", varsuffix, parallel_idx + 1);
    } else {
      sprintf(outname_end, ".mdist%s", varsuffix);
    }
    strcpy(&(tbuf[MAX_POST_EXT * 2]), outname_end);
    if (fopen_checked(outfile3_ptr, outname, mode)) {
      return 1;
    }
  }
  return 0;
}

int distance_open_gz(gzFile* gz_outfile_ptr, gzFile* gz_outfile2_ptr, gzFile* gz_outfile3_ptr, char* outname, char* outname_end, int calculation_type, int parallel_idx, int parallel_tot) {
  if (calculation_type & CALC_DISTANCE_SNPS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".dist.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".dist.gz");
    }
    strcpy(tbuf, outname_end);
    if (gzopen_checked(gz_outfile_ptr, outname, "wb")) {
      return 1;
    }
  }
  if (calculation_type & CALC_DISTANCE_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mibs.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".mibs.gz");
    }
    strcpy(&(tbuf[MAX_POST_EXT]), outname_end);
    if (gzopen_checked(gz_outfile_ptr, outname, "wb")) {
      return 1;
    }
  }
  if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".mdist.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".mdist.gz");
    }
    strcpy(&(tbuf[MAX_POST_EXT * 2]), outname_end);
    if (gzopen_checked(gz_outfile_ptr, outname, "wb")) {
      return 1;
    }
  }
  return 0;
}

void distance_print_done(int format_code, char* outname, char* outname_end) {
  if (!format_code) {
    strcpy(outname_end, tbuf);
    printf("\rDistances (in SNPs) written to %s.\n", outname);
  } else if (format_code == 1) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT]));
    printf("\rIBS matrix written to %s.\n", outname);
  } else if (format_code == 2) {
    strcpy(outname_end, &(tbuf[MAX_POST_EXT * 2]));
    printf("\rDistances (proportions) written to %s.\n", outname);
  }
}

// implementation used in PLINK stats.cpp
double normdist(double zz) {
  double sqrt2pi = 2.50662827463;
  double t0;
  double z1;
  double p0;
  t0 = 1 / (1 + 0.2316419 * fabs(zz));
  z1 = exp(-0.5 * zz * zz) / sqrt2pi;
  p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
 return zz >= 0 ? 1 - p0 : p0;
}

inline int is_founder(unsigned char* founder_info, int ii) {
  return ((!founder_info) || ((founder_info[ii / 8]) & (1 << (ii % 8))));
}

inline void set_founder_info(unsigned char* founder_info, int ii) {
  founder_info[ii / 8] |= (1 << (ii % 8));
}

int calc_genome(pthread_t* threads, FILE* pedfile, int bed_offset, int marker_ct, int unfiltered_marker_ct, unsigned char* marker_exclude, double* mafs, int marker_pos_start, int unfiltered_indiv_ct, unsigned char* indiv_exclude, char* person_id, int max_person_id_len, char* paternal_id, int max_paternal_id_len, char* maternal_id, int max_maternal_id_len, unsigned char* founder_info, int parallel_idx, int parallel_tot, char* outname, char* outname_end, int nonfounders, int calculation_type, int genome_output_gz, int genome_output_full, int genome_ibd_unbounded) {
  FILE* outfile = NULL;
  gzFile gz_outfile = NULL;
  int retval = 0;
  unsigned char* wkspace_mark = wkspace_base;
  int unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  unsigned char* loadbuf; // from file
  unsigned char* gptr;
  char* cptr;
  char* cptr2;
  char* cptr3;
  char* cptr4;
  char* cptr5 = NULL;
  char* cptr6 = NULL;
  char* cptr7 = NULL;
  char* cptr8 = NULL;
  char* tbuf_mid = &(tbuf[64]);
  char* sptr_cur;
  int ii;
  int jj;
  int kk;
  int mm;
  int nn;
  int oo;
  int pp;
  int qq;
  unsigned long* glptr2;
  unsigned long* glptr3;
  unsigned int* giptr;
  unsigned int* giptr2;
  unsigned int* giptr3;
  unsigned int uii;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long ulkk;
  int missing_ct_buf[BITCT];
  int missing_ct_all;
  double marker_recip = 0.5 / (double)marker_ct;
  // assumes GENOME_MULTIPLEX >= MULTIPLEX_COV
  double maf_buf[GENOME_MULTIPLEX];
  double e00 = 0.0;
  double e01 = 0.0;
  double e02 = 0.0;
  double e11 = 0.0;
  double e12 = 0.0;
  int ibd_prect = 0;
  double dpp;
  double dqq;
  double dpp_sq;
  double dqq_sq;
  double dxx;
  double dyy;
  double dxx_recip;
  double dyy_recip;
  double dxx1;
  double dxx2;
  double dyy1;
  double dyy2;
  long long cur_line = 0;
  long long tot_cells;
  long long tot_lines;
  // imitate PLINK behavior (see Plink::prettyPrintLengths() in helper.cpp);
  int max_person_fid_len = 4;
  int max_person_iid_len = 4;
  double num_alleles;
  double num_allelesf2;
  double num_allelesf3;
  int tstc;
  int is_founder_fixed = 0;

  triangle_fill(thread_start, indiv_ct, thread_ct, parallel_tot - parallel_idx - 1, parallel_tot, 1, 1);
  // invert order, for --genome --parallel to naturally work
  for (ii = 0; ii <= thread_ct / 2; ii++) {
    jj = thread_start[ii];
    thread_start[ii] = indiv_ct - thread_start[thread_ct - ii];
    thread_start[thread_ct - ii] = indiv_ct - jj;
  }

  if (!parallel_idx) {
    cur_line = 1;
  }
  tstc = thread_start[thread_ct];
  // f(tstc) - f(thread_start[0])
  // f(0) = 0
  // f(1) = indiv_ct - 1
  // f(2) = 2indiv_ct - 3
  // ...
  // f(n) = nindiv_ct - n(n+1)/2
  tot_cells = (long long)indiv_ct * (tstc - thread_start[0]) - ((long long)tstc * (tstc + 1) - (long long)thread_start[0] * (thread_start[0] + 1)) / 2;
  tot_lines = cur_line + tot_cells;
  for (ii = 0; ii < unfiltered_indiv_ct; ii++) {
    if (!excluded(indiv_exclude, ii)) {
      cptr = &(person_id[ii * max_person_id_len]);
      jj = strlen_se(cptr);
      cptr2 = next_item(cptr);
      cptr[jj] = '\0';
      if (jj > max_person_fid_len) {
	max_person_fid_len = jj + 2;
      }
      jj = strlen_se(cptr2);
      if (jj > max_person_iid_len) {
        max_person_iid_len = jj + 2;
      }
    }
  }
  if (wkspace_alloc_ui_checked(&missing_dbl_excluded, tot_cells * sizeof(int))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&indiv_missing_unwt, indiv_ct * sizeof(int))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ui_checked(&genome_main, tot_cells * 5 * sizeof(int))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_uc_checked(&loadbuf, GENOME_MULTIPLEX * unfiltered_indiv_ct4)) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_uc_checked(&ped_geno, indiv_ct * (GENOME_MULTIPLEX / 4))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&masks, indiv_ct * (GENOME_MULTIPLEX / 4))) {
    goto calc_genome_ret_NOMEM;
  }
  if (wkspace_alloc_ul_checked(&mmasks, indiv_ct * sizeof(long))) {
    goto calc_genome_ret_NOMEM;
  }

  fill_int_zero((int*)missing_dbl_excluded, tot_cells);
  fill_int_zero((int*)indiv_missing_unwt, indiv_ct);
  if (!marker_pos_start) {
    fill_int_zero((int*)genome_main, tot_cells * 5);
  } else {
    ii = marker_pos_start * 2;
    giptr = genome_main;
    for (ulii = 0; ulii < tot_cells; ulii++) {
      *giptr++ = 0;
      *giptr++ = 0;
      *giptr++ = ii;
      *giptr++ = 0;
      *giptr++ = 0;
    }
  }
  if (!excluded(marker_exclude, 0)) {
    fseeko(pedfile, bed_offset, SEEK_SET);
  }

  ii = 0; // raw marker index
  low_ct = 0; // after excluding missing
  do {
    kk = marker_ct - low_ct;
    if (kk > GENOME_MULTIPLEX) {
      kk = GENOME_MULTIPLEX;
    }
    glptr2 = marker_window;
    for (jj = 0; jj < kk; jj++) {
      if (excluded(marker_exclude, ii)) {
	do {
	  ii++;
	} while (excluded(marker_exclude, ii));
	fseeko(pedfile, bed_offset + ii * unfiltered_indiv_ct4, SEEK_SET);
      }
      if (fread(&(loadbuf[jj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	retval = RET_READ_FAIL;
	goto calc_genome_ret_1;
      }
      maf_buf[jj] = mafs[ii];
      ulii = 2 * (marker_pos[jj + low_ct] - low_ct);
      if (ulii < BITCT + (2 * (jj & (~(BITCT2 - 1))))) {
	ulii = ~0LU << (ulii & (BITCT - 1));
      }
      *glptr2++ = ulii;
      *glptr2++ = ulii;
      ii++;
    }
    if (kk < GENOME_MULTIPLEX) {
      memset(&(loadbuf[kk * unfiltered_indiv_ct4]), 0, (GENOME_MULTIPLEX - kk) * unfiltered_indiv_ct4);
      fill_long_zero((long*)ped_geno, indiv_ct * (GENOME_MULTIPLEX / BITCT2));
      fill_long_zero((long*)masks, indiv_ct * (GENOME_MULTIPLEX / BITCT2));
      for (mm = kk * 2; mm < GENOME_MULTIPLEX2; mm++) {
	*glptr2++ = GENOME_MULTIPLEX2;
      }
    }
    high_ct = low_ct + kk;
    for (jj = 0; jj < kk; jj += BITCT) {
      glptr = &(((unsigned long*)ped_geno)[jj / BITCT2]);
      glptr2 = &(masks[jj / BITCT2]);
      glptr3 = mmasks;
      giptr = indiv_missing_unwt;
      nn = 0; // raw indiv index
      fill_int_zero(missing_ct_buf, BITCT);
      missing_ct_all = 0;
      for (mm = 0; mm < indiv_ct; mm++) {
	while (excluded(indiv_exclude, nn)) {
	  nn++;
	}
	oo = (nn % 4) * 2;
	ulii = 0;
	ulkk = 0;
        gptr = &(loadbuf[nn / 4 + jj * unfiltered_indiv_ct4]);
	qq = (nonfounders || is_founder(founder_info, mm));
	if (qq) {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[pp] += 1;
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	} else {
	  missing_ct_all++;
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	}
	ulii ^= FIVEMASK;
	*glptr++ = ulii;
	ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	*glptr2++ = ulii * 3;
	*glptr3 = ulkk;
	ulii = 0;
	ulkk = 0;
	gptr = &(loadbuf[nn / 4 + (jj + BITCT2) * unfiltered_indiv_ct4]);
	if (qq) {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      missing_ct_buf[pp + BITCT2] += 1;
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	} else {
	  for (pp = 0; pp < BITCT2; pp++) {
	    uljj = (gptr[pp * unfiltered_indiv_ct4] >> oo) & 3;
	    ulii |= uljj << (pp * 2);
	    if (uljj == 1) {
	      ulkk |= 1LU << pp;
	      *giptr += 1;
	    }
	  }
	}
	ulii ^= FIVEMASK;
	*glptr = ulii;
	ulii = (ulii | (ulii >> 1)) & FIVEMASK;
	*glptr2 = ulii * 3;
	*glptr3++ |= ulkk << BITCT2;
	glptr = &(glptr[(GENOME_MULTIPLEX2 / BITCT) - 1]);
	glptr2 = &(glptr2[(GENOME_MULTIPLEX2 / BITCT) - 1]);
	giptr++;
	nn++;
      }
      nn = kk - jj;
      if (nn > BITCT) {
	nn = BITCT;
      }
      for (mm = 0; mm < nn; mm++) {
	dpp = maf_buf[jj + mm];
	dqq = 1.0 - dpp;
	num_alleles = (double)(2 * (indiv_ct - missing_ct_buf[mm] - missing_ct_all));
	if ((num_alleles > 3) && (dpp > 0.0) && (dqq > 0.0)) {
	  // update e00, e01, e02, e11, e12, ibd_prect
	  // see Plink::preCalcGenomeIBD() in genome.cpp
	  num_allelesf2 = num_alleles * num_alleles / ((num_alleles - 1) * (num_alleles - 2));
	  num_allelesf3 = num_allelesf2 * num_alleles / (num_alleles - 3);
	  dxx = dpp * num_alleles;
	  dyy = dqq * num_alleles;
	  dxx_recip = 1.0 / dxx;
	  dyy_recip = 1.0 / dyy;
	  dpp_sq = dpp * dpp;
	  dqq_sq = dqq * dqq;
	  dxx1 = (dxx - 1) * dxx_recip;
	  dxx2 = dxx1 * (dxx - 2) * dxx_recip;
	  dyy1 = (dyy - 1) * dyy_recip;
	  dyy2 = dyy1 * (dyy - 2) * dyy_recip;
	  e00 += 2 * dpp_sq * dqq_sq * dxx1 * dyy1 * num_allelesf3;
	  e01 += 4 * dpp * dqq * num_allelesf3 * (dpp_sq * dxx2 + dqq_sq * dyy2);
	  e02 += num_allelesf3 * (dqq_sq * dqq_sq * dyy2 * (dyy - 3) * dyy_recip + dpp_sq * dpp_sq * dxx2 * (dxx - 3) * dxx_recip + 4 * dpp_sq * dqq_sq * dxx1 * dyy1);
	  e11 += 2 * dpp * dqq * num_allelesf2 * (dpp * dxx1 + dqq * dyy1);
	  e12 += num_allelesf2 * (dpp_sq * dpp * dxx2 + dqq_sq * dqq * dyy2 + dpp_sq * dqq * dxx1 + dpp * dqq_sq * dyy1);
	  ibd_prect++;
	}
      }
      for (ulii = 1; ulii < thread_ct; ulii++) {
	if (pthread_create(&(threads[ulii - 1]), NULL, &calc_genomem_thread, (void*)ulii)) {
	  goto calc_genome_ret_THREAD_CREATE_FAIL;
	}
      }
      incr_dists_rm_inv(missing_dbl_excluded, 0);
      for (nn = 0; nn < thread_ct - 1; nn++) {
	pthread_join(threads[nn], NULL);
      }
    }
    for (ulii = 1; ulii < thread_ct; ulii++) {
      if (pthread_create(&(threads[ulii - 1]), NULL, &calc_genome_thread, (void*)ulii)) {
	goto calc_genome_ret_THREAD_CREATE_FAIL;
      }
    }
    incr_genome(genome_main, (unsigned long*)ped_geno, 0);
    for (nn = 0; nn < thread_ct - 1; nn++) {
      pthread_join(threads[nn], NULL);
    }
    low_ct = high_ct;
    printf("\r%d markers complete.", low_ct);
    fflush(stdout);
  } while (low_ct < marker_ct);
  printf("\rIBD calculations complete.  \n");
  dxx = 1.0 / (double)ibd_prect;
  e00 *= dxx;
  e01 *= dxx;
  e02 *= dxx;
  e11 *= dxx;
  e12 *= dxx;

  if (calculation_type & CALC_PLINK_IBS_MATRIX) {
    strcpy(outname_end, ".mibs");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPENFAIL;
    }
    giptr = genome_main;
    giptr2 = missing_dbl_excluded;
    kk = 1;
    for (ii = 0; ii < indiv_ct; ii++) {
      giptr3 = indiv_missing_unwt;
      uii = marker_ct - giptr3[ii];
      uljj = ii - 1;
      for (ulii = 0; ulii < ii; ulii++) {
	if (fprintf(outfile, "%g ", 1.0 - ((double)(genome_main[uljj * 5] + 2 * genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + missing_dbl_excluded[uljj])))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	uljj += indiv_ct - ulii - 2;
      }
      if (fwrite_checked("1 ", 2, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      giptr3++;
      for (jj = ii + 1; jj < indiv_ct; jj++) {
	if (fprintf(outfile, "%g ", 1.0 - ((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++))))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	giptr = &(giptr[5]);
      }
      if (fwrite_checked("\n", 1, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (ii * 100 >= (kk * indiv_ct)) {
	kk = (ii * 100) / indiv_ct;
	printf("\rWriting... %d%%", kk++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    printf("\rIBS matrix written to %s.\n", outname);
    strcpy(outname_end, ".mibs.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }

  if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
    strcpy(outname_end, ".mdist");
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPENFAIL;
    }
    giptr = genome_main;
    giptr2 = missing_dbl_excluded;
    kk = 1;
    for (ii = 0; ii < indiv_ct; ii++) {
      giptr3 = indiv_missing_unwt;
      uii = marker_ct - giptr3[ii];
      uljj = ii - 1;
      for (ulii = 0; ulii < ii; ulii++) {
	if (fprintf(outfile, "%g ", ((double)(genome_main[uljj * 5] + 2 * genome_main[uljj * 5 + 1])) / ((double)(2 * (uii - (*giptr3++) + missing_dbl_excluded[uljj])))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	uljj += indiv_ct - ulii - 2;
      }
      if (fwrite_checked("0 ", 2, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      giptr3++;
      for (jj = ii + 1; jj < indiv_ct; jj++) {
	if (fprintf(outfile, "%g ", ((double)((*giptr) + 2 * giptr[1])) / ((double)(2 * (uii - (*giptr3++) + (*giptr2++))))) < 0) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
	giptr = &(giptr[5]);
      }
      if (fwrite_checked("\n", 1, outfile)) {
	goto calc_genome_ret_WRITE_FAIL;
      }
      if (ii * 100 >= (kk * indiv_ct)) {
	kk = (ii * 100) / indiv_ct;
	printf("\rWriting... %d%%", kk++);
	fflush(stdout);
      }
    }
    if (fclose_null(&outfile)) {
      goto calc_genome_ret_WRITE_FAIL;
    }
    printf("\rDistances (proportions) written to %s.\n", outname);
    strcpy(outname_end, ".mdist.id");
    retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
    if (retval) {
      goto calc_genome_ret_1;
    }
  }

  sprintf(tbuf, " %%%ds %%%ds %%%ds %%%ds ", max_person_fid_len - 1, max_person_iid_len - 1, max_person_fid_len - 1, max_person_iid_len - 1);
  if (!parallel_idx) {
    if (genome_output_full) {
      sprintf(tbuf_mid, "%%%ds%%%ds%%%ds%%%ds RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO    IBS0    IBS1    IBS2  HOMHOM  HETHET\n", max_person_fid_len, max_person_iid_len, max_person_fid_len, max_person_iid_len);
    } else {
      sprintf(tbuf_mid, "%%%ds%%%ds%%%ds%%%ds RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO\n", max_person_fid_len, max_person_iid_len, max_person_fid_len, max_person_iid_len);
    }
  }
  mm = 1;
  ulii = 0;
  uljj = 0;
  if (genome_output_gz) {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d.gz", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome.gz");
    }
    if (gzopen_checked(&gz_outfile, outname, "wb")) {
      goto calc_genome_ret_1;
    }
    if (!parallel_idx) {
      if (!gzprintf(gz_outfile, tbuf_mid, " FID1", " IID1", " FID2", " IID2")) {
	goto calc_genome_ret_WRITE_FAIL;
      }
    }
  } else {
    if (parallel_tot > 1) {
      sprintf(outname_end, ".genome.%d", parallel_idx + 1);
    } else {
      strcpy(outname_end, ".genome");
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto calc_genome_ret_OPENFAIL;
    }
    if (!parallel_idx) {
      if (fprintf(outfile, tbuf_mid, " FID1", " IID1", " FID2", " IID2") < 0) {
	goto calc_genome_ret_WRITE_FAIL;
      }
    }
  }
  for (ii = thread_start[0]; ii < tstc; ii++) {
    cptr = &(person_id[ii * max_person_id_len]);
    cptr2 = &(cptr[strlen(cptr) + 1]);
    if (paternal_id) {
      cptr5 = &(paternal_id[ii * max_paternal_id_len]);
      cptr6 = &(maternal_id[ii * max_maternal_id_len]);
      is_founder_fixed = is_founder(founder_info, ii);
    }
    for (jj = ii + 1; jj < indiv_ct; jj++) {
      cptr3 = &(person_id[jj * max_person_id_len]);
      cptr4 = &(cptr3[strlen(cptr3) + 1]);
      if (paternal_id) {
	cptr7 = &(paternal_id[jj * max_paternal_id_len]);
	cptr8 = &(maternal_id[jj * max_maternal_id_len]);
      }
      sptr_cur = &(tbuf_mid[sprintf(tbuf_mid, tbuf, cptr, cptr2, cptr3, cptr4)]);
      if (!strcmp(cptr, cptr3)) {
	// eventually need to replace this with a recursive pedigree analysis
	while (1) {
	  if (paternal_id) {
	    if (!(is_founder_fixed || is_founder(founder_info, jj))) {
	      if ((!strcmp(cptr5, cptr7)) && (!strcmp(cptr6, cptr8))) {
		sptr_cur += sprintf(sptr_cur, "FS   0.5");
		break;
	      } else if ((!strcmp(cptr5, cptr7)) || (!strcmp(cptr6, cptr8))) {
		sptr_cur += sprintf(sptr_cur, "HS  0.25");
		break;
	      }
	    }
	    if ((!strcmp(cptr5, cptr4)) || (!strcmp(cptr6, cptr4)) || (!strcmp(cptr7, cptr2)) || (!strcmp(cptr8, cptr2))) {
	      sptr_cur += sprintf(sptr_cur, "PO   0.5");
	      break;
	    }
	  }
	  memcpy(sptr_cur, "OT     0", 8);
	  sptr_cur += 8;
	  break;
	}
      } else {
	memcpy(sptr_cur, "UN    NA", 8);
	sptr_cur += 8;
      }
      nn = marker_ct - indiv_missing_unwt[ii] - indiv_missing_unwt[jj] + missing_dbl_excluded[uljj];
      oo = nn - genome_main[ulii] - genome_main[ulii + 1];
      dxx = (double)genome_main[ulii + 1] / (e00 * nn);
      dyy = ((double)genome_main[ulii] - dxx * e01 * nn) / (e11 * nn);
      dxx1 = ((double)oo - nn * (dxx * e02 + dyy * e12)) / ((double)nn);
      if (!genome_ibd_unbounded) {
	if (dxx > 1) {
	  dxx = 1;
	  dyy = 0;
	  dxx1 = 0;
	} else if (dyy > 1) {
	  dyy = 1;
	  dxx = 0;
	  dxx1 = 0;
	} else if (dxx1 > 1) {
	  dxx1 = 1;
	  dyy = 0;
	  dxx1 = 0;
	} else if (dxx < 0) {
	  dxx2 = 1.0 / (dyy + dxx1);
	  dyy *= dxx2;
	  dxx1 *= dxx2;
	  dxx = 0;
	}
	if (dyy < 0) {
	  dxx2 = 1.0 / (dxx + dxx1);
	  dxx *= dxx2;
	  dxx1 *= dxx2;
	  dyy = 0;
	}
	if (dxx1 < 0) {
	  dxx2 = 1.0 / (dxx + dyy);
	  dxx *= dxx2;
	  dyy *= dxx2;
	  dxx1 = 0;
	}
      }

      sptr_cur += sprintf(sptr_cur, "  %1.4f  %1.4f  %1.4f  %1.4f  ", dxx, dyy, dxx1, dyy * 0.5 + dxx1);
      if (pheno_c) {
	if ((pheno_c[ii] != 1) && (pheno_c[jj] != 1)) {
	  memcpy(sptr_cur, "-1", 2);
	} else if ((pheno_c[ii] == 1) && (pheno_c[jj] == 1)) {
	  memcpy(sptr_cur, " 1", 2);
	} else {
	  memcpy(sptr_cur, " 0", 2);
	}
      } else {
	memcpy(sptr_cur, "NA", 2);
      }
      sptr_cur += 2;
      dxx = (double)genome_main[ulii + 4];
      dyy = (double)genome_main[ulii + 3];
      dxx1 = 1.0 / ((double)(genome_main[ulii + 4] + genome_main[ulii + 3]));
      dxx2 = normdist((dxx * dxx1 - 0.666666) / (sqrt(0.2222222 * dxx1)));
      if (genome_main[ulii + 3]) {
	sptr_cur += sprintf(sptr_cur, "  %1.6f  %1.4f  %1.4f", 1.0 - marker_recip * (genome_main[ulii] + 2 * genome_main[ulii + 1]), dxx2, dxx / dyy);
      } else {
	sptr_cur += sprintf(sptr_cur, "  %1.6f  %1.4f      NA", 1.0 - marker_recip * (genome_main[ulii] + 2 * genome_main[ulii + 1]), dxx2);
      }
      if (genome_output_full) {
	sptr_cur += sprintf(sptr_cur, " %7d %7d %7d  %1.4f  %1.4f\n", genome_main[ulii + 1], genome_main[ulii], oo, dyy * dxx1, dxx * dxx1);
      } else {
	*sptr_cur++ = '\n';
      }
      if (genome_output_gz) {
	if (gzwrite_checked(gz_outfile, tbuf_mid, (unsigned long)(sptr_cur - tbuf_mid))) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
      } else {
	if (fwrite_checked(tbuf_mid, (unsigned long)(sptr_cur - tbuf_mid), outfile)) {
	  goto calc_genome_ret_WRITE_FAIL;
	}
      }
      ulii += 5;
      uljj++;
    }
    cur_line += indiv_ct - ii - 1;
    if (cur_line * 100 >= tot_lines * mm) {
      mm = (cur_line * 100) / tot_lines;
      printf("\rWriting... %d%%", mm++);
      fflush(stdout);
    }
  }
  printf("\rFinished writing %s.\n", outname);
  wkspace_reset(wkspace_mark);
  while (0) {
  calc_genome_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  calc_genome_ret_OPENFAIL:
    retval = RET_OPENFAIL;
    break;
  calc_genome_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  calc_genome_ret_THREAD_CREATE_FAIL:
    retval = RET_THREAD_CREATE_FAIL;
    // todo: thread join, etc.
    break;
  }
 calc_genome_ret_1:
  gzclose_cond(gz_outfile);
  fclose_cond(outfile);
  return retval;
}

inline int distance_wt_req(int calculation_type) {
  return ((calculation_type & CALC_DISTANCE_MASK) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && ((calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE))));
}

inline int distance_req(int calculation_type) {
  return ((calculation_type & CALC_DISTANCE_MASK) || ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) && (!(calculation_type & CALC_GENOME))) || ((!(calculation_type & CALC_LOAD_DISTANCES)) && ((calculation_type & CALC_GROUPDIST) || (calculation_type & CALC_REGRESS_DISTANCE))));
}

inline int relationship_req(int calculation_type) {
  return (calculation_type & (CALC_RELATIONSHIP_MASK | CALC_UNRELATED_HERITABILITY | CALC_REL_CUTOFF | CALC_REGRESS_REL));
}

inline int relationship_or_ibc_req(int calculation_type) {
  return (relationship_req(calculation_type) || (calculation_type & CALC_IBC));
}

int wdist(char* outname, char* pedname, char* mapname, char* famname, char* phenoname, char* extractname, char* excludename, char* keepname, char* removename, char* filtername, char* freqname, char* loaddistname, char* makepheno_str, char* filterval, int mfilter_col, int ped_col_1, int ped_col_34, int ped_col_5, int ped_col_6, int ped_col_7, char missing_geno, int missing_pheno, int mpheno_col, char* phenoname_str, int prune, int affection_01, unsigned int chrom_mask, double exponent, double min_maf, double geno_thresh, double mind_thresh, double hwe_thresh, double rel_cutoff, int tail_pheno, double tail_bottom, double tail_top, int calculation_type, int groupdist_iters, int groupdist_d, int regress_iters, int regress_d, int regress_rel_iters, int regress_rel_d, double unrelated_herit_tol, double unrelated_herit_covg, double unrelated_herit_covr, int ibc_type, int parallel_idx, int parallel_tot, int ppc_gap, int nonfounders, int genome_output_gz, int genome_output_full, int genome_ibd_unbounded) {
  FILE* outfile = NULL;
  FILE* outfile2 = NULL;
  FILE* outfile3 = NULL;
  FILE* mapfile = NULL;
  FILE* pedfile = NULL;
  FILE* famfile = NULL;
  gzFile gz_outfile = NULL;
  gzFile gz_outfile2 = NULL;
  gzFile gz_outfile3 = NULL;
  FILE* phenofile = NULL;
  FILE* filterfile = NULL;
  FILE* bedtmpfile = NULL;
  FILE* bimtmpfile = NULL;
  FILE* famtmpfile = NULL;
  FILE* freqfile = NULL;
  FILE* loaddistfile = NULL;
  char* id_buf = NULL;
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
  int plink_maxsnp;
  double* mafs = NULL;
  int ped_recalc_len = 0;
  char* fgets_return;
  int unfiltered_indiv_ct = 0;
  int unfiltered_indiv_ct4 = 0;
  unsigned char* indiv_exclude = NULL;
  int indiv_exclude_ct = 0;
  unsigned char* founder_info = NULL;
  int ii;
  int jj = 0;
  int kk = 0;
  int mm;
  int nn = 0;
  int oo = 0;
  int pp;
  int qq;
  int rr;
  unsigned int uii = 0;
  unsigned int ujj;
  unsigned long ulii;
  unsigned long uljj;
  unsigned long ulkk;
  double dxx;
  double dyy;
  double dzz;
  double dww;
  double dvv;
  double duu;
  long long llxx = 0;
  long long llyy;
  char* marker_ids = NULL;
  int* marker_chroms = NULL;
  int marker_pos_start = -1;
  int* marker_allele_cts = NULL;
  char* marker_alleles = NULL;
  char* marker_alleles_tmp = NULL;
  int* person_missing_cts = NULL;
  char* bufptr;
  int retval = RET_SUCCESS;
  int map_cols = 3;
  int affection = 0;
  double* phenor_d = NULL;
  char* phenor_c = NULL;
  char* person_id = NULL;
  int max_person_id_len = 4;
  char* paternal_id = NULL;
  int max_paternal_id_len = 2;
  char* maternal_id = NULL;
  int max_maternal_id_len = 2;
  unsigned char* wkspace_mark = NULL;
  unsigned int* giptr = NULL;
  unsigned int* giptr2 = NULL;
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
  double* dptr5;
  double* rel_ibc;
  long long last_tell;
  int binary_files = 0;
  int maf_int_thresh;
  int geno_int_thresh;
  int mind_int_thresh;
  int marker_exclude_ct = 0;
  int pheno_line_ct = 0;
  int makepheno_all = 0;
  int filter_line_ct = 0;
  int max_pid_len = 0;
  int max_id_len = 0;
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
  int filter_type = FILTER_NONE;
  int wt_needed = 0;
  int unwt_needed = 0;
  int unwt_needed_full = 0;
  int bed_offset = 3;

  ii = missing_pheno;
  if (ii < 0) {
    ii = -ii;
    missing_pheno_len++;
  }
  while (ii > 9) {
    ii /= 10;
    missing_pheno_len++;
  }
  if (calculation_type & CALC_RELATIONSHIP_COV) {
    var_std = 0;
    ibc_type = -1;
  }

  if (fopen_checked(&pedfile, pedname, "r")) {
    goto wdist_ret_OPENFAIL;
  }
  if (famname[0]) {
    binary_files = 1;
    if (fopen_checked(&famfile, famname, "r")) {
      goto wdist_ret_OPENFAIL;
    }
  }

  if (phenoname[0] && fopen_checked(&phenofile, phenoname, "r")) {
    goto wdist_ret_OPENFAIL;
  }
  outname_end = outname;
  while (*outname_end) {
    outname_end++;
  }

  // note that this actually only allocates up to half of the main workspace to
  // the distance array, because of the difference between n(n-1)/2 and n^2.
  if (exp0) {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(int));
  } else {
    max_people = (int)sqrt((malloc_size_mb * 1048576) / sizeof(double));
  }

  tbuf[MAXLINELEN - 6] = ' ';
  tbuf[MAXLINELEN - 1] = ' ';
  retval = load_map_or_bim(&mapfile, mapname, binary_files, &map_cols, &unfiltered_marker_ct, &marker_exclude_ct, &max_marker_id_len, &plink_maxsnp, &marker_exclude, &mafs, &marker_allele_cts, &marker_alleles, &marker_ids, &marker_chroms, chrom_mask, &marker_pos_start, *extractname, *excludename, *freqname, ppc_gap, calculation_type);
  if (retval) {
    goto wdist_ret_2;
  }

  if (!(calculation_type & (CALC_GDISTANCE_MASK | CALC_LOAD_DISTANCES | CALC_RELATIONSHIP_MASK | CALC_IBC | CALC_FREQ | CALC_WRITE_SNPLIST | CALC_GENOME))) {
    prune = 1;
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
      ii = strlen_se(tbuf);
      cptr = next_item(tbuf);
      if (no_more_items(cptr)) {
        printf(errstr_phenotype_format);
        goto wdist_ret_INVALID_FORMAT;
      }
      bufptr = cptr;
      jj = strlen_se(cptr);
      if (ii + jj + 2 > max_pid_len) {
        max_pid_len = ii + jj + 2;
      }
      if (!nn) {
        nn = 1;
        if ((ii == 3) && (jj == 3) && (!memcmp("FID", tbuf, 3)) && (!memcmp("IID", cptr, 3))) {
          pheno_line_ct = -1;
          if (phenoname_str) {
            jj = strlen(phenoname_str);
            do {
              cptr = next_item(cptr);
              if (no_more_items(cptr)) {
		printf("Error: --pheno-name column not found.\n");
		goto wdist_ret_INVALID_FORMAT;
              }
              mpheno_col++;
              ii = strlen_se(cptr);
            } while ((ii != jj) || memcmp(cptr, phenoname_str, ii));
            phenoname_str = NULL;
          }
        } else if (phenoname_str) {
          printf("Error: No header row in phenotype file.\n");
          goto wdist_ret_INVALID_FORMAT;
        }
        if (!mpheno_col) {
          mpheno_col = 1;
        }
      }
      if (pheno_line_ct == -1) {
        pheno_line_ct = 0;
      } else {
        for (ii = 0; ii < mpheno_col; ii++) {
          bufptr = next_item(bufptr);
        }
        if (no_more_items(bufptr)) {
	  printf(errstr_phenotype_format);
	  goto wdist_ret_INVALID_FORMAT;
        }
        if (makepheno_str && (!makepheno_all)) {
          if (!strncmp(makepheno_str, bufptr, kk)) {
            if (is_space_or_eoln(bufptr[kk])) {
              pheno_line_ct++;
            }
          }
        } else {
          if (!mm) {
            affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
            mm = 1;
          }
          if (affection) {
            if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
              pheno_line_ct++;
            }
          } else {
	    if (sscanf(bufptr, "%lg", &dxx) != 1) {
	      printf(errstr_phenotype_format);
	      goto wdist_ret_INVALID_FORMAT;
	    }
            if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
              pheno_line_ct++;
            }
          }
        }
      }
      if (!tbuf[MAXLINELEN - 1]) {
        printf("Error: Excessively long line in phenotype file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_INVALID_FORMAT;
      }
    }
    if (!feof(phenofile)) {
      goto wdist_ret_READ_FAIL;
    }
    id_buf = (char*)malloc(max_pid_len * sizeof(char));
    if (!id_buf) {
      goto wdist_ret_NOMEM;
    }
    if (pheno_line_ct) {
      pid_list = (char*)malloc(max_pid_len * pheno_line_ct * sizeof(char));
      if (!pid_list) {
	goto wdist_ret_NOMEM;
      }
      if (affection || tail_pheno) {
	phenor_c = (char*)malloc(pheno_line_ct * sizeof(char));
	if (!phenor_c) {
	  goto wdist_ret_NOMEM;
	}
      } else if (!makepheno_str) {
	phenor_d = (double*)malloc(pheno_line_ct * sizeof(double));
	if (!phenor_d) {
	  goto wdist_ret_NOMEM;
	}
      }

      rewind(phenofile);
      ii = 0;
      // ----- phenotype file load, second pass -----
      while (fgets(tbuf, MAXLINELEN, phenofile) != NULL) {
	if (ii == pheno_line_ct) {
	  break;
	}
	if (*tbuf == '\n') {
	  continue;
	}
        cptr = read_next_unsafe(&(pid_list[ii * max_pid_len]), tbuf);
	*cptr++ = '\t';
        bufptr = next_item(tbuf);
        read_next_terminate(cptr, bufptr);
	for (jj = 0; jj < mpheno_col; jj++) {
	  bufptr = next_item(bufptr);
	}
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
	if (qsort_ext(pid_list, pheno_line_ct, max_pid_len, strcmp_deref, phenor_c, sizeof(char)) == -1) {
	  goto wdist_ret_NOMEM;
        }
      } else if (makepheno_str) {
        qsort(pid_list, pheno_line_ct, max_pid_len, strcmp_casted);
      } else {
	if (qsort_ext(pid_list, pheno_line_ct, max_pid_len, strcmp_deref, (char*)phenor_d, sizeof(double)) == -1) {
          goto wdist_ret_NOMEM;
        }
      }
    } else {
      printf("Note: No valid entries in phenotype file.\n");
    }
    fclose_null(&phenofile);
  }

  if (extractname[0] || excludename[0]) {
    wkspace_mark = wkspace_base;
    ii = unfiltered_marker_ct - marker_exclude_ct;
    if (sort_marker_ids(&cptr, &iptr, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len)) {
      goto wdist_ret_NOMEM;
    }
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

  if (keepname[0] || filtername[0]) {
    if (keepname[0]) {
      cptr = keepname;
    } else {
      cptr = filtername;
      jj = strlen(filterval);
      if (!mfilter_col) {
        mfilter_col = 1;
      }
    }
    if (fopen_checked(&filterfile, cptr, "r")) {
      goto wdist_ret_OPENFAIL;
    }
    filter_type = FILTER_KEEP;
  }

  if (filterfile) {
    // ----- --keep/--filter load, two passes -----
    max_id_len = determine_max_id_len(filterfile, filterval, mfilter_col, &filter_line_ct);
    if (max_id_len < 0) {
      goto wdist_ret_INVALID_FORMAT;
    }
    id_buf = resize_id_buf(id_buf, max_id_len, max_pid_len);
    if (!id_buf) {
      goto wdist_ret_NOMEM;
    }
    if (filter_line_ct) {
      id_list = (char*)malloc(max_id_len * filter_line_ct * sizeof(char));
      if (!id_list) {
	goto wdist_ret_NOMEM;
      }
      rewind(filterfile);
      ii = 0;
      while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
	if (*tbuf == '\n') {
	  continue;
	}
        cptr = read_next_unsafe(&(id_list[ii * max_id_len]), tbuf);
	*cptr++ = '\t';
        bufptr = next_item(tbuf);
        read_next_terminate(cptr, bufptr);
        if (filtername[0]) {
	  for (kk = 0; kk < mfilter_col; kk++) {
	    bufptr = next_item(bufptr);
	  }
	  if (strncmp(filterval, bufptr, jj)) {
            continue;
	  }
        }
	ii++;
	if (ii == filter_line_ct) {
	  break;
	}
      }
      qsort(id_list, filter_line_ct, max_id_len, strcmp_casted);
    } else {
      printf("Note: No valid entries in filter file.\n");
    }
    fclose_null(&filterfile);
  }
  if (removename[0]) {
    if (fopen_checked(&filterfile, removename, "r")) {
      goto wdist_ret_OPENFAIL;
    }
    if (filter_type == FILTER_KEEP) {
      // instead of building new list, just remove items from --keep/--filter
      // list
      jj = 0; // items removed
      while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
        if (*tbuf == '\n') {
          continue;
        }
        cptr = next_item(tbuf);
        if (no_more_items(cptr)) {
          continue;
        }
        ii = bsearch_fam_indiv(id_buf, id_list, max_id_len, filter_line_ct, tbuf, cptr);
        if (ii != -1) {
          // mark item for removal
          id_list[ii * max_id_len] = '\0';
          jj++;
        }
      }
      if (jj) {
        ii = 0;
        while (id_list[ii * max_id_len]) {
          ii++;
        }
        filter_line_ct -= jj;
        kk = ii + 1;
        while (ii < jj) {
          if (id_list[kk * max_id_len]) {
            memcpy(&(id_list[ii * max_id_len]), &(id_list[kk * max_id_len]), max_id_len);
            ii++;
          }
          kk++;
        }
      }
    } else {
      filter_type = FILTER_REMOVE;
      // ----- --remove load without --keep/--filter, two passes -----
      max_id_len = determine_max_id_len(filterfile, NULL, 0, &filter_line_ct);
      if (max_id_len < 0) {
        goto wdist_ret_INVALID_FORMAT;
      }
      id_buf = resize_id_buf(id_buf, max_id_len, max_pid_len);
      if (!id_buf) {
	goto wdist_ret_NOMEM;
      }
      if (filter_line_ct) {
	id_list = (char*)malloc(max_id_len * filter_line_ct * sizeof(char));
	if (!id_list) {
	  goto wdist_ret_NOMEM;
	}
        rewind(filterfile);
        ii = 0;
	// ----- --remove load without --keep/--filter, second pass -----
	while (fgets(tbuf, MAXLINELEN, filterfile) != NULL) {
	  if (*tbuf == '\n') {
	    continue;
	  }
          cptr = read_next_unsafe(&(id_list[ii * max_id_len]), tbuf);
	  *cptr++ = '\t';
          read_next_terminate(cptr, next_item(tbuf));
	  ii++;
	  if (ii == filter_line_ct) {
	    break;
          }
	}
	qsort(id_list, filter_line_ct, max_id_len, strcmp_casted);
      } else {
	printf("Note: No valid entries in --remove file.\n");
      }
    }
    fclose_null(&filterfile);
  }

  line_locs = (long long*)malloc(max_people * sizeof(long long));
  if (!line_locs) {
    goto wdist_ret_NOMEM;
  }

  // TODO: remove use of malloc'd pedbuf, do potentially large reads in wkspace
  // instead
  if (binary_files) {
    unfiltered_marker_ct4 = (unfiltered_marker_ct + 3) / 4;
    pedbuflen = unfiltered_marker_ct4;
  } else {
    pedbuflen = unfiltered_marker_ct * 4 + PEDBUFBASE;
  }
  pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
  if (!pedbuf) {
    goto wdist_ret_NOMEM;
  }

  mind_int_thresh = (int)(mind_thresh * (unfiltered_marker_ct - marker_exclude_ct));
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
            ii = is_contained(id_buf, id_list, max_id_len, filter_line_ct, tbuf, cptr);
            if (filter_type == FILTER_REMOVE) {
              ii = 1 - ii;
            }
          } else {
            ii = 1;
          }
	  if (ii && ped_col_34) {
	    cptr2 = next_item(cptr);
	    jj = strlen_se(cptr2) + 1;
	    if (jj > max_paternal_id_len) {
	      max_paternal_id_len = jj;
	    }
	    cptr2 = next_item(cptr2);
	    jj = strlen_se(cptr2) + 1;
	    if (jj > max_maternal_id_len) {
	      max_maternal_id_len = jj;
	    }
	  }
          if (phenoname[0]) {
            if (ii) {
	      if (makepheno_str || (!prune)) {
		line_locs[nn++] = unfiltered_indiv_ct;
	      } else if (is_contained(id_buf, pid_list, max_pid_len, pheno_line_ct, tbuf, cptr)) {
		line_locs[nn++] = unfiltered_indiv_ct;
	      }
            }
          } else if (ii || (!unfiltered_indiv_ct)) {
	    if (ped_col_1) {
	      bufptr = next_item(bufptr);
	    }
	    if (ped_col_34) {
	      bufptr = next_item(next_item(bufptr));
	    }
	    if (ped_col_5) {
	      bufptr = next_item(bufptr);
	    }
	    if (!bufptr) {
	      printf(errstr_fam_format);
	      goto wdist_ret_INVALID_FORMAT;
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
                  printf(errstr_fam_format);
                  goto wdist_ret_INVALID_FORMAT;
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
        printf("Error: Excessively long line in .fam file (max %d chars).\n", MAXLINELEN - 3);
        goto wdist_ret_INVALID_FORMAT;
      }
    }
    if (nn < 2) {
      printf("Error: Less than two valid people in .fam file after initial filter.\n");
      goto wdist_ret_INVALID_FORMAT;
    } else if (nn > max_people) {
      printf("Error: Too many people in .fam file for this algorithm.\n");
      goto wdist_ret_INVALID_FORMAT;
    }

    bin_pheno = (makepheno_str || affection || tail_pheno);
    if (bin_pheno) {
      pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
      if (!pheno_c) {
	goto wdist_ret_NOMEM;
      }
    } else {
      pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
      if (!pheno_d) {
	goto wdist_ret_NOMEM;
      }
    }
    person_id = (char*)malloc(unfiltered_indiv_ct * max_person_id_len * sizeof(char));
    if (!person_id) {
      goto wdist_ret_NOMEM;
    }
    if (ped_col_34) {
      paternal_id = (char*)malloc(unfiltered_indiv_ct * max_paternal_id_len * sizeof(char));
      if (!paternal_id) {
	goto wdist_ret_NOMEM;
      }
      maternal_id = (char*)malloc(unfiltered_indiv_ct * max_maternal_id_len * sizeof(char));
      if (!maternal_id) {
	goto wdist_ret_NOMEM;
      }
      founder_info = (unsigned char*)malloc(((unfiltered_indiv_ct + 7) / 8) * sizeof(char));
      if (!founder_info) {
	goto wdist_ret_NOMEM;
      }
      memset(founder_info, 0, (unfiltered_indiv_ct + 7) / 8);
    }

    indiv_exclude = (unsigned char*)calloc(sizeof(char), ((unfiltered_indiv_ct + 7) / 8));
    if (!indiv_exclude) {
      goto wdist_ret_NOMEM;
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
	  if (ped_col_34) {
	    cptr2 = next_item(cptr);
	    kk = strlen_se(cptr2);
	    memcpy(&(paternal_id[ii * max_paternal_id_len]), cptr2, kk);
	    paternal_id[ii * max_paternal_id_len + kk] = '\0';
	    cptr2 = next_item(cptr2);
	    mm = strlen_se(cptr2);
	    if ((kk == 1) && (mm == 1) && (paternal_id[ii * max_paternal_id_len] == '0') && (*cptr2 == '0')) {
	      set_founder_info(founder_info, ii);
	    }
	    memcpy(&(maternal_id[ii * max_maternal_id_len]), cptr2, mm);
	    maternal_id[ii * max_maternal_id_len + kk] = '\0';
	  }
          if (phenoname[0]) {
	    kk = bsearch_fam_indiv(id_buf, pid_list, max_pid_len, pheno_line_ct, tbuf, cptr);
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
	      bufptr = next_item(next_item(bufptr));
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
      retval = read_freq(freqname, &freqfile, unfiltered_marker_ct, marker_exclude, marker_exclude_ct, marker_ids, max_marker_id_len, marker_chroms, marker_alleles, mafs);
      if (retval) {
	goto wdist_ret_1;
      }
    }

    maf_int_thresh = 2 * (nn - indiv_exclude_ct) - (int)((1.0 - min_maf) * (nn - indiv_exclude_ct) * 2);
    geno_int_thresh = 2 * (nn - indiv_exclude_ct) - (int)(geno_thresh * 2 * (nn - indiv_exclude_ct));

    if (fread(pedbuf, 1, 3, pedfile) < 3) {
      goto wdist_ret_READ_FAIL;
    }
    if ((pedbuf[0] != 'l') || (pedbuf[1] != '\x1b') || (pedbuf[2] != '\x01')) {
      if (pedbuf[0] == '\x01') {
	bed_offset = 1;
      } else {
	strcpy(outname_end, ".bed.tmp");
	printf("Individual-major .bed file detected.  Converting to SNP-major.\n");
	fclose(pedfile);
	retval = indiv_major_to_snp_major(pedname, outname, &outfile, unfiltered_marker_ct);
	if (retval) {
	  goto wdist_ret_1;
	}
	strcpy(pedname, outname);
	if (fopen_checked(&pedfile, pedname, "rb")) {
	  goto wdist_ret_OPENFAIL;
	}
	fseeko(pedfile, 3, SEEK_SET);
      }
    }

    // ----- snp-major .bed load, first pass -----
    person_missing_cts = (int*)calloc(sizeof(int), unfiltered_indiv_ct);
    if (!person_missing_cts) {
      goto wdist_ret_NOMEM;
    }
    unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
    if (pedbuflen < (4 * unfiltered_indiv_ct4)) {
      free(pedbuf);
      pedbuflen = 4 * unfiltered_indiv_ct4;
      pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
      if (!pedbuf) {
	goto wdist_ret_NOMEM;
      }
    }
    for (ii = 0; ii < unfiltered_marker_ct; ii++) {
      if (excluded(marker_exclude, ii)) {
	fseeko(pedfile, (off_t)unfiltered_indiv_ct4, SEEK_CUR);
	continue;
      }
      if (fread(pedbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	goto wdist_ret_READ_FAIL;
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
      goto wdist_ret_NOMEM;
    }
    dptr2 = marker_weights;

    fseeko(pedfile, bed_offset, SEEK_SET);
    // ----- .snp-major .bed load, second pass -----
    for (ii = 0; ii < unfiltered_marker_ct; ii++) {
      if (excluded(marker_exclude, ii)) {
	fseeko(pedfile, (off_t)unfiltered_indiv_ct4, SEEK_CUR);
	continue;
      }
      if (fread(pedbuf, 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	goto wdist_ret_READ_FAIL;
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
	if ((!excluded(indiv_exclude, jj)) && (nonfounders || is_founder(founder_info, jj))) {
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
    if (calculation_type & CALC_GENOME) {
      marker_pos_start = collapse_marker_pos(marker_pos, marker_exclude, unfiltered_marker_ct, marker_exclude_ct, marker_pos_start);
    }
  } else {
    if (freqname[0]) {
      printf("--read-freq currently ignored on text file load.\n");
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
	ii = strlen_se((char*)pedbuf) + strlen_se(cptr) + 2;
	if (ii > max_person_id_len) {
	  max_person_id_len = ii;
	}
	if (filter_type) {
	  ii = is_contained(id_buf, id_list, max_id_len, filter_line_ct, (char*)pedbuf, cptr);
	  if (filter_type == FILTER_REMOVE) {
	    ii = 1 - ii;
	  }
	} else {
	  ii = 1;
	}
	if (ii && ped_col_34) {
	  cptr2 = next_item(cptr);
	  jj = strlen_se(cptr2) + 1;
	  if (jj > max_paternal_id_len) {
	    max_paternal_id_len = jj;
	  }
	  cptr2 = next_item(cptr2);
	  jj = strlen_se(cptr2) + 1;
	  if (jj > max_maternal_id_len) {
	    max_maternal_id_len = jj;
	  }
	}
	if (phenoname[0]) {
	  if (ii) {
	    if (makepheno_str || (!prune)) {
	      line_locs[unfiltered_indiv_ct] = 1;
	      nn++;
	    } else {
	      if (is_contained(id_buf, pid_list, max_pid_len, pheno_line_ct, (char*)pedbuf, cptr)) {
		line_locs[unfiltered_indiv_ct] = 1;
		nn++;
	      }
	    }
	  }
	} else if (ii || (!unfiltered_indiv_ct)) {
	  if (ped_col_1) { // actually column 2
	    bufptr = next_item(bufptr);
	  }
	  if (ped_col_34) {
	    bufptr = next_item(next_item(bufptr));
	  }
	  if (ped_col_5) {
	    bufptr = next_item(bufptr);
	  }
	  if (!bufptr) {
	    printf(errstr_ped_format);
	    goto wdist_ret_INVALID_FORMAT;
	  }
	  if (!unfiltered_indiv_ct) {
	    affection = eval_affection(bufptr, missing_pheno, missing_pheno_len, affection_01);
	  }
	  if (ii) {
	    if (affection) {
	      if ((!prune) || (!is_missing(bufptr, missing_pheno, missing_pheno_len, affection_01))) {
		line_locs[unfiltered_indiv_ct] = 1;
		nn++;
	      }
	    } else {
	      if (sscanf(bufptr, "%lg", &dxx) != 1) {
		printf(errstr_ped_format);
		goto wdist_ret_INVALID_FORMAT;
	      }
	      if ((!prune) || ((dxx != missing_phenod) && ((!tail_pheno) || ((dxx <= tail_bottom) || (dxx > tail_top))))) {
		line_locs[unfiltered_indiv_ct] = 1;
		nn++;
	      }
	    }
	  }
	}
        unfiltered_indiv_ct++;
      }
      if (!pedbuf[pedbuflen - 1]) {
	// determine extended read buffer length needed to handle unexpectedly
        // long line
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
      printf("Error: Less than two valid people in .ped file.\n");
      goto wdist_ret_INVALID_FORMAT;
    } else if (unfiltered_indiv_ct > max_people) {
      printf("Error: Too many people in .ped file for this algorithm.\n");
      goto wdist_ret_INVALID_FORMAT;
    }
    rewind(pedfile);
    if (ped_recalc_len) {
      free(pedbuf);
      pedbuflen = ped_recalc_len;
      pedbuf = (unsigned char*)malloc(pedbuflen * sizeof(char));
      if (!pedbuf) {
	goto wdist_ret_NOMEM;
      }
    }

    person_id = (char*)malloc(unfiltered_indiv_ct * max_person_id_len * sizeof(char));
    if (!person_id) {
      goto wdist_ret_NOMEM;
    }
    if (ped_col_34) {
      paternal_id = (char*)malloc(unfiltered_indiv_ct * max_paternal_id_len * sizeof(char));
      if (!paternal_id) {
	goto wdist_ret_NOMEM;
      }
      maternal_id = (char*)malloc(unfiltered_indiv_ct * max_maternal_id_len * sizeof(char));
      if (!maternal_id) {
	goto wdist_ret_NOMEM;
      }
      founder_info = (unsigned char*)malloc(((unfiltered_indiv_ct + 7) / 8) * sizeof(char));
      if (!founder_info) {
	goto wdist_ret_NOMEM;
      }
      memset(founder_info, 0, (unfiltered_indiv_ct + 7) / 8);
    }

    marker_alleles = (char*)calloc(sizeof(char), unfiltered_marker_ct * 4);
    if (!marker_alleles) {
      goto wdist_ret_NOMEM;
    }
    marker_allele_cts = (int*)calloc(sizeof(int), unfiltered_marker_ct * 4);
    if (!marker_allele_cts) {
      goto wdist_ret_NOMEM;
    }
    bin_pheno = (makepheno_str || affection || tail_pheno);
    if (bin_pheno) {
      pheno_c = (char*)malloc(unfiltered_indiv_ct * sizeof(char));
      if (!pheno_c) {
	goto wdist_ret_NOMEM;
      }
    } else {
      pheno_d = (double*)malloc(unfiltered_indiv_ct * sizeof(double));
      if (!pheno_d) {
	goto wdist_ret_NOMEM;
      }
    }
    indiv_exclude = (unsigned char*)calloc(sizeof(char), ((unfiltered_indiv_ct + 7) / 8));
    if (!indiv_exclude) {
      goto wdist_ret_NOMEM;
    }

    // ----- .ped load, second pass -----
    kk = 0;
    for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
      do {
        last_tell = ftello(pedfile);
        fgets((char*)pedbuf, pedbuflen, pedfile);
      } while ((pedbuf[0] <= ' ') || (pedbuf[0] == '#'));
      if (!line_locs[ii]) {
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
      if (ped_col_34) {
        cptr = next_item(bufptr);
	oo = strlen_se(cptr);
	memcpy(&(paternal_id[ii * max_paternal_id_len]), cptr, oo);
	paternal_id[ii * max_paternal_id_len + oo] = '\0';
	cptr2 = next_item(cptr);
	pp = strlen_se(cptr2);
	if ((oo == 1) && (pp == 1) && (*cptr == '0') && (*cptr2 == '0')) {
	  set_founder_info(founder_info, ii);
	}
	memcpy(&(maternal_id[ii * max_maternal_id_len]), cptr2, pp);
	maternal_id[ii * max_maternal_id_len + pp] = '\0';
      }
      if (phenoname[0]) {
        jj = bsearch_fam_indiv(id_buf, pid_list, max_pid_len, pheno_line_ct, (char*)pedbuf, bufptr);
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
        bufptr = next_item(bufptr);
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
        printf(errstr_ped_format);
        goto wdist_ret_INVALID_FORMAT;
      }
      line_locs[ii] = last_tell + (bufptr - (char*)pedbuf);
      mm = 0; // number of missing markers not previously filtered out
      for (jj = 0; jj < unfiltered_marker_ct; jj += 1) {
        for (nn = 0; nn < 2; nn++) {
          cc = *bufptr;
	  if (cc == '\0') {
	    printf(errstr_ped_format);
	    goto wdist_ret_INVALID_FORMAT;
	  }
          if (!excluded(marker_exclude, jj)) {
	    if (cc == missing_geno) {
	      if (nn) {
		printf("Error: Second allele marked as missing when first was not.\n(individual %d, marker %d)\n", ii, jj);
		goto wdist_ret_INVALID_FORMAT;
	      }
	      mm++;
	      nn++;
	      bufptr++;
	      while ((*bufptr == ' ') || (*bufptr == '\t')) {
		bufptr++;
	      }
	    } else if (nonfounders || is_founder(founder_info, ii)) {
	      if (cc == marker_alleles[jj * 4]) {
		marker_allele_cts[jj * 4] += 1;
	      } else if (cc == marker_alleles[jj * 4 + 1]) {
		marker_allele_cts[jj * 4 + 1] += 1;
	      } else if (cc == marker_alleles[jj * 4 + 2]) {
		marker_allele_cts[jj * 4 + 2] += 1;
	      } else if (cc == marker_alleles[jj * 4 + 3]) {
		marker_allele_cts[jj * 4 + 3] += 1;
	      } else if (marker_alleles[jj * 4 + 3]) {
		printf("Error: More than 4 different allele types at marker %d.\n", jj + 1);
		goto wdist_ret_INVALID_FORMAT;
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
    // to match PLINK, --geno applies after --mind filter
    geno_int_thresh = (unfiltered_indiv_ct - indiv_exclude_ct) * 2 - (int)(geno_thresh * (unfiltered_indiv_ct - indiv_exclude_ct) * 2);
    for (ii = 0; ii < unfiltered_marker_ct; ii++) {
      for (jj = 1; jj < 4; jj++) {
        for (kk = jj - 1; kk >= 0; kk--) {
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
      maf_int_thresh = mm - (int)((1.0 - min_maf) * mm);

      if (marker_allele_cts[ii * 4] + marker_allele_cts[ii * 4 + 1] < geno_int_thresh) {
        exclude(marker_exclude, ii, &marker_exclude_ct);
      } else if (marker_allele_cts[ii * 4 + 1] < maf_int_thresh) {
        exclude(marker_exclude, ii, &marker_exclude_ct);
      } else if (!marker_allele_cts[ii * 4 + 1]) {
        marker_alleles[ii * 4 + 1] = missing_geno;
      } else if (marker_allele_cts[ii * 4 + 1] == marker_allele_cts[ii * 4 + 2]) {
        printf("Error: Ambiguous minor allele at marker %d.\n", ii + 1);
        goto wdist_ret_INVALID_FORMAT;
      }
    }
    if (indiv_exclude_ct > (unfiltered_indiv_ct - 2)) {
      printf("Error: Too many people fail QC.\n");
      goto wdist_ret_INVALID_FORMAT;
    } else if (marker_exclude_ct == unfiltered_marker_ct) {
      printf("Error: All markers fail QC.\n");
      goto wdist_ret_INVALID_FORMAT;
    }

    het_probs = (double*)malloc(((unfiltered_marker_ct - marker_exclude_ct) / 2 + 1) * sizeof(double));
    if (!het_probs) {
      goto wdist_ret_NOMEM;
    }
    dptr2 = marker_weights;

    hwe_ll = (int*)wkspace_alloc(unfiltered_marker_ct * 9 * sizeof(int));
    if (!hwe_ll) {
      goto wdist_ret_NOMEM;
    }
    // unfortunately, this requires a third pass...
    hwe_lh = &(hwe_ll[unfiltered_marker_ct]);
    hwe_hh = &(hwe_ll[unfiltered_marker_ct * 2]);
    hwe_u_ll = &(hwe_ll[unfiltered_marker_ct * 3]);
    hwe_u_lh = &(hwe_ll[unfiltered_marker_ct * 4]);
    hwe_u_hh = &(hwe_ll[unfiltered_marker_ct * 5]);
    hwe_a_ll = &(hwe_ll[unfiltered_marker_ct * 6]);
    hwe_a_lh = &(hwe_ll[unfiltered_marker_ct * 7]);
    hwe_a_hh = &(hwe_ll[unfiltered_marker_ct * 8]);
    fill_int_zero(hwe_ll, unfiltered_marker_ct * 9);

    for (ii = 0; ii < unfiltered_indiv_ct; ii += 1) {
      if (excluded(indiv_exclude, ii) || ((!nonfounders) && (!is_founder(founder_info, ii)))) {
	continue;
      }
      fseeko(pedfile, line_locs[ii], SEEK_SET);
      fgets((char*)pedbuf, pedbuflen, pedfile);
      bufptr = (char*)pedbuf;
      for (jj = 0; jj < unfiltered_marker_ct; jj += 1) {
	if (excluded(marker_exclude, jj)) {
          bufptr = next_item(next_item(bufptr));
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
    wkspace_reset(hwe_ll);
    rewind(pedfile);

    ulii = unfiltered_indiv_ct - indiv_exclude_ct;
    ulii = (ulii * (ulii - 1)) / 2;
    dists_alloc = ulii * (sizeof(int) + sizeof(double));
    
    llxx = malloc_size_mb * 1048576 - dists_alloc;
    geno_window_size = llxx / (unfiltered_indiv_ct - indiv_exclude_ct);
    unfiltered_marker_ct4 = (unfiltered_marker_ct - marker_exclude_ct + 3) / 4;
    uljj = (unfiltered_indiv_ct - indiv_exclude_ct + 3) / 4;
    ped_geno = wkspace_base;
    if (uljj <= (wkspace_left / (unfiltered_marker_ct - marker_exclude_ct))) {
      printf("Writing binary files...");
      strcpy(outname_end, ".bed");
      if (fopen_checked(&bedtmpfile, outname, "wb")) {
	goto wdist_ret_OPENFAIL;
      }
      strcpy(outname_end, ".fam");
      if (fopen_checked(&famtmpfile, outname, "wb")) {
	goto wdist_ret_OPENFAIL;
      }
      if (fwrite_checked("l\x1b\x01", 3, bedtmpfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
      if (wkspace_left < uljj * (unfiltered_marker_ct - marker_exclude_ct)) {
        goto wdist_ret_NOMEM;
      }
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
        if (fwrite_checked(pedbuf, jj, famtmpfile)) {
          goto wdist_ret_WRITE_FAIL;
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
      if (fwrite_checked(ped_geno, ulii, bedtmpfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
      if (fclose_null(&bedtmpfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
    } else {
      printf("Using two-step .ped -> .bed conversion, since .ped is very large...");
      strcpy(outname_end, ".fam");
      if (fopen_checked(&famtmpfile, outname, "wb")) {
	goto wdist_ret_OPENFAIL;
      }
      strcpy(outname_end, ".bed.tmp");
      if (fopen_checked(&bedtmpfile, outname, "wb")) {
	goto wdist_ret_OPENFAIL;
      }
      if (fwrite_checked("l\x1b", 3, bedtmpfile)) {
	goto wdist_ret_WRITE_FAIL;
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
        if (fwrite_checked(pedbuf, jj, famtmpfile)) {
          goto wdist_ret_WRITE_FAIL;
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
	if (fwrite_checked(ped_geno, unfiltered_marker_ct4, bedtmpfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
      }
      if (fclose_null(&bedtmpfile)) {
	goto wdist_ret_WRITE_FAIL;
      }

      strcpy(tbuf, outname);
      strcpy(outname_end, ".bed");
      retval = indiv_major_to_snp_major(tbuf, outname, &outfile, unfiltered_marker_ct);
      if (retval) {
	goto wdist_ret_1;
      }
    }

    strcpy(outname_end, ".bim");
    if (fopen_checked(&bimtmpfile, outname, "wb")) {
      goto wdist_ret_OPENFAIL;
    }

    marker_alleles_tmp = (char*)malloc((unfiltered_marker_ct - marker_exclude_ct) * 2 * sizeof(char));
    if (!marker_alleles_tmp) {
      goto wdist_ret_NOMEM;
    }
    cptr = marker_alleles_tmp;
    iwptr = (int*)malloc((unfiltered_marker_ct - marker_exclude_ct) * 2 * sizeof(int));
    iptr = iwptr;
    dptr2 = mafs;
    for (ii = 0; ii < unfiltered_marker_ct; ii += 1) {
      if (excluded(marker_exclude, ii)) {
	continue;
      }
      *cptr++ = marker_alleles[ii * 4 + 1];
      *cptr++ = marker_alleles[ii * 4];
      *iptr++ = marker_allele_cts[ii * 4 + 1];
      *iptr++ = marker_allele_cts[ii * 4];
      *dptr2++ = mafs[ii];
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
      *bufptr++ = marker_alleles[kk * 2];
      *bufptr++ = '\t';
      *bufptr++ = marker_alleles[kk * 2 + 1];
      *bufptr++ = '\n';
      kk++;
      jj = (int)(bufptr - tbuf);
      if (fwrite_checked(tbuf, jj, bimtmpfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
    }
    printf(" done.\n");
    fclose(pedfile);
    fclose_null(&bimtmpfile);
    fclose_null(&famtmpfile);
    strcpy(outname_end, ".bed");
    if (fopen_checked(&pedfile, outname, "rb")) {
      goto wdist_ret_OPENFAIL;
    }
    binary_files = 1;
    if (pheno_c) {
      collapse_arr(pheno_c, sizeof(char), indiv_exclude, unfiltered_indiv_ct);
    } else if (pheno_d) {
      collapse_arr((char*)pheno_d, sizeof(double), indiv_exclude, unfiltered_indiv_ct);
    }
    if (calculation_type & (CALC_FREQ | CALC_WRITE_SNPLIST)) {
      collapse_arr(marker_ids, max_marker_id_len, marker_exclude, unfiltered_marker_ct);
      collapse_arr((char*)marker_chroms, sizeof(int), marker_exclude, unfiltered_marker_ct);
    }
    collapse_arr(person_id, max_person_id_len, indiv_exclude, unfiltered_indiv_ct);
    if (ped_col_34) {
      collapse_arr(paternal_id, max_paternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_arr(maternal_id, max_maternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_bitarr(founder_info, indiv_exclude, unfiltered_indiv_ct);
    }
    collapse_arr((char*)marker_weights, sizeof(double), marker_exclude, unfiltered_marker_ct);
    if (calculation_type & CALC_GENOME) {
      marker_pos_start = collapse_marker_pos(marker_pos, marker_exclude, unfiltered_marker_ct, marker_exclude_ct, marker_pos_start);
    }
    unfiltered_marker_ct -= marker_exclude_ct;
    marker_exclude_ct = 0;
    memset(marker_exclude, 0, (unfiltered_marker_ct + 7) / 8);
    unfiltered_indiv_ct -= indiv_exclude_ct;
    unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
    indiv_exclude_ct = 0;
    memset(indiv_exclude, 0, (unfiltered_indiv_ct + 7) / 8);
  }

  if ((calculation_type & CALC_GROUPDIST) && (!bin_pheno)) {
    printf("Error: --groupdist calculation requires binary phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_REGRESS_DISTANCE) && (!pheno_d)) {
    printf("Error: --regress-distance calculation requires scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  } else if ((calculation_type & CALC_UNRELATED_HERITABILITY) && (!pheno_d)) {
    printf("Error: --unrelated-heritability requires scalar phenotype.\n");
    goto wdist_ret_INVALID_CMDLINE;
  }

  marker_ct = unfiltered_marker_ct - marker_exclude_ct;
  indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;

  printf("%d markers and %d people pass filters and QC%s.\n", marker_ct, indiv_ct, (calculation_type & CALC_REL_CUTOFF)? " (before --rel-cutoff)": "");
  if (parallel_tot > indiv_ct / 2) {
    printf("Error: Too many --parallel jobs (maximum %d/2 = %d).\n", indiv_ct, indiv_ct / 2);
    goto wdist_ret_INVALID_CMDLINE;
  }
  if (thread_ct > 1) {
    printf("Using %d threads (change this with --threads).\n", thread_ct);
  }

  if (calculation_type & CALC_WRITE_SNPLIST) {
    strcpy(outname_end, ".snplist");
    if (fopen_checked(&outfile, outname, "w")) {
      goto wdist_ret_OPENFAIL;
    }
    for (ii = 0; ii < unfiltered_marker_ct; ii++) {
      if (!excluded(marker_exclude, ii)) {
        if (fprintf(outfile, "%s\n", &(marker_ids[ii * max_marker_id_len])) < 0) {
          goto wdist_ret_WRITE_FAIL;
        }
      }
    }
    if (fclose_null(&outfile)) {
      goto wdist_ret_WRITE_FAIL;
    }
  }

  if (calculation_type & CALC_FREQ) {
    strcpy(outname_end, ".frq");
    if (fopen_checked(&outfile, outname, "w")) {
      goto wdist_ret_OPENFAIL;
    }
    if (plink_maxsnp < 5) {
      if (fprintf(outfile, " CHR  SNP   A1   A2          MAF  NCHROBS\n") < 0) {
	goto wdist_ret_WRITE_FAIL;
      }
      strcpy(tbuf, "%4d %4s    %c    %c       %6.4g  %7d\n");
    } else {
      sprintf(tbuf, " CHR %%%ds   A1   A2       MAF  NCHROBS\n", plink_maxsnp);
      if (fprintf(outfile, tbuf, "SNP") < 0) {
	goto wdist_ret_WRITE_FAIL;
      }
      sprintf(tbuf, "%%4d %%%ds    %%c    %%c    %%6.4g  %%7d\n", plink_maxsnp);
    }
    for (ii = 0; ii < unfiltered_marker_ct; ii++) {
      if (excluded(marker_exclude, ii)) {
	continue;
      }
      if (mafs[ii] < 0.5) {
	if (fprintf(outfile, tbuf, marker_chroms[ii], &(marker_ids[ii * max_marker_id_len]), marker_alleles[ii * 2 + 1], marker_alleles[ii * 2], mafs[ii], marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]) < 0) {
	  goto wdist_ret_WRITE_FAIL;
	}
      } else {
	if (fprintf(outfile, tbuf, marker_chroms[ii], &(marker_ids[ii * max_marker_id_len]), marker_alleles[ii * 2], marker_alleles[ii * 2 + 1], 1.0 - mafs[ii], marker_allele_cts[ii * 2] + marker_allele_cts[ii * 2 + 1]) < 0) {
	  goto wdist_ret_WRITE_FAIL;
	}
      }
    }
    if (fclose_null(&outfile)) {
      goto wdist_ret_WRITE_FAIL;
    }
    printf("Allele frequencies written to %s.\n", outname);
  }

  wkspace_reset(marker_weights);
  wt_needed = distance_wt_req(calculation_type);
  if (wt_needed) {
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
    marker_weights_i = (unsigned int*)wkspace_alloc(marker_ct * sizeof(int));
  }

  if (calculation_type & CALC_RELATIONSHIP_COV) {
    retval = calc_cov(pedfile, bed_offset, unfiltered_marker_ct, marker_exclude, marker_ct, mafs, unfiltered_indiv_ct, indiv_exclude, parallel_idx, parallel_tot, outname, outname_end, person_id, max_person_id_len, calculation_type, threads);
    if (retval) {
      goto wdist_ret_2;
    }
  } else if (relationship_or_ibc_req(calculation_type)) {
    indiv_missing_unwt = (unsigned int*)wkspace_alloc(indiv_ct * sizeof(int));
    if (!indiv_missing_unwt) {
      goto wdist_ret_NOMEM;
    }
    fill_int_zero((int*)indiv_missing_unwt, indiv_ct);
;
    triangle_fill(thread_start, indiv_ct, thread_ct, parallel_idx, parallel_tot, 1, 1);
    if (relationship_req(calculation_type)) {
      llxx = thread_start[thread_ct];
      llxx = ((llxx * (llxx - 1)) - (long long)thread_start[0] * (thread_start[0] - 1)) / 2;
      dists_alloc = llxx * sizeof(double);
      if (!(calculation_type & CALC_UNRELATED_HERITABILITY)) {
        // if the memory isn't needed for CALC_UNRELATED_HERITABILITY,
	// positioning the missingness matrix here will let us avoid
	// recalculating it if --distance-matrix or --matrix is requested
        missing_dbl_excluded = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
        if (!missing_dbl_excluded) {
          goto wdist_ret_NOMEM;
        }
        fill_int_zero((int*)missing_dbl_excluded, llxx);
      }
      rel_dists = (double*)wkspace_alloc(dists_alloc);
      if (!rel_dists) {
	goto wdist_ret_NOMEM;
      }
      fill_double_zero(rel_dists, llxx);
    }
    if (calculation_type & CALC_IBC) {
      ii = indiv_ct * 3;
    } else {
      ii = indiv_ct;
    }
    rel_ibc = (double*)wkspace_alloc(ii * sizeof(double));
    if (!rel_ibc) {
      goto wdist_ret_NOMEM;
    }
    fill_double_zero(rel_ibc, ii);
    wkspace_mark = wkspace_base;
    if (relationship_req(calculation_type) && (!missing_dbl_excluded)) {
      missing_dbl_excluded = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
      if (!missing_dbl_excluded) {
	goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)missing_dbl_excluded, llxx);
    }
    if (binary_files) {
      fseeko(pedfile, bed_offset, SEEK_SET);
      ii = 0;
      pp = 0;
      ped_geno = wkspace_alloc(indiv_ct * sizeof(long));
      if (!ped_geno) {
        goto wdist_ret_NOMEM;
      }
      mmasks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
      if (!mmasks) {
	goto wdist_ret_NOMEM;
      }
      gptr = wkspace_alloc(MULTIPLEX_REL * unfiltered_indiv_ct4);
      if (!gptr) {
        goto wdist_ret_NOMEM;
      }
      masks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
      if (!masks) {
	goto wdist_ret_NOMEM;
      }

      // See comments at the beginning of this file, and those in the main
      // CALC_DISTANCE loop.  The main difference between this calculation and
      // the (nonzero exponent) distance calculation is that we have to pad
      // each marker to 3 bits and use + instead of XOR to distinguish the
      // cases.
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
            goto wdist_ret_READ_FAIL;
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
		masks[oo] |= 7LU << (mm * 3);
		mmasks[oo] |= 1LU << (nn + mm);
		indiv_missing_unwt[oo] += 1;
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
		goto wdist_ret_THREAD_CREATE_FAIL;
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
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_missing_thread, (void*)ulii)) {
	      goto wdist_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  incr_dists_rm(missing_dbl_excluded, 0, thread_start);
	  for (oo = 0; oo < thread_ct - 1; oo++) {
	    pthread_join(threads[oo], NULL);
	  }
	}
        printf("\r%d markers complete.", pp);
        fflush(stdout);
      }
      if (relationship_req(calculation_type)) {
        printf("\rRelationship matrix calculation complete.\n");
        dist_ptr = rel_dists;
      } else {
        printf("\n");
      }
      dptr2 = rel_ibc;
      if (calculation_type & CALC_IBC) {
        dptr3 = &(rel_ibc[indiv_ct]);
        dptr4 = &(rel_ibc[indiv_ct * 2]);
      }
      giptr2 = missing_dbl_excluded;
      for (ii = 0; ii < indiv_ct; ii++) {
        uii = marker_ct - indiv_missing_unwt[ii];
	if ((ii >= thread_start[0]) && (ii < thread_start[thread_ct])) {
	  if (relationship_req(calculation_type)) {
	    giptr = indiv_missing_unwt;
	    for (jj = 0; jj < ii; jj++) {
	      *dist_ptr /= uii - (*giptr++) + (*giptr2++);
	      dist_ptr++;
	    }
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
	    if (*dist_ptr++ > rel_cutoff) {
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
              if (*dist_ptr > rel_cutoff) {
                *dist_ptr = 0.0;
                break;
              }
              dist_ptr++;
            }
            if (mm == kk) {
              do {
                mm++;
	        dist_ptr = &(rel_dists[(mm * (mm - 1)) / 2 + kk]);
              } while (*dist_ptr <= rel_cutoff);
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
	    if (*dist_ptr > rel_cutoff) {
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
	    if (*dist_ptr > rel_cutoff) {
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
            giptr = indiv_missing_unwt;
            giptr2 = indiv_missing_unwt;
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
        if (fopen_checked(&outfile, outname, "w")) {
	  goto wdist_ret_OPENFAIL;
	}
        dptr2 = rel_ibc;
        dptr3 = &(rel_ibc[indiv_ct]);
        dptr4 = &(rel_ibc[indiv_ct * 2]);
        if (fprintf(outfile, "FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3\n") < 0) {
          goto wdist_ret_WRITE_FAIL;
        }
	for (ii = 0; ii < indiv_ct; ii++) {
          if (fprintf(outfile, "%d\t%d\t%d\t%g\t%g\t%g\n", ii + 1, ii + 1, marker_ct - indiv_missing_unwt[ii], *dptr3++ - 1.0, *dptr4++ - 1.0, *dptr2++ - 1.0) < 0) {
            goto wdist_ret_WRITE_FAIL;
          }
	}
	if (fclose_null(&outfile)) {
          goto wdist_ret_WRITE_FAIL;
        }
        printf("%s written.\n", outname);
      }
      if (calculation_type & CALC_RELATIONSHIP_MASK) {
        mm = 1;
        nn = calculation_type & CALC_RELATIONSHIP_SHAPEMASK;
	oo = thread_start[thread_ct];
	pp = thread_start[0];
	if (pp == 1) {
	  pp = 0;
	}
        if (calculation_type & CALC_IBC) {
          dptr2 = &(rel_ibc[ibc_type * indiv_ct + pp]);
        } else {
          dptr2 = &(rel_ibc[pp]);
        }
	llxx = ((long long)pp * (pp - 1)) / 2;
	llyy = (((long long)oo * (oo + 1)) / 2) - llxx;
	if (calculation_type & CALC_RELATIONSHIP_BIN) {
          if (nn == CALC_RELATIONSHIP_SQ0) {
            fill_double_zero((double*)ped_geno, indiv_ct - 1);
          }
          strcpy(outname_end, ".rel.bin");
	  if (parallel_tot > 1) {
	    sprintf(&(outname_end[8]), ".%d", parallel_idx + 1);
	  }
          if (fopen_checked(&outfile, outname, "wb")) {
            goto wdist_ret_OPENFAIL;
          }
	  for (ii = pp; ii < oo; ii++) {
	    if (fwrite_checked(&(rel_dists[((long long)ii * (ii - 1)) / 2 - llxx]), ii * sizeof(double), outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (fwrite_checked(dptr2++, sizeof(double), outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
            if (nn == CALC_RELATIONSHIP_TRI) {
              if ((((long long)ii + 1) * (ii + 2) / 2 - llxx) * 100 >= llyy * mm) {
		mm = (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100) / llyy;
                printf("\rWriting... %d%%", mm++);
                fflush(stdout);
              }
            } else {
	      if (nn == CALC_RELATIONSHIP_SQ0) {
		if (fwrite_checked(ped_geno, (indiv_ct - ii - 1) * sizeof(double), outfile)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      } else {
		for (jj = ii + 1; jj < indiv_ct; jj++) {
		  if (fwrite_checked(&(rel_dists[((long long)jj * (jj - 1) / 2) + ii - llxx]), sizeof(double), outfile)) {
		    goto wdist_ret_WRITE_FAIL;
		  }
		}
	      }
	      if ((ii + 1 - pp) * 100 >= mm * (oo - pp)) {
		mm = ((ii + 1 - pp) * 100) / (oo - pp);
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
            }
	  }
          if (fclose_null(&outfile)) {
            goto wdist_ret_WRITE_FAIL;
          }
        } else if (calculation_type & CALC_RELATIONSHIP_GRM) {
          giptr2 = missing_dbl_excluded;
	  if (calculation_type & CALC_RELATIONSHIP_GZ) {
	    if (parallel_tot > 1) {
	      sprintf(outname_end, ".grm.%d.gz", parallel_idx + 1);
	    } else {
	      strcpy(outname_end, ".grm.gz");
	    }
	    if (gzopen_checked(&gz_outfile, outname, "wb")) {
	      goto wdist_ret_OPENFAIL;
	    }
	    dist_ptr = rel_dists;
	    for (ii = pp; ii < oo; ii++) {
	      if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * mm) {
		mm = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	      for (jj = 0; jj < ii; jj += 1) {
		kk = marker_ct - *giptr2++;
		if (!gzprintf(gz_outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dist_ptr++)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
              kk = marker_ct - indiv_missing_unwt[ii];
              if (!gzprintf(gz_outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dptr2++)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    gzclose(gz_outfile);
	    gz_outfile = NULL;
	  } else {
	    strcpy(outname_end, ".grm");
	    if (parallel_tot > 1) {
	      sprintf(&(outname_end[4]), ".%d", parallel_idx + 1);
	    }
            if (fopen_checked(&outfile, outname, "w")) {
	      goto wdist_ret_OPENFAIL;
	    }
	    dist_ptr = rel_dists;
	    for (ii = pp; ii < oo; ii++) {
	      if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * mm) {
		mm = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
		printf("\rWriting... %d%%", mm++);
		fflush(stdout);
	      }
	      for (jj = 0; jj < ii; jj++) {
		kk = marker_ct - *giptr2++;
		if (fprintf(outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dist_ptr++) < 0) {
                  goto wdist_ret_WRITE_FAIL;
                }
	      }
              kk = marker_ct - indiv_missing_unwt[ii];
              if (fprintf(outfile, "%d\t%d\t%d\t%g\n", ii + 1, jj + 1, kk, *dptr2++) < 0) {
                goto wdist_ret_WRITE_FAIL;
              }
	    }
	    if (fclose_null(&outfile)) {
              goto wdist_ret_WRITE_FAIL;
            }
	  }
        } else {
	  if (nn == CALC_RELATIONSHIP_SQ0) {
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
	  if (calculation_type & CALC_RELATIONSHIP_GZ) {
	    if (parallel_tot > 1) {
	      sprintf(outname_end, ".rel.%d.gz", parallel_idx + 1);
	    } else {
	      strcpy(outname_end, ".rel.gz");
	    }
	    if (gzopen_checked(&gz_outfile, outname, "wb")) {
	      goto wdist_ret_OPENFAIL;
	    }
	    dist_ptr = rel_dists;
	    if (pp) {
	      ii = pp;
	    } else {
	      if (!gzprintf(gz_outfile, "%g", *dptr2++)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      if (nn == CALC_RELATIONSHIP_SQ0) {
		if (gzwrite_checked(gz_outfile, ped_geno, (indiv_ct - 1) * 2)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      } else if (nn == CALC_RELATIONSHIP_SQ) {
		for (jj = 1; jj < indiv_ct; jj++) {
		  if (!gzprintf(gz_outfile, "\t%g", rel_dists[(jj * (jj - 1)) / 2])) {
		    goto wdist_ret_WRITE_FAIL;
		  }
		}
	      }
	      if (gzwrite_checked(gz_outfile, "\n", 1)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      ii = 1;
	    }
	    for (; ii < oo; ii++) {
	      if (!gzprintf(gz_outfile, "%g", *dist_ptr++)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = 1; jj < ii; jj++) {
		if (!gzprintf(gz_outfile, "\t%g", *dist_ptr++)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
              if (!gzprintf(gz_outfile, "\t%g", *dptr2++)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      if (nn == CALC_RELATIONSHIP_SQ0) {
		if (gzwrite_checked(gz_outfile, ped_geno, (indiv_ct - jj - 1) * 2)) {
		  goto wdist_ret_WRITE_FAIL;
		}
		if ((ii + 1 - pp) * 100 >= mm * oo) {
		  mm = ((ii + 1 - pp) * 100) / oo;
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
	      } else {
                if (nn == CALC_RELATIONSHIP_SQ) {
                  for (jj = ii + 1; jj < indiv_ct; jj++) {
                    if (!gzprintf(gz_outfile, "\t%g", rel_dists[((jj * (jj - 1)) / 2) + ii])) {
		      goto wdist_ret_WRITE_FAIL;
		    }
                  }
                }
		if (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100 >= llyy * mm) {
		  mm = (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100) / llyy;
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
              }
	      if (gzwrite_checked(gz_outfile, "\n", 1)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    gzclose(gz_outfile);
	    gz_outfile = NULL;
	  } else {
	    strcpy(outname_end, ".rel");
	    if (parallel_tot > 1) {
	      sprintf(&(outname_end[4]), ".%d", parallel_idx + 1);
	    }
            if (fopen_checked(&outfile, outname, "w")) {
	      goto wdist_ret_OPENFAIL;
	    }
	    dist_ptr = rel_dists;
	    if (pp) {
	      ii = pp;
	    } else {
	      if (fprintf(outfile, "%g", *dptr2++) < 0) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      if (nn == CALC_RELATIONSHIP_SQ0) {
		if (fwrite_checked(ped_geno, (indiv_ct - 1) * 2, outfile)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      } else if (nn == CALC_RELATIONSHIP_SQ) {
		for (jj = 1; jj < indiv_ct; jj++) {
		  if (fprintf(outfile, "\t%g", rel_dists[(jj * (jj - 1)) / 2]) < 0) {
		    goto wdist_ret_WRITE_FAIL;
		  }
		}
	      }
	      if (fwrite_checked("\n", 1, outfile)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      ii = 1;
	    }
	    for (; ii < oo; ii++) {
              if (fprintf(outfile, "%g", *dist_ptr++) < 0) {
                goto wdist_ret_WRITE_FAIL;
              }
	      for (jj = 1; jj < ii; jj++) {
		if (fprintf(outfile, "\t%g", *dist_ptr++) < 0) {
                  goto wdist_ret_WRITE_FAIL;
                }
	      }
              if (fprintf(outfile, "\t%g", *dptr2++) < 0) {
                goto wdist_ret_WRITE_FAIL;
              }
	      if (nn == CALC_RELATIONSHIP_SQ0) {
		if (fwrite_checked(ped_geno, (indiv_ct - jj - 1) * 2, outfile)) {
		  goto wdist_ret_WRITE_FAIL;
		}
		if ((ii + 1 - pp) * 100 >= mm * (oo - pp)) {
		  mm = ((ii + 1 - pp) * 100) / (oo - pp);
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
	      } else {
                if (nn == CALC_RELATIONSHIP_SQ) {
                  for (jj = ii + 1; jj < indiv_ct; jj++) {
                    if (fprintf(outfile, "\t%g", rel_dists[((jj * (jj - 1)) / 2) + ii]) < 0) {
                      goto wdist_ret_WRITE_FAIL;
                    }
                  }
                }
		if (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100 >= llyy * mm) {
		  mm = (((long long)(ii + 1) * (ii + 2) / 2 - llxx) * 100) / llyy;
		  printf("\rWriting... %d%%", mm++);
		  fflush(stdout);
		}
              }
	      if (fwrite_checked("\n", 1, outfile)) {
                goto wdist_ret_WRITE_FAIL;
              }
	    }
	    if (fclose_null(&outfile)) {
              goto wdist_ret_WRITE_FAIL;
            }
	  }
	}
	printf("\rRelationship matrix written to %s.\n", outname);
	if (!parallel_idx) {
	  strcpy(&(outname_end[4]), ".id");
	  retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
	  if (retval) {
	    goto wdist_ret_2;
	  }
	}
      }
      wkspace_reset(wkspace_mark);
    }

    if (calculation_type & CALC_REGRESS_REL) {
      retval = regress_rel_main(indiv_exclude, unfiltered_indiv_ct, regress_rel_iters, regress_rel_d, threads);
      if (retval) {
	goto wdist_ret_2;
      }
    }
#ifndef NOLAPACK
    if (calculation_type & CALC_UNRELATED_HERITABILITY) {
      missing_dbl_excluded = NULL;
      ulii = indiv_ct;
      ulii = CACHEALIGN_DBL(ulii * ulii);
      dptr4 = &(rel_dists[ulii]);
      ulii = ulii * 3 + CACHEALIGN_DBL(indiv_ct) * 3;
      wkspace_reset(rel_dists);
      rel_dists = (double*)wkspace_alloc(ulii * sizeof(double));
      if (!rel_dists) {
        goto wdist_ret_NOMEM;
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
      reml_em_one_trait(rel_dists, dptr2, &unrelated_herit_covg, &unrelated_herit_covr, unrelated_herit_tol, calculation_type & CALC_UNRELATED_HERITABILITY_STRICT);
      printf("h^2 estimate: %g\n", unrelated_herit_covg);
    }
#endif
    if ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX | CALC_GENOME)) && (!(calculation_type & CALC_REL_CUTOFF))) {
      wkspace_reset(rel_dists);
    } else {
      wkspace_reset(indiv_missing_unwt);
      indiv_missing_unwt = NULL;
      missing_dbl_excluded = NULL;
    }
  }

  if (distance_req(calculation_type)) {
    triangle_fill(thread_start, indiv_ct, thread_ct, parallel_idx, parallel_tot, 1, 1);
    llxx = thread_start[thread_ct];
    llxx = ((llxx * (llxx - 1)) - (long long)thread_start[0] * (thread_start[0] - 1)) / 2;
    dists_alloc = llxx * sizeof(double);
    // additional + CACHELINE is to fix weird-ass aliasing bug that shows up
    // with -O2 in some cases
    if ((calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) && (!missing_dbl_excluded)) {
      missing_dbl_excluded = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
      if (!missing_dbl_excluded) {
        goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)missing_dbl_excluded, llxx);
      unwt_needed = 1;
      if (!indiv_missing_unwt) {
        indiv_missing_unwt = (unsigned int*)wkspace_alloc(indiv_ct * sizeof(int));
        if (!indiv_missing_unwt) {
          goto wdist_ret_NOMEM;
        }
	unwt_needed_full = 1;
        fill_int_zero((int*)indiv_missing_unwt, indiv_ct);
      }
    }
    dists = (double*)wkspace_alloc(dists_alloc + CACHELINE);
    if (!dists) {
      goto wdist_ret_NOMEM;
    }
    wkspace_mark = wkspace_base;
    if (wt_needed) {
      missing_tot_weights = (unsigned int*)wkspace_alloc(llxx * sizeof(int));
      if (!missing_tot_weights) {
	goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)missing_tot_weights, llxx);
      indiv_missing = (unsigned int*)wkspace_alloc(indiv_ct * sizeof(int));
      if (!indiv_missing) {
	goto wdist_ret_NOMEM;
      }
      fill_int_zero((int*)indiv_missing, indiv_ct);
    }

    if (exp0) {
      idists = (int*)(((char*)wkspace_mark) - CACHEALIGN(llxx * sizeof(int)));
      fill_int_zero(idists, llxx);
      masks = (unsigned long*)wkspace_alloc(indiv_ct * (MULTIPLEX_2DIST / 8));
    } else {
      fill_double_zero(dists, llxx);
      masks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
    }
    if (!masks) {
      goto wdist_ret_NOMEM;
    }
    mmasks = (unsigned long*)wkspace_alloc(indiv_ct * sizeof(long));
    if (!mmasks) {
      goto wdist_ret_NOMEM;
    }

    if (exp0) {
      multiplex = MULTIPLEX_DIST;
      ped_geno = wkspace_alloc(indiv_ct * (MULTIPLEX_2DIST / 8));
    } else {
      multiplex = MULTIPLEX_DIST_EXP;
      ped_geno = wkspace_alloc(indiv_ct * sizeof(long));
    }
    if (!ped_geno) {
      goto wdist_ret_NOMEM;
    }

    if (binary_files) {
      if (pedbuflen < multiplex * unfiltered_indiv_ct4) {
	free(pedbuf);
	pedbuf = (unsigned char*)malloc(multiplex * unfiltered_indiv_ct4 * sizeof(char));
	if (!pedbuf) {
	  goto wdist_ret_NOMEM;
	}
      }
      fseeko(pedfile, bed_offset, SEEK_SET);
      ii = 0; // current SNP index
      pp = 0; // after subtracting out excluded
      while (pp < marker_ct) {
	for (jj = 0; jj < multiplex; jj++) {
	  maf_buf[jj] = 0.5;
	}
	fill_int_zero((int*)wtbuf, multiplex);
	jj = 0; // actual SNPs read

	// For each pair (g_j, g_k) of 2-bit PLINK genotypes, we perform the
	// following operations:
	//
	// 1. XOR each genotype with 01.  This shuffles the genotype
	// representation to:
	//    00 = missing (important for simplifying step 2)
	//    01 = homozygote 1
	//    10 = homozygote 2
	//    11 = heterozygote
	//
	// 2. Next, compute
	//    mask_i := ((g_i | (g_i >> 1)) & 01) * 11
	// which is 00 whenever g_i is missing, and 11 otherwise.
	//
	// 3. Then, (g_j ^ g_k) & (mask_j & mask_k) distinguishes the
	// possible distances between the genotypes:
	//    - It's equal to zero iff either g_j == g_k or one/both is
	//      missing.  It's fine for these cases to overlap since either
	//      way we do not want to increment the numerator of our final
	//      distance.  (We can handle the effect of missingness on the
	//      denominator outside the main loop.)
	//    - It's equal to 01 or 10 iff neither is missing and exactly
	//      one is a heterozygote.
	//    - It's equal to 11 iff one is homozygote rare and the other is
	//      homozygote common.
	//
	// 4. Finally, we perform the update
	//    A_{jk} := A_{jk} + f_0(x_0) + f_1(x_1) + ... + f_31(x_31)
	// in the nonzero exponent case, or
	//    A_{jk} := A_{jk} + f(x_0) + f(x_1) + ... + f(x_959)
	// in the zero exponent case.
	//
	// For nonzero exponents, we structure the update as
	//    A_{jk} := A_{jk} + f_{0-6}(x_{0-6}) + f_{7-13}(x_{7-13}) +
	//              f_{14-19}(x_{14-19}) + f_{20-25}(x_{20-25}) +
	//              f_{26-31}(x_{26-31})
	// which requires 352 KB of table space.  (This is a conservative
	// choice; the 2 MB 8-8-8-8 table would work better on some newer
	// systems.)
	//
	// See the comments at the beginning of this file for discussion of
	// the zero exponent special case.

	while ((jj < multiplex) && (pp < marker_ct)) {
	  while (excluded(marker_exclude, ii)) {
	    ii++;
	    fseeko(pedfile, (off_t)unfiltered_indiv_ct4, SEEK_CUR);
	  }
	  if (fread(&(pedbuf[jj * unfiltered_indiv_ct4]), 1, unfiltered_indiv_ct4, pedfile) < unfiltered_indiv_ct4) {
	    goto wdist_ret_READ_FAIL;
	  }
	  maf_buf[jj] = mafs[ii];
	  if (wt_needed) {
	    wtbuf[jj++] = marker_weights_i[pp++];
	  } else {
	    jj++;
	    pp++;
	  }
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
	    if (wt_needed) {
	      giptr = indiv_missing;
	    }
	    if (unwt_needed_full) {
	      giptr2 = indiv_missing_unwt;
	    }
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
		    ulkk |= 1LU << mm;
		    if (wt_needed) {
		      *giptr += wtbuf[mm + nn];
		    }
		    if (unwt_needed_full) {
		      *giptr2 += 1;
		    }
		  }
		}
		// use xor to convert representation to 0 = missing,
		// 1 or 2 = homozygote, 3 = heterozygote
		ulii ^= FIVEMASK;
		*glptr++ = ulii;
		ulii = (ulii | (ulii >> 1)) & FIVEMASK;
		*glptr2++ = ulii * 3;
		*glptr3 = ulkk;
		ulii = 0;
		ulkk = 0;
		gptr = &(pedbuf[oo / 4 + (nn + BITCT2) * unfiltered_indiv_ct4]);
		for (mm = 0; mm < BITCT2; mm++) {
		  uljj = (gptr[mm * unfiltered_indiv_ct4] >> kk) & 3;
		  ulii |= uljj << (mm * 2);
		  if (uljj == 1) {
		    ulkk |= 1LU << mm;
		    if (wt_needed) {
		      *giptr += wtbuf[mm + nn + BITCT2];
		    }
		    if (unwt_needed_full) {
		      *giptr2 += 1;
		    }
		  }
		}
		ulii ^= FIVEMASK;
		*glptr = ulii;
		ulii = (ulii | (ulii >> 1)) & FIVEMASK;
		*glptr2 = ulii * 3;
		*glptr3++ |= ulkk << BITCT2;
		glptr = &(glptr[(MULTIPLEX_2DIST / BITCT) - 1]);
		glptr2 = &(glptr2[(MULTIPLEX_2DIST / BITCT) - 1]);
		if (wt_needed) {
		  giptr++;
		}
		if (unwt_needed_full) {
		  giptr2++;
		}
	      }
	    }

	    if (wt_needed) {
	      weights_i = &(wtbuf[nn]);
	      for (ulii = 1; ulii < thread_ct; ulii++) {
		if (pthread_create(&(threads[ulii - 1]), NULL, &calc_distm_thread, (void*)ulii)) {
		  goto wdist_ret_THREAD_CREATE_FAIL;
		}
	      }
	      incr_wt_dist_missing(missing_tot_weights, 0);
	      for (oo = 0; oo < thread_ct - 1; oo++) {
		pthread_join(threads[oo], NULL);
	      }
	    }
	    if (unwt_needed) {
	      for (ulii = 1; ulii < thread_ct; ulii++) {
		if (pthread_create(&(threads[ulii - 1]), NULL, &calc_missing_thread, (void*)ulii)) {
		  goto wdist_ret_THREAD_CREATE_FAIL;
		}
	      }
	      incr_dists_rm(missing_dbl_excluded, 0, thread_start);
	      for (oo = 0; oo < thread_ct - 1; oo++) {
		pthread_join(threads[oo], NULL);
	      }
	    }
	  }
	  for (ulii = 1; ulii < thread_ct; ulii++) {
	    if (pthread_create(&(threads[ulii - 1]), NULL, &calc_idist_thread, (void*)ulii)) {
	      goto wdist_ret_THREAD_CREATE_FAIL;
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
		    ulkk |= 1LU << mm;
		    *giptr3 += wtbuf[mm + nn];
		  }
		}
#if __LP64__
		ulii ^= FIVEMASK;
		*glptr++ = ulii;
		ulii = (ulii | (ulii >> 1)) & FIVEMASK;
#else
		// note that FIVEMASK does NOT work here because the mask is
		// only 28 bits
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
		goto wdist_ret_THREAD_CREATE_FAIL;
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
	      goto wdist_ret_THREAD_CREATE_FAIL;
	    }
	  }
	  incr_wt_dist_missing(missing_tot_weights, 0);
	  for (oo = 0; oo < thread_ct - 1; oo++) {
	    pthread_join(threads[oo], NULL);
	  }
	}
	printf("\r%d markers complete.", pp);
	fflush(stdout);
      }
    } else {
      printf("text distance calculation not done (use --make-bed).\n");
      retval = RET_CALC_NOT_YET_SUPPORTED;
      goto wdist_ret_2;
    }
    printf("\rDistance matrix calculation complete.\n");
    if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
      strcpy(outname_end, ".mdist");
      if (fopen_checked(&outfile, outname, "w")) {
	goto wdist_ret_OPENFAIL;
      }
      iptr = idists;
      giptr = missing_dbl_excluded;
      kk = 1;
      for (ii = 0; ii < indiv_ct; ii++) {
	giptr2 = indiv_missing_unwt;
	uii = marker_ct - giptr2[ii];
	for (jj = 0; jj < ii; jj++) {
	  if (fprintf(outfile, "%g ", ((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++)))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("0 ", 2, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	giptr2++;
	for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	  uljj = ulii * (ulii - 1) / 2 + ii;
	  if (fprintf(outfile, "%g ", ((double)idists[uljj]) / (2 * (uii - (*giptr2++) + missing_dbl_excluded[uljj]))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("\n", 1, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	if (ii * 100 >= (kk * indiv_ct)) {
	  kk = (ii * 100) / indiv_ct;
	  printf("\rWriting... %d%%", kk++);
	  fflush(stdout);
	}
      }
      if (fclose_null(&outfile)) {
	goto wdist_ret_WRITE_FAIL;
      }
      printf("\rDistances (proportions) written to %s.\n", outname);
      if (!parallel_idx) {
	strcpy(outname_end, ".mdist.id");
	retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
	if (retval) {
	  goto wdist_ret_2;
	}
      }
    }
    if (calculation_type & CALC_PLINK_IBS_MATRIX) {
      strcpy(outname_end, ".mibs");
      if (fopen_checked(&outfile, outname, "w")) {
	goto wdist_ret_OPENFAIL;
      }
      iptr = idists;
      giptr = missing_dbl_excluded;
      kk = 1;
      for (ii = 0; ii < indiv_ct; ii++) {
	giptr2 = indiv_missing_unwt;
	uii = marker_ct - giptr2[ii];
	for (jj = 0; jj < ii; jj++) {
	  if (fprintf(outfile, "%g ", 1.0 - (((double)(*iptr++)) / (2 * (uii - (*giptr2++) + (*giptr++))))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("1 ", 2, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	giptr2++;
	for (ulii = ii + 1; ulii < indiv_ct; ulii++) {
	  uljj = (ulii * (ulii - 1)) / 2 + ii;
	  if (fprintf(outfile, "%g ", 1.0 - (((double)idists[uljj]) / (2 * (uii - (*giptr2++) + missing_dbl_excluded[uljj])))) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fwrite_checked("\n", 1, outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	if (ii * 100 >= (kk * indiv_ct)) {
	  kk = (ii * 100) / indiv_ct;
	  printf("\rWriting... %d%%", kk++);
	  fflush(stdout);
	}
      }
      fclose_null(&outfile);
      printf("\rIBS matrix written to %s.\n", outname);
      strcpy(outname_end, ".mibs.id");
      retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
      if (retval) {
	goto wdist_ret_2;
      }
    }
    if (wt_needed) {
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
    }
  }

  if (calculation_type & CALC_DISTANCE_MASK) {
    if (!parallel_idx) {
      if (calculation_type & CALC_DISTANCE_SNPS) {
	strcpy(outname_end, ".dist.id");
	retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
	if (retval) {
	  goto wdist_ret_2;
	}
      }
      if (calculation_type & CALC_DISTANCE_IBS) {
	strcpy(outname_end, ".mibs.id");
	retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
	if (retval) {
	  goto wdist_ret_2;
	}
      }
      if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
	strcpy(outname_end, ".mdist.id");
	retval = write_ids(outname, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len);
	if (retval) {
	  goto wdist_ret_2;
	}
      }
    }
    mm = calculation_type & CALC_DISTANCE_SHAPEMASK;
    oo = thread_start[thread_ct];
    pp = thread_start[0];
    nn = calculation_type & CALC_DISTANCE_SNPS;
    qq = calculation_type & CALC_DISTANCE_IBS;
    rr = calculation_type & CALC_DISTANCE_1_MINUS_IBS;
    dyy = 0.5 / (double)marker_ct;
    if (pp == 1) {
      pp = 0;
    }
    llxx = ((long long)pp * (pp - 1)) / 2;
    llyy = (((long long)oo * (oo - 1)) / 2) - llxx;
    if (mm == CALC_DISTANCE_SQ0) {
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
    if (calculation_type & CALC_DISTANCE_GZ) {
      if (distance_open_gz(&gz_outfile, &gz_outfile2, &gz_outfile3, outname, outname_end, calculation_type, parallel_idx, parallel_tot)) {
        goto wdist_ret_OPENFAIL;
      }
      if (pp) {
	ii = pp;
      } else {
	if (mm == CALC_DISTANCE_SQ0) {
	  if (nn) {
	    if (gzwrite_checked(gz_outfile, &(ped_geno[1]), indiv_ct * 2 - 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (gzwrite_checked(gz_outfile, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (qq) {
	    ped_geno[1] = '1';
	    if (gzwrite_checked(gz_outfile2, &(ped_geno[1]), indiv_ct * 2 - 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
            if (gzwrite_checked(gz_outfile2, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    ped_geno[1] = '0';
	  }
	  if (rr) {
	    if (gzwrite_checked(gz_outfile3, &(ped_geno[1]), indiv_ct * 2 - 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
            if (gzwrite_checked(gz_outfile3, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	} else if (mm == CALC_DISTANCE_SQ) {
	  if (nn) {
	    if (gzwrite_checked(gz_outfile, "0", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < indiv_ct; jj++) {
	      if (!gzprintf(gz_outfile, "\t%g", dists[(jj * (jj - 1)) / 2])) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (gzwrite_checked(gz_outfile, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (qq) {
	    if (gzwrite_checked(gz_outfile2, "1", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < indiv_ct; jj++) {
	      if (!gzprintf(gz_outfile2, "\t%g", 1.0 - dists[(jj * (jj - 1)) / 2] * dyy)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (gzwrite_checked(gz_outfile2, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (rr) {
	    if (gzwrite_checked(gz_outfile3, "0", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < indiv_ct; jj++) {
	      if (!gzprintf(gz_outfile3, "\t%g", dists[(jj * (jj - 1)) / 2] * dyy)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (gzwrite_checked(gz_outfile3, "\n", 1)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	}
	ii = 1;
      }
      if (nn) {
        dist_ptr = dists;
	for (; ii < oo; ii++) {
	  if (!gzprintf(gz_outfile, "%g", *dist_ptr++)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  for (jj = 1; jj < ii; jj++) {
	    if (!gzprintf(gz_outfile, "\t%g", *dist_ptr++)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (gzwrite_checked(gz_outfile, ped_geno, (indiv_ct - jj) * 2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if ((ii - pp) * 100 >= kk * (oo - pp)) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  } else {
	    if (mm == CALC_DISTANCE_SQ) {
	      if (gzwrite_checked(gz_outfile, "\t0", 2)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		if (!gzprintf(gz_outfile, "\t%g", dists[((jj * (jj - 1)) / 2) + ii])) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * kk) {
	      kk = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (gzwrite_checked(gz_outfile, "\n", 1)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	gzclose(gz_outfile);
	distance_print_done(0, outname, outname_end);
	kk = 1;
	gz_outfile = NULL;
      }
      if (pp) {
	ii = pp;
      } else {
	ii = 1;
      }
      if (qq) {
        dist_ptr = dists;
	ped_geno[1] = '1';
	for (; ii < oo; ii++) {
	  if (!gzprintf(gz_outfile2, "%g", 1.0 - (*dist_ptr++) * dyy)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  for (jj = 1; jj < ii; jj++) {
	    if (!gzprintf(gz_outfile2, "\t%g", 1.0 - (*dist_ptr++) * dyy)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (gzwrite_checked(gz_outfile2, ped_geno, (indiv_ct - jj) * 2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if ((ii - pp) * 100 >= kk * (oo - pp)) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  } else {
	    if (mm == CALC_DISTANCE_SQ) {
	      if (gzwrite_checked(gz_outfile2, "\t1", 2)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		if (!gzprintf(gz_outfile2, "\t%g", 1.0 - dists[((jj * (jj - 1)) / 2) + ii] * dyy)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * kk) {
	      kk = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (gzwrite_checked(gz_outfile2, "\n", 1)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	ped_geno[1] = '0';
	gzclose(gz_outfile2);
	distance_print_done(1, outname, outname_end);
	kk = 1;
	gz_outfile2 = NULL;
      }
      if (pp) {
	ii = pp;
      } else {
	ii = 1;
      }
      if (rr) {
        dist_ptr = dists;
	for (; ii < oo; ii++) {
	  if (!gzprintf(gz_outfile3, "%g", (*dist_ptr++) * dyy)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  for (jj = 1; jj < ii; jj++) {
	    if (!gzprintf(gz_outfile3, "\t%g", (*dist_ptr++) * dyy)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (gzwrite_checked(gz_outfile3, ped_geno, (indiv_ct - jj) * 2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if ((ii - pp) * 100 >= kk * (oo - pp)) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  } else {
	    if (mm == CALC_DISTANCE_SQ) {
	      if (gzwrite_checked(gz_outfile3, "\t0", 2)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		if (!gzprintf(gz_outfile3, "\t%g", dists[((jj * (jj - 1)) / 2) + ii] * dyy)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * kk) {
	      kk = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (gzwrite_checked(gz_outfile3, "\n", 1)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	gzclose(gz_outfile3);
	distance_print_done(2, outname, outname_end);
	gz_outfile3 = NULL;
      }
    } else if (calculation_type & CALC_DISTANCE_BIN) {
      if (distance_open(&outfile, &outfile2, &outfile3, outname, outname_end, ".bin", "wb", calculation_type, parallel_idx, parallel_tot)) {
	goto wdist_ret_OPENFAIL;
      }
      if (mm == CALC_DISTANCE_TRI) {
	if (nn) {
          if (fwrite_checked(dists, llyy * sizeof(double), outfile)) {
            goto wdist_ret_WRITE_FAIL;
          }
	}
	if (qq) {
	  dist_ptr = dists;
          for (ii = 0; ii < llyy; ii++) {
	    dxx = 1.0 - (*dist_ptr++) * dyy;
	    if (fwrite_checked(&dxx, sizeof(double), outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	}
	if (rr) {
	  dist_ptr = dists;
          for (ii = 0; ii < llyy; ii++) {
	    dxx = (*dist_ptr++) * dyy;
	    if (fwrite_checked(&dxx, sizeof(double), outfile3)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	}
      } else {
	if (mm == CALC_DISTANCE_SQ0) {
	  fill_double_zero((double*)ped_geno, indiv_ct);
	}
	if (nn) {
	  dxx = 0.0;
	  dist_ptr = dists;
	  for (ii = pp; ii < oo; ii++) {
	    if (fwrite_checked(dist_ptr, ii * sizeof(double), outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    dist_ptr = &(dist_ptr[ii]);
	    if (mm == CALC_DISTANCE_SQ0) {
	      if (fwrite_checked(ped_geno, (indiv_ct - ii) * sizeof(double), outfile)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    } else {
	      // square matrix, no need to handle parallel case
	      if (fwrite_checked(&dxx, sizeof(double), outfile)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		if (fwrite_checked(&(dists[(jj * (jj - 1)) / 2 + ii]), sizeof(double), outfile)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if ((ii - pp) * 100 >= kk * (oo - pp)) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (fclose_null(&outfile)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  distance_print_done(0, outname, outname_end);
	  kk = 1;
	}
	if (qq) {
	  dist_ptr = dists;
	  dzz = 1.0;
	  ped_geno[1] = '1';
	  for (ii = pp; ii < oo; ii++) {
	    for (jj = 0; jj < ii; jj++) {
	      dxx = 1.0 - (*dist_ptr++) * dyy;
	      if (fwrite_checked(&dxx, sizeof(double), outfile2)) {
	        goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (mm == CALC_DISTANCE_SQ0) {
	      if (fwrite_checked(ped_geno, (indiv_ct - ii) * sizeof(double), outfile2)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    } else {
	      // square matrix
	      if (fwrite_checked(&dzz, sizeof(double), outfile2)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		dxx = 1.0 - dists[(jj * (jj - 1)) / 2 + ii] * dyy;
		if (fwrite_checked(&dxx, sizeof(double), outfile2)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if ((ii - pp) * 100 >= kk * (oo - pp)) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  ped_geno[1] = '0';
	  if (fclose_null(&outfile2)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  distance_print_done(1, outname, outname_end);
	  kk = 1;
	}
	if (rr) {
	  dist_ptr = dists;
	  dzz = 0.0;
	  for (ii = pp; ii < oo; ii++) {
	    for (jj = 0; jj < ii; jj++) {
	      dxx = (*dist_ptr++) * dyy;
	      if (fwrite_checked(&dxx, sizeof(double), outfile3)) {
	        goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (mm == CALC_DISTANCE_SQ0) {
	      if (fwrite_checked(ped_geno, (indiv_ct - ii) * sizeof(double), outfile3)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    } else {
	      // square matrix
	      if (fwrite_checked(&dzz, sizeof(double), outfile3)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		dxx = dists[(jj * (jj - 1)) / 2 + ii] * dyy;
		if (fwrite_checked(&dxx, sizeof(double), outfile3)) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if ((ii - pp) * 100 >= kk * (oo - pp)) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (fclose_null(&outfile3)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  distance_print_done(2, outname, outname_end);
	}
      }
    } else {
      if (distance_open(&outfile, &outfile2, &outfile3, outname, outname_end, "", "w", calculation_type, parallel_idx, parallel_tot)) {
	goto wdist_ret_OPENFAIL;
      }
      if (nn) {
	if (pp) {
	  ii = pp;
	} else {
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (fwrite_checked(&(ped_geno[1]), indiv_ct * 2 - 1, outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (fwrite_checked("\n", 1, outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  } else if (mm == CALC_DISTANCE_SQ) {
	    if (fwrite_checked("0", 1, outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < indiv_ct; jj++) {
	      if (fprintf(outfile, "\t%g", dists[(jj * (jj - 1)) / 2]) < 0) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (fwrite_checked("\n", 1, outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  ii = 1;
	}
        dist_ptr = dists;
	for (; ii < oo; ii++) {
	  if (fprintf(outfile, "%g", *dist_ptr++) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  for (jj = 1; jj < ii; jj++) {
	    if (fprintf(outfile, "\t%g", *dist_ptr++) < 0) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (fwrite_checked(ped_geno, (indiv_ct - jj) * 2, outfile)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if ((ii - pp) * 100 >= (kk * (oo - pp))) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  } else {
	    if (mm == CALC_DISTANCE_SQ) {
	      if (fwrite_checked("\t0", 2, outfile)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		if (fprintf(outfile, "\t%g", dists[((jj * (jj - 1)) / 2) + ii]) < 0) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * kk) {
	      kk = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (fwrite_checked("\n", 1, outfile)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	distance_print_done(0, outname, outname_end);
	kk = 1;
      }
      if (qq) {
	ped_geno[1] = '1';
	if (pp) {
	  ii = pp;
	} else {
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (fwrite_checked(&(ped_geno[1]), indiv_ct * 2 - 1, outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (fwrite_checked("\n", 1, outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  } else if (mm == CALC_DISTANCE_SQ) {
	    if (fwrite_checked("1", 1, outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < indiv_ct; jj++) {
	      if (fprintf(outfile2, "\t%g", 1.0 - dists[(jj * (jj - 1)) / 2] * dyy) < 0) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (fwrite_checked("\n", 1, outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  } else {
	    if (fwrite_checked("1\n", 2, outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  ii = 1;
	}
        dist_ptr = dists;
	for (; ii < oo; ii++) {
	  if (fprintf(outfile2, "%g", 1.0 - (*dist_ptr++) * dyy) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  for (jj = 1; jj < ii; jj++) {
	    if (fprintf(outfile2, "\t%g", 1.0 - (*dist_ptr++) * dyy) < 0) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (fwrite_checked(ped_geno, (indiv_ct - jj) * 2, outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if ((ii - pp) * 100 >= (kk * (oo - pp))) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  } else {
	    if (fwrite_checked("\t1", 2, outfile2)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (mm == CALC_DISTANCE_SQ) {
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		if (fprintf(outfile2, "\t%g", 1.0 - dists[((jj * (jj - 1)) / 2) + ii] * dyy) < 0) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * kk) {
	      kk = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (fwrite_checked("\n", 1, outfile2)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	ped_geno[1] = '0';
	if (fclose_null(&outfile2)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	distance_print_done(1, outname, outname_end);
	kk = 1;
      }
      if (rr) {
	if (pp) {
	  ii = pp;
	} else {
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (fwrite_checked(&(ped_geno[1]), indiv_ct * 2 - 1, outfile3)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if (fwrite_checked("\n", 1, outfile3)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  } else if (mm == CALC_DISTANCE_SQ) {
	    if (fwrite_checked("0", 1, outfile3)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    for (jj = 1; jj < indiv_ct; jj++) {
	      if (fprintf(outfile3, "\t%g", dists[(jj * (jj - 1)) / 2] * dyy) < 0) {
		goto wdist_ret_WRITE_FAIL;
	      }
	    }
	    if (fwrite_checked("\n", 1, outfile3)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  ii = 1;
	}
        dist_ptr = dists;
	for (; ii < oo; ii++) {
	  if (fprintf(outfile3, "%g", (*dist_ptr++) * dyy) < 0) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	  for (jj = 1; jj < ii; jj++) {
	    if (fprintf(outfile3, "\t%g", (*dist_ptr++) * dyy) < 0) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	  }
	  if (mm == CALC_DISTANCE_SQ0) {
	    if (fwrite_checked(ped_geno, (indiv_ct - jj) * 2, outfile3)) {
	      goto wdist_ret_WRITE_FAIL;
	    }
	    if ((ii - pp) * 100 >= (kk * (oo - pp))) {
	      kk = ((ii - pp) * 100) / (oo - pp);
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  } else {
	    if (mm == CALC_DISTANCE_SQ) {
	      if (fwrite_checked("\t0", 2, outfile3)) {
		goto wdist_ret_WRITE_FAIL;
	      }
	      for (jj = ii + 1; jj < indiv_ct; jj++) {
		if (fprintf(outfile3, "\t%g", dists[((jj * (jj - 1)) / 2) + ii] * dyy) < 0) {
		  goto wdist_ret_WRITE_FAIL;
		}
	      }
	    }
	    if (((long long)ii * (ii + 1) / 2 - llxx) * 100 >= llyy * kk) {
	      kk = (((long long)ii * (ii + 1) / 2 - llxx) * 100) / llyy;
	      printf("\rWriting... %d%%", kk++);
	      fflush(stdout);
	    }
	  }
	  if (fwrite_checked("\n", 1, outfile3)) {
	    goto wdist_ret_WRITE_FAIL;
	  }
	}
	if (fclose_null(&outfile3)) {
	  goto wdist_ret_WRITE_FAIL;
	}
	distance_print_done(2, outname, outname_end);
      }
    }
    wkspace_reset(wkspace_mark);
  } else if (calculation_type & CALC_LOAD_DISTANCES) {
    dists_alloc = ((long)indiv_ct * (indiv_ct - 1)) * (sizeof(double) / 2);
    dists = (double*)wkspace_alloc(dists_alloc);
    if (!dists) {
      goto wdist_ret_NOMEM;
    }
    if (fopen_checked(&loaddistfile, loaddistname, "rb")) {
      goto wdist_ret_OPENFAIL;
    }
    if (fseeko(loaddistfile, 0, SEEK_END)) {
      goto wdist_ret_READ_FAIL;
    }
    if (ftello(loaddistfile) != dists_alloc) {
      printf("Invalid --load-dists file.  (Triangular binary of size %lld expected.)\n", dists_alloc);
      goto wdist_ret_INVALID_FORMAT;
    }
    rewind(loaddistfile);
    if (fread(dists, 1, dists_alloc, loaddistfile) < dists_alloc) {
      goto wdist_ret_READ_FAIL;
    }
    fclose(loaddistfile);
  }

  if (calculation_type & CALC_GROUPDIST) {
    collapse_arr(pheno_c, 1, indiv_exclude, unfiltered_indiv_ct);
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
    ped_geno = wkspace_alloc(thread_ct * CACHEALIGN(high_ct + low_ct + (jackknife_d + 1) * sizeof(int)));
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
      // this can be sped up using the same method used in regress-distance,
      // if it's important
      jackknife_iters = (groupdist_iters + thread_ct - 1) / thread_ct;

      if (groupdist_d) {
	jackknife_d = groupdist_d;
      } else {
	jackknife_d = set_default_jackknife_d(high_ct + low_ct);
      }

      for (ulii = 1; ulii < thread_ct; ulii++) {
	if (pthread_create(&(threads[ulii - 1]), NULL, &groupdist_jack_thread, (void*)ulii)) {
	  goto wdist_ret_THREAD_CREATE_FAIL;
	}
      }
      ulii = 0;
      groupdist_jack_thread((void*)ulii);
      for (ii = 1; ii < thread_ct; ii++) {
        pthread_join(threads[ii - 1], NULL);
        for (jj = 0; jj < 9; jj++) {
	  calc_result[jj][0] += calc_result[jj][ii];
        }
      }
      dxx = 1.0 / thread_ct;
      calc_result[0][0] *= dxx;
      calc_result[1][0] *= dxx;
      calc_result[2][0] *= dxx;
      dxx /= (jackknife_iters - 1) * thread_ct;
      for (ii = 3; ii < 9; ii++) {
        calc_result[ii][0] *= dxx;
      }
      printf("\r  AA mean - AU mean avg difference (s.e.): %g (%g)\n", calc_result[0][0] - calc_result[1][0], sqrt(((high_ct + low_ct) / ((double)jackknife_d)) * (calc_result[3][0] + calc_result[4][0] - 2 * calc_result[6][0])));
      printf("  AA mean - UU mean avg difference (s.e.): %g (%g)\n", calc_result[0][0] - calc_result[2][0], sqrt(((high_ct + low_ct) / ((double)jackknife_d)) * (calc_result[3][0] + calc_result[5][0] - 2 * calc_result[7][0])));
      printf("  AU mean - UU mean avg difference (s.e.): %g (%g)\n", calc_result[1][0] - calc_result[2][0], sqrt(((high_ct + low_ct) / ((double)jackknife_d)) * (calc_result[4][0] + calc_result[5][0] - 2 * calc_result[8][0])));
    }
  }
  if (calculation_type & CALC_REGRESS_DISTANCE) {
    // beta = (mean(xy) - mean(x)*mean(y)) / (mean(x^2) - mean(x)^2)
    if (unfiltered_indiv_ct != indiv_ct) {
      collapse_arr((char*)pheno_d, sizeof(double), indiv_exclude, unfiltered_indiv_ct);
    }
    if (!(calculation_type & CALC_REGRESS_REL)) {
      print_pheno_stdev(pheno_d);
    }
    ulii = indiv_ct;
    ulii = ulii * (ulii - 1) / 2;
    reg_tot_xy = 0.0;
    reg_tot_x = 0.0;
    reg_tot_y = 0.0;
    reg_tot_xx = 0.0;
    reg_tot_yy = 0.0;
    dptr4 = dists;
    dist_ptr = pheno_d;
    // Linear regression slope is a function of sum(xy), sum(x), sum(y),
    // sum(x^2), and n.  To speed up the jackknife calculation, we precompute
    // (i) the global xy, x, y, x^2, and
    // (ii) the xy, x, y, x^2 for each row.
    // Then for each delete-d jackknife iteration, we take the global sums,
    // subtract the partial row sums corresponding to the deleted individuals,
    // and then add back the elements in the intersection of two deletions.
    jackknife_precomp = (double*)wkspace_alloc(indiv_ct * JACKKNIFE_VALS_DIST * sizeof(double));
    if (!jackknife_precomp) {
      goto wdist_ret_NOMEM;
    }
    fill_double_zero(jackknife_precomp, indiv_ct * JACKKNIFE_VALS_DIST);
    for (ii = 1; ii < indiv_ct; ii++) {
      dzz = *(++dist_ptr);
      dptr2 = pheno_d;
      dptr3 = &(jackknife_precomp[ii * JACKKNIFE_VALS_DIST]);
      dptr5 = jackknife_precomp;
      while (dptr2 < dist_ptr) {
	dxx = (dzz + *dptr2++) * 0.5;
	dyy = (*dptr4++);
	dww = dxx * dyy;
	dvv = dxx * dxx;
	duu = dyy * dyy;
	reg_tot_xy += dww;
	*dptr3 += dww;
	*dptr5 += dww;
	dptr5++;
	reg_tot_x += dxx;
	dptr3[1] += dxx;
	*dptr5 += dxx;
	dptr5++;
	reg_tot_y += dyy;
	dptr3[2] += dyy;
	*dptr5 += dyy;
	dptr5++;
	reg_tot_xx += dvv;
	dptr3[3] += dvv;
	*dptr5 += dvv;
	dptr5++;
	reg_tot_yy += duu;
	dptr3[4] += duu;
	*dptr5 += duu;
	dptr5++;
      }
    }

    dxx = ulii;
    printf("Regression slope (y = genomic distance, x = avg phenotype): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y / dxx) / (reg_tot_xx - reg_tot_x * reg_tot_x / dxx));
    printf("Regression slope (y = avg phenotype, x = genomic distance): %g\n", (reg_tot_xy - reg_tot_x * reg_tot_y / dxx) / (reg_tot_yy - reg_tot_y * reg_tot_y / dxx));

    jackknife_iters = (regress_iters + thread_ct - 1) / thread_ct;
    if (regress_d) {
      jackknife_d = regress_d;
    } else {
      jackknife_d = set_default_jackknife_d(indiv_ct);
    }
    ped_geno = wkspace_alloc(thread_ct * CACHEALIGN(indiv_ct + (jackknife_d + 1) * sizeof(int)));
    if (!ped_geno) {
      goto wdist_ret_NOMEM;
    }
    for (ulii = 1; ulii < thread_ct; ulii++) {
      if (pthread_create(&(threads[ulii - 1]), NULL, &regress_jack_thread, (void*)ulii)) {
	goto wdist_ret_THREAD_CREATE_FAIL;
      }
    }
    ulii = 0;
    regress_jack_thread((void*)ulii);
    dyy = calc_result[0][0]; // sum
    dzz = calc_result[1][0]; // sum of squares
    dww = calc_result[2][0]; // reverse regression sum
    dvv = calc_result[3][0]; // reverse regression sum of squares
    for (ii = 0; ii < thread_ct - 1; ii++) {
      pthread_join(threads[ii], NULL);
      dyy += calc_result[0][ii + 1];
      dzz += calc_result[1][ii + 1];
      dww += calc_result[2][ii + 1];
      dvv += calc_result[3][ii + 1];
    }
    regress_iters = jackknife_iters * thread_ct;
    printf("\rJackknife s.e.: %g\n", sqrt((indiv_ct / ((double)jackknife_d)) * (dzz - dyy * dyy / regress_iters) / (regress_iters - 1)));
    printf("Jackknife s.e. (y = avg phenotype): %g\n", sqrt((indiv_ct / ((double)jackknife_d)) * (dvv - dww * dww / regress_iters) / (regress_iters - 1)));
  }

  if (calculation_type & CALC_GENOME) {
    wkspace_reset((char*)CACHEALIGN((unsigned long)(&(marker_pos[unfiltered_marker_ct]))));
    retval = calc_genome(threads, pedfile, bed_offset, marker_ct, unfiltered_marker_ct, marker_exclude, mafs, marker_pos_start, unfiltered_indiv_ct, indiv_exclude, person_id, max_person_id_len, paternal_id, max_paternal_id_len, maternal_id, max_maternal_id_len, founder_info, parallel_idx, parallel_tot, outname, outname_end, nonfounders, calculation_type, genome_output_gz, genome_output_full, genome_ibd_unbounded);
    if (retval) {
      goto wdist_ret_2;
    }
  }

  while (0) {
  wdist_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  wdist_ret_OPENFAIL:
    retval = RET_OPENFAIL;
    break;
  wdist_ret_INVALID_FORMAT:
    retval = RET_INVALID_FORMAT;
    break;
  wdist_ret_INVALID_CMDLINE:
    retval = RET_INVALID_CMDLINE;
    break;
  wdist_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  wdist_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  wdist_ret_THREAD_CREATE_FAIL:
    printf("Error: Could not create thread.\n");
    retval = RET_THREAD_CREATE_FAIL;
    for (uljj = 0; uljj < ulii - 1; uljj++) {
      pthread_join(threads[uljj], NULL);
    }
    break;
  }
 wdist_ret_2:
  gzclose_cond(gz_outfile3);
  gzclose_cond(gz_outfile2);
  gzclose_cond(gz_outfile);
  free_cond(person_missing_cts);
  free_cond(marker_alleles_tmp);
  fclose_cond(bedtmpfile);
  fclose_cond(bimtmpfile);
  fclose_cond(famtmpfile);
 wdist_ret_1:
  free_cond(het_probs);
  free_cond(founder_info);
  free_cond(maternal_id);
  free_cond(paternal_id);
  free_cond(person_id);
  free_cond(indiv_exclude);
  free_cond(pheno_d);
  free_cond(pheno_c);
  free_cond(phenor_d);
  free_cond(phenor_c);
  if (!binary_files) {
    free_cond(marker_allele_cts);
    free_cond(marker_alleles);
  }
  free_cond(line_locs);
  free_cond(id_buf);
  free_cond(id_list);
  free_cond(pid_list);
  free_cond(pedbuf);
  fclose_cond(loaddistfile);
  fclose_cond(freqfile);
  fclose_cond(filterfile);
  fclose_cond(phenofile);
  fclose_cond(famfile);
  fclose_cond(mapfile);
  fclose_cond(pedfile);
  fclose_cond(outfile3);
  fclose_cond(outfile2);
  fclose_cond(outfile);
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
  char outname[FNAMESIZE];
  char mapname[FNAMESIZE];
  char pedname[FNAMESIZE];
  char famname[FNAMESIZE];
  char phenoname[FNAMESIZE];
  char extractname[FNAMESIZE];
  char excludename[FNAMESIZE];
  char keepname[FNAMESIZE];
  char removename[FNAMESIZE];
  char filtername[FNAMESIZE];
  char freqname[FNAMESIZE];
  char loaddistname[FNAMESIZE];
  char* makepheno_str = NULL;
  char* filterval = NULL;
  char* argptr;
  char* sptr;
  int retval;
  int load_params = 0; // describes what file parameters have been provided
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
  double rel_cutoff = 0.025;
  int cur_arg = 1;
  int calculation_type = 0;
  char* bubble;
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
  int regress_rel_iters = ITERS_DEFAULT;
  int regress_rel_d = 0;
  double unrelated_herit_tol = 0.0000001;
  double unrelated_herit_covg = 0.45;
  double unrelated_herit_covr = 0.55;
  int ibc_type = 0;
  int parallel_idx = 0;
  int parallel_tot = 1;
  int nonfounders = 0;
  int ii;
  int jj;
  int kk;
  unsigned int chrom_mask = DEFAULT_CHROM_MASK;
  int ppc_gap = DEFAULT_PPC_GAP;
  unsigned long int rseed = 0;
  int genome_output_gz = 0;
  int genome_output_full = 0;
  int genome_ibd_unbounded = 0;
  FILE* scriptfile;
  int num_params;
  int in_param;
  printf(ver_str);
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
  keepname[0] = '\0';
  removename[0] = '\0';
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
      if (fopen_checked(&scriptfile, argv[cur_arg + 1], "rb")) {
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
      printf("Note: --make-bed is effectively always on.\n");
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
      if (filtername[0]) {
        printf("Error: --keep and --filter flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (keepname[0]) {
        printf("Error: Duplicate --keep flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii == 0) {
        printf("Error: Missing --keep parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 1) {
        printf("Error: Too many --keep parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--keep filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(keepname, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--remove")) {
      if (removename[0]) {
        printf("Error: Duplicate --remove flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii == 0) {
        printf("Error: Missing --remove parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 1) {
        printf("Error: Too many --remove parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("--remove filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(removename, argv[cur_arg + 1]);
      cur_arg += 2;
    } else if (!strcmp(argptr, "--filter")) {
      if (keepname[0]) {
        printf("Error: --keep and --filter flags cannot coexist.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (filtername[0]) {
        printf("Error: Duplicate --filter flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii < 2) {
        printf("Error: Not enough --filter parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 2) {
        printf("Error: Too many --filter parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
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
#ifndef __LP64__
      if (malloc_size_mb > 2944) {
	printf("Error: --memory parameter too large for 32-bit version.\n");
	return dispmsg(RET_INVALID_CMDLINE);
      }
#endif
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
      if (calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) {
        printf("Error: --exponent flag cannot be used with --distance-matrix or --matrix.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
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
    } else if (!strcmp(argptr, "--nonfounders")) {
      nonfounders = 1;
      cur_arg += 1;
    } else if ((!strcmp(argptr, "--rel-cutoff")) || (!strcmp(argptr, "--grm-cutoff"))) {
      if (parallel_tot > 1) {
	printf("Error: --parallel cannot be used with --rel-cutoff.  (Use a combination of\n--make-rel, --keep/--remove, and a filtering script.)%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
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
	if (sscanf(argv[cur_arg + 1], "%lg", &rel_cutoff) != 1) {
	  printf("Error: Invalid --rel-cutoff parameter.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
	if ((rel_cutoff <= 0.0) || (rel_cutoff >= 1.0)) {
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
      if (strlen(argv[cur_arg + 1]) > (FNAMESIZE - MAX_POST_EXT)) {
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
      if (ii > 3) {
        printf("Error: --make-rel accepts at most 3 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      jj = 0;
      for (kk = 1; kk <= ii; kk++) {
        if (!strcmp(argv[cur_arg + kk], "cov")) {
          if (calculation_type & CALC_IBC) {
            printf("Error: --make-rel 'cov' modifier cannot coexist with --ibc flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (calculation_type & CALC_UNRELATED_HERITABILITY) {
            printf("Error: --make-rel 'cov' modifier cannot coexist with\n--unrelated-heritability flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_COV;
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
          if ((calculation_type & CALC_RELATIONSHIP_SHAPEMASK) == CALC_RELATIONSHIP_SQ0) {
            printf("Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if ((calculation_type & CALC_RELATIONSHIP_SHAPEMASK) == CALC_RELATIONSHIP_TRI) {
            printf("Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if (parallel_tot > 1) {
            printf("Error: --parallel cannot be used with '--make-rel square'.  Use '--make-rel\nsquare0' or plain '--make-rel' instead.%s", errstr_append);
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
          calculation_type |= CALC_RELATIONSHIP_SQ;
        } else if (!strcmp(argv[cur_arg + kk], "square0")) {
          if ((calculation_type & CALC_RELATIONSHIP_SHAPEMASK) == CALC_RELATIONSHIP_SQ) {
            printf("Error: --make-rel 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if ((calculation_type & CALC_RELATIONSHIP_SHAPEMASK) == CALC_RELATIONSHIP_TRI) {
            printf("Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_SQ0;
	} else if (!strcmp(argv[cur_arg + kk], "triangle")) {
          if ((calculation_type & CALC_RELATIONSHIP_SHAPEMASK) == CALC_RELATIONSHIP_SQ) {
            printf("Error: --make-rel 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if ((calculation_type & CALC_RELATIONSHIP_SHAPEMASK) == CALC_RELATIONSHIP_SQ0) {
            printf("Error: --make-rel 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_TRI;
        } else if ((!strcmp(argv[cur_arg + kk], "ibc1")) || (!strcmp(argv[cur_arg + kk], "ibc2"))) {
          if (calculation_type & CALC_RELATIONSHIP_COV) {
            printf("Error: --make-rel 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
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
      if (((calculation_type & CALC_RELATIONSHIP_BIN) && (!(calculation_type & CALC_RELATIONSHIP_SHAPEMASK))) && (parallel_tot > 1)) {
        printf("Error: --parallel cannot be used with plain '--make-rel bin'.  Use '--make-rel\nbin square0' or '--make-rel bin triangle' instead.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += ii + 1;
      if (!(calculation_type & CALC_RELATIONSHIP_SHAPEMASK)) {
        calculation_type |= CALC_RELATIONSHIP_TRI;
      }
    } else if (!strcmp(argptr, "--make-grm")) {
      if (calculation_type & CALC_RELATIONSHIP_MASK) {
        printf("Error: --make-grm cannot coexist with another relationship matrix file\ncreation flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 2) {
        printf("Error: --make-grm accepts at most 2 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      calculation_type |= CALC_RELATIONSHIP_GZ;
      for (kk = 1; kk <= ii; kk++) {
        if (!strcmp(argv[cur_arg + kk], "cov")) {
          if (calculation_type & CALC_IBC) {
            printf("Error: --make-grm 'cov' modifier cannot coexist with --ibc flag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (calculation_type & CALC_UNRELATED_HERITABILITY) {
            printf("Error: --make-grm 'cov' modifier cannot coexist with --unrelated-heritability\nflag.\n");
            return dispmsg(RET_INVALID_CMDLINE);
          }
          if (ibc_type) {
            printf("Error: --make-grm 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_RELATIONSHIP_COV;
        } else if (!strcmp(argv[cur_arg + kk], "no-gz")) {
          calculation_type &= ~CALC_RELATIONSHIP_GZ;
        } else if ((!strcmp(argv[cur_arg + kk], "ibc1")) || (!strcmp(argv[cur_arg + kk], "ibc2"))) {
          if (calculation_type & CALC_RELATIONSHIP_COV) {
            printf("Error: --make-grm 'cov' modifier cannot coexist with an IBC modifier.%s", errstr_append);
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
    } else if (!strcmp(argptr, "--ibc")) {
      if (calculation_type & CALC_RELATIONSHIP_COV) {
        printf("Error: --ibc flag cannot coexist with a covariance matrix calculation.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 1;
      calculation_type |= CALC_IBC;
    } else if (!strcmp(argptr, "--distance")) {
      if (calculation_type & CALC_LOAD_DISTANCES) {
        printf("Error: --load-dists cannot coexist with a distance matrix calculation flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_DISTANCE_MASK) {
        printf("Error: Duplicate --distance flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 5) {
        printf("Error: --distance accepts at most 5 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      jj = 0;
      for (kk = 1; kk <= ii; kk++) {
        if (!strcmp(argv[cur_arg + kk], "square")) {
          if ((calculation_type & CALC_DISTANCE_SHAPEMASK) == CALC_DISTANCE_SQ0) {
            printf("Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if ((calculation_type & CALC_DISTANCE_SHAPEMASK) == CALC_DISTANCE_TRI) {
            printf("Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if (parallel_tot > 1) {
            printf("Error: --parallel cannot be used with '--distance square'.  Use '--distance\nsquare0' or plain --distance instead.%s", errstr_append);
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
          calculation_type |= CALC_DISTANCE_SQ;
        } else if (!strcmp(argv[cur_arg + kk], "square0")) {
          if ((calculation_type & CALC_DISTANCE_SHAPEMASK) == CALC_DISTANCE_SQ) {
            printf("Error: --distance 'square' and 'square0' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if ((calculation_type & CALC_DISTANCE_SHAPEMASK) == CALC_DISTANCE_TRI) {
            printf("Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_DISTANCE_SQ0;
        } else if (!strcmp(argv[cur_arg + kk], "triangle")) {
          if ((calculation_type & CALC_DISTANCE_SHAPEMASK) == CALC_DISTANCE_SQ) {
            printf("Error: --distance 'square' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          } else if ((calculation_type & CALC_DISTANCE_SHAPEMASK) == CALC_DISTANCE_SQ0) {
            printf("Error: --distance 'square0' and 'triangle' modifiers cannot coexist.%s", errstr_append);
            return dispmsg(RET_INVALID_CMDLINE);
          }
          calculation_type |= CALC_DISTANCE_TRI;
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
	} else if (!strcmp(argv[cur_arg + kk], "ibs")) {
	  if (calculation_type & CALC_PLINK_IBS_MATRIX) {
            printf("Error: --matrix flag cannot be used with '--distance ibs'.%s", errstr_append);
	    return dispmsg(RET_INVALID_CMDLINE);
	  } else if (calculation_type & CALC_DISTANCE_IBS) {
	    printf("Error: Duplicate 'ibs' modifier.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
	  calculation_type |= CALC_DISTANCE_IBS;
	} else if (!strcmp(argv[cur_arg + kk], "1-ibs")) {
	  if (calculation_type & CALC_PLINK_DISTANCE_MATRIX) {
            printf("Error: --matrix flag cannot be used with '--distance 1-ibs'.%s", errstr_append);
	    return dispmsg(RET_INVALID_CMDLINE);
	  } else if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
	    printf("Error: Duplicate '1-ibs' modifier.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
	  calculation_type |= CALC_DISTANCE_1_MINUS_IBS;
	} else if (!strcmp(argv[cur_arg + kk], "snps")) {
	  if (calculation_type & CALC_DISTANCE_SNPS) {
	    printf("Error: Duplicate 'snps' modifier.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
	  calculation_type |= CALC_DISTANCE_SNPS;
        } else {
          printf("Error: Invalid --distance parameter.%s", errstr_append);
          return dispmsg(RET_INVALID_CMDLINE);
        }
      }
      if (((calculation_type & CALC_DISTANCE_BIN) && (!(calculation_type & CALC_DISTANCE_SHAPEMASK))) && (parallel_tot > 1)) {
        printf("Error: --parallel cannot be used with plain '--distance bin'.  Use '--distance\nbin square0' or '--distance bin triangle' instead.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += ii + 1;
      if (!(calculation_type & CALC_DISTANCE_FORMATMASK)) {
        calculation_type |= CALC_DISTANCE_SNPS;
      }
    } else if ((!strcmp(argptr, "--distance-matrix")) || (!strcmp(argptr, "--matrix"))) {
      if (calculation_type & CALC_LOAD_DISTANCES) {
        printf("Error: --load-dists cannot coexist with a distance matrix calculation flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (exponent != 0.0) {
        printf("Error: --exponent flag cannot be used with --distance-matrix or --matrix.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (parallel_tot > 1) {
	printf("Error: --parallel and --distance-matrix/--matrix cannot be used together.  Use\n--distance instead.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 1;
      if (argptr[2] == 'd') {
	if (calculation_type & CALC_DISTANCE_1_MINUS_IBS) {
	  printf("Error: --distance-matrix flag cannot be used with '--distance 1-ibs'.%s", errstr_append);
	  return dispmsg(RET_INVALID_CMDLINE);
	}
        calculation_type |= CALC_PLINK_DISTANCE_MATRIX;
      } else {
	if (calculation_type & CALC_DISTANCE_IBS) {
	  printf("Error: --matrix flag cannot be used with '--distance ibs'.%s", errstr_append);
	  return dispmsg(RET_INVALID_CMDLINE);
	}
        calculation_type |= CALC_PLINK_IBS_MATRIX;
      }
    } else if (!strcmp(argptr, "--load-dists")) {
      if (calculation_type & CALC_GDISTANCE_MASK) {
        printf("Error: --load-dists cannot coexist with a distance matrix calculation flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_LOAD_DISTANCES) {
        printf("Error: Duplicate --load-dists flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (!ii) {
        printf("Error: Missing --load-dists parameter.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 1) {
        printf("Error: Too many --load-dists parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("Error: --load-dists filename too long.\n");
        return dispmsg(RET_OPENFAIL);
      }
      strcpy(loaddistname, argv[cur_arg + 1]);
      cur_arg += 2;
      calculation_type |= CALC_LOAD_DISTANCES;
    } else if (!strcmp(argptr, "--parallel")) {
      if ((calculation_type & CALC_DISTANCE_SHAPEMASK) == CALC_DISTANCE_SQ) {
        printf("Error: --parallel cannot be used with '--distance square'.  Use '--distance\nsquare0' or plain --distance instead.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if ((calculation_type & CALC_DISTANCE_BIN) && (!(calculation_type & CALC_DISTANCE_SHAPEMASK))) {
        printf("Error: --parallel cannot be used with plain '--distance bin'.  Use '--distance\nbin square0' or '--distance bin triangle' instead.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      } else if ((calculation_type & CALC_RELATIONSHIP_SHAPEMASK) == CALC_RELATIONSHIP_SQ) {
        printf("Error: --parallel cannot be used with '--make-rel square'.  Use '--make-rel\nsquare0' or plain '--make-rel' instead.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if ((calculation_type & CALC_RELATIONSHIP_BIN) && (!(calculation_type & CALC_RELATIONSHIP_SHAPEMASK))) {
        printf("Error: --parallel cannot be used with plain '--make-rel bin'.  Use '--make-rel\nbin square0' or '--make-rel bin triangle' instead.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & (CALC_PLINK_DISTANCE_MATRIX | CALC_PLINK_IBS_MATRIX)) {
	printf("Error: --parallel and --distance-matrix/--matrix cannot be used together.  Use\n--distance instead.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_GROUPDIST) {
	printf("Error: --parallel and --groupdist cannot be used together.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_REGRESS_DISTANCE) {
	printf("Error: --parallel and --regress-distance cannot be used together.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_REGRESS_REL) {
	printf("Error: --parallel and --regress-rel cannot be used together.%s", errstr_append);
      } else if (calculation_type & CALC_UNRELATED_HERITABILITY) {
	printf("Error: --parallel and --unrelated-heritability cannot be used together.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_REL_CUTOFF) {
	printf("Error: --parallel cannot be used with --rel-cutoff.  (Use a combination of\n--make-rel, --keep/--remove, and a filtering script.)%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      if (parallel_tot > 1) {
	printf("Error: Duplicate --parallel flag.\n");
	return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii < 2) {
        printf("Error: Not enough --parallel parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 2) {
        printf("Error: Too many --parallel parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      parallel_idx = atoi(argv[cur_arg + 1]);
      if ((parallel_idx < 1) || (parallel_idx > PARALLEL_MAX)) {
        printf("Error: Invalid --parallel job index.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      parallel_idx -= 1; // internal 0..(n-1) indexing
      parallel_tot = atoi(argv[cur_arg + 2]);
      if ((parallel_tot < 2) || (parallel_tot > PARALLEL_MAX) || (parallel_tot < parallel_idx)) {
        printf("Error: Invalid --parallel total job count.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 3;
    } else if (!strcmp(argptr, "--genome")) {
      if (calculation_type & CALC_GENOME) {
	printf("Error: Duplicate --genome flag.\n");
	return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 3) {
	printf("Error: --genome accepts at most 3 parameters.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      for (kk = 1; kk <= ii; kk++) {
	if (!strcmp(argv[cur_arg + kk], "full")) {
	  genome_output_full = 1;
	} else if (!strcmp(argv[cur_arg + kk], "unbounded")) {
	  genome_ibd_unbounded = 1;
	} else if (!strcmp(argv[cur_arg + kk], "gz")) {
	  genome_output_gz = 1;
	} else {
	  printf("Error: Invalid --genome parameter.%s", errstr_append);
	  return dispmsg(RET_INVALID_CMDLINE);
	}
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_GENOME;
    } else if (!strcmp(argptr, "--Z-genome")) {
      if (calculation_type & CALC_GENOME) {
	printf("Error: Duplicate --genome flag.\n");
	return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 2) {
	printf("Error: --Z-genome accepts at most 2 parameters.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      for (kk = 1; kk <= ii; kk++) {
	if (!strcmp(argv[cur_arg + kk], "full")) {
	  genome_output_full = 1;
	} else if (!strcmp(argv[cur_arg + kk], "unbounded")) {
	  genome_ibd_unbounded = 1;
	} else {
	  printf("Error: Invalid --Z-genome parameter.%s", errstr_append);
	  return dispmsg(RET_INVALID_CMDLINE);
	}
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_GENOME;
      genome_output_gz = 1;
    } else if (!strcmp(argptr, "--ppc-gap")) {
      if (ppc_gap != DEFAULT_PPC_GAP) {
	printf("Error: Duplicate --ppc-gap flag.\n");
	return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (!ii) {
	printf("Error: Missing --ppc-gap parameter.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 1) {
	printf("Error: Too many --ppc-gap parameters.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      if ((argv[cur_arg + 1][0] < '0') || (argv[cur_arg + 1][0] > '9')) {
	printf("Error: Invalid --ppc-gap parameter.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      ppc_gap = atoi(argv[cur_arg + 1]);
      if (ppc_gap > 2147483) {
	printf("Error: Invalid --ppc-gap parameter.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      }
      ppc_gap *= 1000;
      cur_arg += 2;
    } else if (!strcmp(argptr, "--groupdist")) {
      if (parallel_tot > 1) {
	printf("Error: --parallel and --groupdist cannot be used together.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_GROUPDIST) {
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
    } else if (!strcmp(argptr, "--regress-distance")) {
      if (parallel_tot > 1) {
	printf("Error: --parallel and --regress-distance cannot be used together.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_REGRESS_DISTANCE) {
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
    } else if (!strcmp(argptr, "--regress-rel")) {
      if (parallel_tot > 1) {
	printf("Error: --parallel and --regress-rel flags cannot coexist.%s", errstr_append);
	return dispmsg(RET_INVALID_CMDLINE);
      } else if (calculation_type & CALC_REGRESS_REL) {
        printf("Error: Duplicate --regress-rel flag.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (ii > 2) {
        printf("Error: --regress-rel accepts at most 2 parameters.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (ii) {
	regress_rel_iters = atoi(argv[cur_arg + 1]);
	if (regress_rel_iters < 2) {
	  printf("Error: Invalid --regress-rel jackknife iteration count.\n");
	  return dispmsg(RET_INVALID_CMDLINE);
	}
        if (ii == 2) {
	  regress_rel_d = atoi(argv[cur_arg + 2]);
	  if (regress_rel_d <= 0) {
	    printf("Error: Invalid --regress-rel jackknife delete parameter.\n");
	    return dispmsg(RET_INVALID_CMDLINE);
	  }
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_REGRESS_REL;
    } else if (!strcmp(argptr, "--unrelated-heritability")) {
      if (calculation_type & CALC_RELATIONSHIP_COV) {
        printf("Error: --unrelated-heritability flag cannot coexist with a covariance\nmatrix calculation.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (parallel_tot > 1) {
	printf("Error: --parallel and --unrelated-heritability cannot be used together.%s", errstr_append);
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
	      printf("Error: Invalid --unrelated-heritability genomic covariance prior.\n");
	      return dispmsg(RET_INVALID_CMDLINE);
	    }
	    if ((unrelated_herit_covg <= 0.0) || (unrelated_herit_covg > 1.0)) {
	      printf("Error: Invalid --unrelated-heritability genomic covariance prior.\n");
	      return dispmsg(RET_INVALID_CMDLINE);
	    }
	    if (ii == jj + 2) {
	      if (sscanf(argv[cur_arg + jj + 2], "%lg", &unrelated_herit_covr) != 1) {
		printf("Error: Invalid --unrelated-heritability residual covariance prior.\n");
		return dispmsg(RET_INVALID_CMDLINE);
	      }
	      if ((unrelated_herit_covr <= 0.0) || (unrelated_herit_covr > 1.0)) {
		printf("Error: Invalid --unrelated-heritability residual covariance prior.\n");
		return dispmsg(RET_INVALID_CMDLINE);
	      }
	    } else {
	      unrelated_herit_covr = 1.0 - unrelated_herit_covg;
	    }
          }
        }
      }
      cur_arg += ii + 1;
      calculation_type |= CALC_UNRELATED_HERITABILITY;
    } else if (!strcmp(argptr, "--freq")) {
      if (freqname[0]) {
        printf("Error: --freq and --read-freq flags cannot coexist.%s", errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (calculation_type & CALC_FREQ) {
        printf("Error: Duplicate --freq flag.\n");
        return dispmsg(RET_INVALID_CMDLINE);
      }
      cur_arg += 1;
      calculation_type |= CALC_FREQ;
    } else if ((!strcmp(argptr, "--read-freq")) || (!strcmp(argptr, "--update-freq"))) {
      if (calculation_type & CALC_FREQ) {
        printf("Error: --freq and %s flags cannot coexist.%s", argptr, errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (freqname[0]) {
        printf("Error: Duplicate %s flag.\n", argptr);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      ii = param_count(argc, argv, cur_arg);
      if (!ii) {
        printf("Error: Missing %s parameter.%s", argptr, errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      } else if (ii > 1) {
        printf("Error: %s accepts only one parameter.%s", argptr, errstr_append);
        return dispmsg(RET_INVALID_CMDLINE);
      }
      if (strlen(argv[cur_arg + 1]) > FNAMESIZE - 1) {
        printf("Error: %s filename too long.\n", argptr);
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
    printf("Error: No output requested.%s", errstr_append);
    return dispmsg(RET_INVALID_CMDLINE);
  }
  if (prune && (!phenoname[0]) && (!ped_col_6)) {
    printf("Error: --prune and --no-pheno cannot coexist without an alternate phenotype\nfile.\n");
    return dispmsg(RET_INVALID_CMDLINE);
  }
  if ((calculation_type & CALC_LOAD_DISTANCES) && (!(calculation_type & (CALC_GROUPDIST | CALC_REGRESS_DISTANCE)))) {
    printf("Error: --load-distances cannot be used without either --groupdist or\n--regress-distance.\n");
    return dispmsg(RET_INVALID_CMDLINE);
  }
  free_cond(subst_argv);

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

  unsigned long buf1[48];
  unsigned long buf2[48];
  unsigned long* bptr = buf1;
  for (ii = 0; ii < 48; ii++) {
    buf1[ii] = 0xaaaaaaaaUL;
    buf2[ii] = 0xaaff55ffUL;
  }
  printf("%u\n", vec_dot_product_biased(&bptr, buf2));

  // famname[0] indicates binary vs. text
  // extractname[0], excludename[0], keepname[0], and removename[0] indicate
  // the presence of their respective flags
  // filtername[0] indicates existence of filter
  // freqname[0] signals --read-freq
  retval = wdist(outname, pedname, mapname, famname, phenoname, extractname, excludename, keepname, removename, filtername, freqname, loaddistname, makepheno_str, filterval, mfilter_col, ped_col_1, ped_col_34, ped_col_5, ped_col_6, ped_col_7, (char)missing_geno, missing_pheno, mpheno_col, phenoname_str, prune, affection_01, chrom_mask, exponent, min_maf, geno_thresh, mind_thresh, hwe_thresh, rel_cutoff, tail_pheno, tail_bottom, tail_top, calculation_type, groupdist_iters, groupdist_d, regress_iters, regress_d, regress_rel_iters, regress_rel_d, unrelated_herit_tol, unrelated_herit_covg, unrelated_herit_covr, ibc_type, parallel_idx, parallel_tot, ppc_gap, nonfounders, genome_output_gz, genome_output_full, genome_ibd_unbounded);
  free(wkspace_ua);
  return dispmsg(retval);
}
