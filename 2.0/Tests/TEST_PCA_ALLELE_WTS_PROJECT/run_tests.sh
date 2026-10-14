#!/bin/bash

# PCA projection with --pca allele-wts and --score variance-standardize, on a
# dataset with multiallelic variants.
#
# Projecting the reference samples themselves must give their PCs back, up
# to one constant factor per PC.  The .eigenvec.allele weights of
# multiallelic variants used to be sqrt(2) too large and of the opposite
# sign relative to the biallelic ones, so as soon as a multiallelic variant
# was present the projection stopped being a multiple of the PCs.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# 160 samples, 300 variants, every fourth one triallelic, no missing calls.
awk 'BEGIN {
  srand(29); OFS = "\t"; n = 160
  print "##fileformat=VCFv4.2"
  print "##contig=<ID=1,length=100000000>"
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  hdr = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
  for (s = 1; s <= n; ++s) { hdr = hdr "\ts" s }
  print hdr
  for (v = 1; v <= 300; ++v) {
    multi = (v % 4 == 0)
    f1 = 0.1 + 0.3 * rand(); f2 = multi? (0.1 + 0.2 * rand()) : 0
    line = "1" OFS (v * 1000) OFS "v" v OFS "A" OFS (multi? "C,G" : "C") OFS "." OFS "." OFS "." OFS "GT"
    for (s = 1; s <= n; ++s) {
      # two subpopulations, so that the PCs carry some structure
      shift = (s <= n / 2)? 0.1 : -0.1
      g = ""
      for (h = 0; h != 2; ++h) {
        r = rand()
        a = (r < f1 + shift * (v % 3 == 0))? 1 : ((r < f1 + f2 + shift * (v % 3 == 0))? 2 : 0)
        g = g ((h)? "/" : "") a
      }
      line = line OFS g
    }
    print line
  }
}' > tmp_data.vcf
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

$plink2 --pfile tmp_data --freq --pca 3 allele-wts --out tmp_ref > /dev/null
test "$(grep -v '^#' tmp_ref.eigenvec.allele | cut -f 2 | uniq -c | awk '$1 == 3' | wc -l)" -eq 75
$plink2 --pfile tmp_data --read-freq tmp_ref.afreq --score tmp_ref.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize cols=+scoresums --score-col-nums 6-8 --out tmp_proj > /dev/null

# |corr(PCk, PCk_SUM)| must be 1 to within printing precision.
awk 'function abs(x) { return (x < 0)? -x : x }
     FNR == 1 { next }
     NR == FNR { for (k = 1; k <= 3; ++k) { e[FNR, k] = $(k + 1) }; next }
     { for (k = 1; k <= 3; ++k) { s[FNR, k] = $(NF - 3 + k) }; n = FNR - 1 }
     END {
       for (k = 1; k <= 3; ++k) {
         se = ss = see = sss = ses = 0
         for (i = 2; i <= n + 1; ++i) {
           x = e[i, k]; y = s[i, k]
           se += x; ss += y; see += x * x; sss += y * y; ses += x * y
         }
         c = (ses - se * ss / n) / sqrt((see - se * se / n) * (sss - ss * ss / n))
         print "PC" k " correlation " c
         if (abs(c) < 0.99999) { exit 1 }
       }
     }' tmp_ref.eigenvec tmp_proj.sscore
