#!/bin/bash

# --gxe, checked against PLINK 1.9.

set -exo pipefail

# Deterministic fixture, built with a plain LCG so that every awk produces the
# same files.  The low bits of an LCG have short periods, so values are taken
# from bits 9-30.  ALT frequencies span 0.1 to 0.9, so that ALT is the major
# allele for some variants: that is where PLINK 1.9's minor-allele coding and
# plink2's ALT coding disagree on the sign of BETA.
awk 'BEGIN {
  n = 300
  m = 250
  printf "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
  printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
  for (i = 0; i < n; ++i) { printf "\tS%d", i }
  printf "\n"
  s = 90210
  for (j = 0; j < m; ++j) {
    s = (1103515245 * s + 12345) % 2147483648
    p = 0.1 + (int(s / 512) % 800000) / 1000000.0
    printf "1\t%d\trs%d\tA\tG\t.\t.\t.\tGT", 1000 + j * 100, j
    for (i = 0; i < n; ++i) {
      s = (1103515245 * s + 12345) % 2147483648
      if ((int(s / 512) % 100) < 3) { printf "\t./." ; continue }
      g = 0
      s = (1103515245 * s + 12345) % 2147483648
      if ((int(s / 512) % 1000000) / 1000000.0 < p) { g = g + 1 }
      s = (1103515245 * s + 12345) % 2147483648
      if ((int(s / 512) % 1000000) / 1000000.0 < p) { g = g + 1 }
      if (g == 0) { printf "\t0/0" } else if (g == 1) { printf "\t0/1" } else { printf "\t1/1" }
    }
    printf "\n"
  }
}' > tmp_in.vcf

$1/plink2 $2 $3 --vcf tmp_in.vcf --double-id --make-bed --out tmp_data

# A quantitative phenotype, and a covariate splitting the samples in two.
awk 'BEGIN { s = 271828; print "#FID\tIID\tQT" }
{
  s = (1103515245 * s + 12345) % 2147483648
  u1 = (int(s / 512) % 1000000 + 1) / 1000001.0
  s = (1103515245 * s + 12345) % 2147483648
  u2 = (int(s / 512) % 1000000) / 1000000.0
  printf "%s\t%s\t%.6f\n", $1, $2, sqrt(-2 * log(u1)) * cos(6.283185307 * u2)
}' tmp_data.fam > tmp_pheno.txt

awk 'BEGIN { print "#FID\tIID\tGRP" } { printf "%s\t%s\t%d\n", $1, $2, (NR % 2) + 1 }' tmp_data.fam > tmp_covar.txt

awk 'NR > 1 {print $1, $2, $3}' tmp_pheno.txt > tmp_pheno19.txt
awk 'NR > 1 {print $1, $2, $3}' tmp_covar.txt > tmp_covar19.txt

plink --bfile tmp_data --pheno tmp_pheno19.txt --covar tmp_covar19.txt --gxe \
      --allow-no-sex --out plink19
$1/plink2 $2 $3 --bfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar.txt \
          --gxe --out plink2_gxe
python3 compare_gxe.py plink19.qassoc.gxe plink2_gxe.gxe

# Naming the covariate explicitly must give the same thing as taking the first.
$1/plink2 $2 $3 --bfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar.txt \
          --gxe covar-name=GRP --out plink2_named
diff -q plink2_named.gxe plink2_gxe.gxe

# Zstd output round-trips.
$1/plink2 $2 $3 --bfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar.txt \
          --gxe zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.gxe.zst > plink2_zs.gxe
diff -q plink2_zs.gxe plink2_gxe.gxe

# A covariate with more than two values is an error rather than a guess.
awk 'BEGIN { print "#FID\tIID\tGRP" } { printf "%s\t%s\t%d\n", $1, $2, (NR % 3) + 1 }' tmp_data.fam > tmp_covar3.txt
if $1/plink2 $2 $3 --bfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar3.txt \
        --gxe --out plink2_bad 2> tmp_err.txt; then
    echo "expected --gxe to reject a three-valued covariate"
    exit 1
fi
grep -q "more than two distinct values" tmp_err.txt
