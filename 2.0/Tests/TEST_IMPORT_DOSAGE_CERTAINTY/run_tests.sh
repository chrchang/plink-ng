#!/bin/bash

# --import-dosage with --import-dosage-certainty.
#
# A call with no probability above the certainty threshold was marked
# missing, and then overwritten with a dosage anyway.

set -exo pipefail

plink2="$1/plink2 $2 $3"

printf 'S1 S1 0 0 1 -9\nS2 S2 0 0 2 -9\nS3 S3 0 0 1 -9\n' > tmp_data.fam
printf '1\tv1\t0\t100\n1\tv2\t0\t200\n' > tmp_data.map
# format=3: P(hom A1), P(het), P(hom A2).  S1 on v1 has no probability
# >= 0.6.
cat > tmp_data.dosage <<'DOS'
SNP A1 A2 S1 S1 S2 S2 S3 S3
v1 G A 0.5 0.4 0.1 0.9 0.1 0 0 0 1
v2 G A 0 0.95 0.05 1 0 0 0.1 0.2 0.7
DOS
$plink2 --import-dosage tmp_data.dosage format=3 --fam tmp_data.fam --map tmp_data.map --import-dosage-certainty 0.6 --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --export vcf vcf-dosage=DS --out tmp_out > /dev/null
# S1/v1 must be fully missing; everything else is unaffected.
awk '!/^#/ && $3 == "v1" {print $10}' tmp_out.vcf | grep -qx '\./\.'
[[ "$(awk '!/^#/ && $3 == "v1" {print $11}' tmp_out.vcf)" == "1/1:1.9" || "$(awk '!/^#/ && $3 == "v1" {print $11}' tmp_out.vcf)" == "1/1" ]] || exit 1
