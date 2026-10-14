#!/bin/bash

# A VCF data line with an empty POS field (or one starting with a space) must
# be rejected with the usual "Invalid POS" error.
#
# VcfToPgen() handed POS straight to ScanUintDefcap(), which assumes its first
# character is not a space or delimiter; an empty field tripped its assertion
# and aborted instead.

set -exo pipefail

$1/plink2 $2 $3 --dummy 4 5 --seed 1 --export vcf --out tmp_data

# Sanity check: the unmodified file imports.
$1/plink2 $2 $3 --vcf tmp_data.vcf --make-pgen --out plink2_ok

bad_line=$(awk '!/^#/ { print NR; exit }' tmp_data.vcf)
awk -v n=$bad_line 'BEGIN { FS = OFS = "\t" } NR == n { $2 = "" } { print }' tmp_data.vcf > tmp_empty.vcf
awk -v n=$bad_line 'BEGIN { FS = OFS = "\t" } NR == n { $2 = " " $2 } { print }' tmp_data.vcf > tmp_space.vcf

for f in empty space; do
    if $1/plink2 $2 $3 --vcf tmp_$f.vcf --make-pgen --out plink2_$f 2> tmp_err_$f.txt; then
        echo "expected --vcf to reject tmp_$f.vcf"
        exit 1
    fi
    grep -q "Invalid POS on line $bad_line of --vcf file" tmp_err_$f.txt
done
