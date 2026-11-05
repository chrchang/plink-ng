#!/bin/bash

# PLINK 1.9 --bcf: invalid GT allele indices and out-of-order contig IDX.

set -exo pipefail

# plink's error exit codes are small; a crash would be 128 or more.
rejects() {
    set +e
    plink --bcf "$1.bcf" --make-bed --allow-extra-chr --out "tmp_$1" > "tmp_$1.log" 2>&1
    local rc=$?
    set -e
    [[ $rc -ne 0 && $rc -lt 128 ]]
    grep -q "$2" "tmp_$1.log"
}

python3 make_bcf.py

plink --bcf ok.bcf --make-bed --allow-extra-chr --out tmp_ok
plink --bcf idx_ok.bcf --make-bed --allow-extra-chr --out tmp_idx_ok
[[ $(cut -f1 tmp_idx_ok.bim) == 1 ]]

rejects gt8 'GT allele index'
rejects gt32 'GT allele index'
rejects idx_swapped 'contig IDX='

plink --bcf wide_ok.bcf --make-bed --allow-extra-chr --out tmp_wide_ok
rejects wide_gt8 'GT allele index'
rejects wide_gt16 'GT allele index'
plink --bcf many_alleles_ok.bcf --make-bed --allow-extra-chr --out tmp_many_alleles_ok
