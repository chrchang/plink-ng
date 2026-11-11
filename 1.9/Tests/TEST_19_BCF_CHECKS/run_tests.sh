#!/bin/bash

# PLINK 1.9 --bcf: invalid GT allele indices and out-of-order contig IDX.
#
# Usage: ./run_tests.sh {plink 1.9 dir}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.

set -exo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}

# plink's error exit codes are small; a crash would be 128 or more.
rejects() {
    set +e
    plink --bcf "$1.bcf" --make-bed --allow-extra-chr --out "tmp_$1" > "tmp_$1.log" 2>&1
    local rc=$?
    set -e
    [[ $rc -ne 0 && $rc -lt 128 ]] || exit 1
    grep -q "$2" "tmp_$1.log"
}

python3 make_bcf.py

plink --bcf ok.bcf --make-bed --allow-extra-chr --out tmp_ok
plink --bcf idx_ok.bcf --make-bed --allow-extra-chr --out tmp_idx_ok
[[ $(cut -f1 tmp_idx_ok.bim) == 1 ]] || exit 1

rejects gt8 'GT allele index'
rejects gt32 'GT allele index'
rejects idx_swapped 'contig IDX='

plink --bcf wide_ok.bcf --make-bed --allow-extra-chr --out tmp_wide_ok
rejects wide_gt8 'GT allele index'
rejects wide_gt16 'GT allele index'
plink --bcf many_alleles_ok.bcf --make-bed --allow-extra-chr --out tmp_many_alleles_ok
