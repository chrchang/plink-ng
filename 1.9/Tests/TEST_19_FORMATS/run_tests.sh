#!/bin/bash

# PLINK 1.9 text-format export and import, checked against plink2.  (Ported
# from the old 1.9/tests_old suite, which compared against PLINK 1.07.)
#
# Usage: ./run_tests.sh {plink 1.9 dir} {plink2 build dir} {up to 2 plink2 args}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  {plink2 build dir} defaults to
# ../../../2.0/build_dynamic.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}
plink2="${2:-../../../2.0/build_dynamic}/plink2 ${3:-} ${4:-}"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
# Prints the genotype matrix of a PLINK 1.9-readable fileset as A1 counts of
# the original data's A1 alleles, so that the comparison doesn't depend on
# allele order.  Samples are identified by IID alone.  $1: output prefix,
# remaining args: input flags.
geno_matrix() {
    local out=$1
    shift
    plink "$@" --recode A --recode-allele tmp_a1.txt --out "$out" > /dev/null 2>&1 || return 1
    awk '{ $1 = ""; $3 = $4 = $5 = $6 = ""; print }' "$out.raw" | tr -s ' ' | sed 's/^ //' | sort
}

# 387 samples (not a multiple of 4, so the last .bed byte of each variant is
# partly padding), 300 variants, 3% missing calls.
plink --dummy 387 300 0.03 acgt --seed 3 --make-bed --out tmp_data > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_data.bim > tmp_data.bim.new
mv tmp_data.bim.new tmp_data.bim
awk '{ print $2, $5 }' tmp_data.bim > tmp_a1.txt
geno_matrix tmp_orig --bfile tmp_data > tmp_orig.txt

# --recode (.ped/.map).  plink2 writes tabs in the .ped.
plink --bfile tmp_data --recode --out tmp_ped19 > /dev/null
$plink2 --bfile tmp_data --export ped --out tmp_ped2 > /dev/null
cmp -s tmp_ped19.ped <(tr '\t' ' ' < tmp_ped2.ped) && cmp -s tmp_ped19.map tmp_ped2.map
check "--recode matches plink2 --export ped" $?

# --recode transpose, and --tfile reading it back.  plink2 writes tabs.
plink --bfile tmp_data --recode transpose --out tmp_t19 > /dev/null
$plink2 --bfile tmp_data --export tped --out tmp_t2 > /dev/null
cmp -s tmp_t19.tped <(tr '\t' ' ' < tmp_t2.tped) && cmp -s tmp_t19.tfam <(tr '\t' ' ' < tmp_t2.tfam)
check "--recode transpose matches plink2 --export tped" $?
geno_matrix tmp_t19_rt --tfile tmp_t19 | cmp -s - tmp_orig.txt
check "--tfile reads back --recode transpose output" $?

# --recode vcf, and --vcf reading it back.
plink --bfile tmp_data --recode vcf --out tmp_v19 > /dev/null
$plink2 --bfile tmp_data --export vcf --out tmp_v2 > /dev/null
cmp -s <(grep -v '^##' tmp_v19.vcf) <(grep -v '^##' tmp_v2.vcf)
check "--recode vcf matches plink2 --export vcf" $?
geno_matrix tmp_v19_rt --vcf tmp_v19.vcf | cmp -s - tmp_orig.txt
check "--vcf reads back --recode vcf output" $?

# --recode oxford, and --data reading it back.  (The .sample phenotype column
# name differs.)
plink --bfile tmp_data --recode oxford --out tmp_o19 > /dev/null
$plink2 --bfile tmp_data --export oxford --out tmp_o2 > /dev/null
cmp -s tmp_o19.gen tmp_o2.gen && cmp -s <(tail -n +2 tmp_o19.sample) <(tail -n +2 tmp_o2.sample)
check "--recode oxford matches plink2 --export oxford" $?
geno_matrix tmp_o19_rt --data tmp_o19 | cmp -s - tmp_orig.txt
check "--data reads back --recode oxford output" $?

exit $failed
