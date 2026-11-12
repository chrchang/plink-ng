#!/bin/bash

# Column-value sample/variant filters with no other test coverage:
# --keep-col-match, --keep-col-match-name, --keep-col-match-num, --mfilter,
# --extract-col-cond-match, --extract-col-cond-mismatch, and
# --extract-col-cond-substr (plus --extract-col-cond-min/-max alongside them).

set -exo pipefail

plink2="$1/plink2 $2 $3"

{
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\n'
    for i in 1 2 3 4 5 6; do
        printf "1\t${i}00\tv$i\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\t1/1\t0/1\t0/0\n"
    done
} > tmp_data.vcf
printf '#FID\tIID\tSEX\nF1\tS1\t1\nF1\tS2\t2\nF2\tS3\t1\nF2\tS4\t2\nF3\tS5\t1\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

# $1: expected space-separated IIDs; remaining args: plink2 flags.
check_samples() {
    local expected=$1
    shift
    $plink2 --pfile tmp_data "$@" --make-just-psam --out tmp_s > /dev/null
    test "$(tail -n +2 tmp_s.psam | cut -f 2 | tr '\n' ' ')" = "$expected "
}

# $1: expected space-separated variant IDs; remaining args: plink2 flags.
check_variants() {
    local expected=$1
    shift
    $plink2 --pfile tmp_data "$@" --write-snplist --out tmp_v > /dev/null
    test "$(tr '\n' ' ' < tmp_v.snplist)" = "$expected "
}

# Remaining args: plink2 flags.  The run must fail.
check_fails() {
    if $plink2 --pfile tmp_data "$@" --write-snplist --out tmp_fail > /dev/null 2>&1; then
        exit 1
    fi
}

# --- --keep-col-match ---
# Headerless: FID, IID, then values.  F9/S9 isn't in the dataset, and F3/S3
# matches no sample (S3 is in F2), so neither has any effect.  S5 isn't
# listed at all.
printf 'F1\tS1\tred\tA\nF1\tS2\tblue\tB\nF2\tS3\tred\tB\nF2\tS4\tgreen\tA\nF9\tS9\tred\tA\nF3\tS3\tred\tA\n' > tmp_clu.txt
check_samples "S1 S3" --keep-col-match tmp_clu.txt red
check_samples "S1 S3 S4" --keep-col-match tmp_clu.txt red green
# Exact, case-sensitive matching.
check_fails --keep-col-match tmp_clu.txt RED
check_fails --keep-col-match tmp_clu.txt re
check_samples "S1 S3" --keep-col-match tmp_clu.txt red --keep-col-match-num 3
check_samples "S1 S4" --keep-col-match tmp_clu.txt A --keep-col-match-num 4
# --mfilter n refers to the nth column after the sample IDs.
check_samples "S1 S4" --keep-col-match tmp_clu.txt A --mfilter 2
check_samples "S1 S3" --keep-col-match tmp_clu.txt red --mfilter 1
# Column 2 is a sample ID column.
check_fails --keep-col-match tmp_clu.txt A --keep-col-match-num 2
check_fails --keep-col-match tmp_clu.txt A --keep-col-match-num 2 --mfilter 1

# With a header line.
printf '#FID\tIID\tCOLOR\tGRP\nF1\tS1\tred\tA\nF1\tS2\tblue\tB\nF2\tS3\tred\tB\nF2\tS4\tgreen\tA\n' > tmp_cluh.txt
check_samples "S1 S3" --keep-col-match tmp_cluh.txt red
check_samples "S2 S3" --keep-col-match tmp_cluh.txt B --keep-col-match-name GRP
check_fails --keep-col-match tmp_cluh.txt B --keep-col-match-name NOPE
check_fails --keep-col-match tmp_cluh.txt B --keep-col-match-name GRP --mfilter 2

# --- --extract-col-cond ---
# v6 isn't listed, so every condition excludes it; v7 isn't in the dataset.
printf 'v1\tmissense\t0.5\nv2\tsynonymous\t2\nv3\tmissense_variant\t-1\nv4\tstop_gained\t10\nv5\tSYN\tNA\nv7\tmissense\t1\n' > tmp_ann.txt
# Whole-value matches.
check_variants "v1" --extract-col-cond tmp_ann.txt --extract-col-cond-match missense
check_variants "v1 v4" --extract-col-cond tmp_ann.txt --extract-col-cond-match missense stop_gained
check_variants "v2 v3 v4 v5" --extract-col-cond tmp_ann.txt --extract-col-cond-mismatch missense
# Substring matches (case-sensitive).
check_variants "v1 v3" --extract-col-cond tmp_ann.txt --extract-col-cond-match missense --extract-col-cond-substr
check_variants "v1 v3 v4 v5" --extract-col-cond tmp_ann.txt --extract-col-cond-mismatch syn --extract-col-cond-substr
# Both: must contain a match substring and no mismatch substring.
check_variants "v1 v3 v4" --extract-col-cond tmp_ann.txt --extract-col-cond-match s --extract-col-cond-mismatch nony --extract-col-cond-substr
check_fails --extract-col-cond tmp_ann.txt --extract-col-cond-substr
check_fails --extract-col-cond tmp_ann.txt --extract-col-cond-match missense --extract-col-cond-mismatch missense
check_fails --extract-col-cond-match missense
# Numeric conditions on column 3 (default min 0; NA never passes).
check_variants "v1 v2 v4" --extract-col-cond tmp_ann.txt 3 1
check_variants "v2 v4" --extract-col-cond tmp_ann.txt 3 1 --extract-col-cond-min 1
check_variants "v1 v3" --extract-col-cond tmp_ann.txt 3 1 --extract-col-cond-min -5 --extract-col-cond-max 1
# ID column after the value column, with a skipped header line.
printf 'VAL\tID\nmissense\tv1\nstop_gained\tv4\n' > tmp_ann2.txt
check_variants "v4" --extract-col-cond tmp_ann2.txt 1 2 1 --extract-col-cond-match stop_gained
