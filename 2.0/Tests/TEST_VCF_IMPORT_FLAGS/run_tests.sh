#!/bin/bash

# VCF import flags with no other test coverage: --vcf-min-gq, --vcf-min-dp,
# --vcf-max-dp, --var-min-qual, --vcf-ref-n-missing, --import-max-alleles,
# --vcf-half-call, --polyploid-mode, --vcf-allow-no-nonvar, and
# --import-overlong-var-ids.

set -exo pipefail

plink2="$1/plink2 $2 $3"

header() {
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
}

# $1: input VCF, $2: expected "ID REF ALT GT..." rows joined with '|'
# (tab-separated within a row); remaining args: import flags.  Imports,
# re-exports as VCF, and compares.
check_import() {
    local vcf=$1
    local expected=$2
    shift 2
    $plink2 --vcf $vcf "$@" --export vcf --out tmp_e > /dev/null
    test "$(grep -v '^#' tmp_e.vcf | cut -f 3-5,10- | tr '\n' '|')" = "$(printf "$expected")|"
}

# $1: input VCF; remaining args: import flags.  The import must fail.
check_import_fails() {
    local vcf=$1
    shift
    if $plink2 --vcf $vcf "$@" --make-pgen --out tmp_fail > /dev/null 2>&1; then
        exit 1
    fi
}

# --- Per-genotype GQ/DP filters, QUAL, N REF, and allele count. ---
# q2/S1 has GQ='.', and q2/S2 has no DP field at all; neither should be
# filtered on the missing value.
{
    header
    printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
    printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n'
    printf '1\t100\tq1\tA\tG\t50\t.\t.\tGT:GQ:DP\t0/1:30:10\t1/1:10:10\t0/0:30:50\n'
    printf '1\t200\tq2\tC\tT\t5\t.\t.\tGT:GQ:DP\t0/1:.:3\t1/1:25\t0/1:40:20\n'
    printf '1\t300\tq3\tN\tA\t.\t.\t.\tGT:GQ:DP\t0/1:30:10\t1/1:30:10\t0/0:30:10\n'
    printf '1\t400\tq4\tA\tC,G,T\t40\t.\t.\tGT:GQ:DP\t0/3:30:10\t1/2:30:10\t0/0:30:10\n'
} > tmp_q.vcf
q1='q1\tA\tG\t0/1\t1/1\t0/0'
q2='q2\tC\tT\t0/1\t1/1\t0/1'
q3='q3\tN\tA\t0/1\t1/1\t0/0'
q4='q4\tA\tC,G,T\t0/3\t1/2\t0/0'
check_import tmp_q.vcf "$q1|$q2|$q3|$q4"
check_import tmp_q.vcf 'q1\tA\tG\t0/1\t./.\t0/0'"|$q2|$q3|$q4" --vcf-min-gq 20
check_import tmp_q.vcf "$q1|"'q2\tC\tT\t./.\t1/1\t0/1'"|$q3|$q4" --vcf-min-dp 5
check_import tmp_q.vcf 'q1\tA\tG\t0/1\t1/1\t./.'"|$q2|$q3|$q4" --vcf-max-dp 30
check_import tmp_q.vcf 'q1\tA\tG\t0/1\t1/1\t./.|q2\tC\tT\t./.\t1/1\t0/1'"|$q3|$q4" --vcf-min-dp 5 --vcf-max-dp 30
# Missing QUAL counts as below the threshold.
check_import tmp_q.vcf "$q1|$q4" --var-min-qual 10
check_import tmp_q.vcf "$q1|$q2|$q3" --import-max-alleles 3
check_import tmp_q.vcf "$q1|$q2|$q3|$q4" --import-max-alleles 4

# --vcf-ref-n-missing is only valid on VCFs with the INFO/PR header line that
# .ped-derived VCFs have.
check_import_fails tmp_q.vcf --vcf-ref-n-missing
awk '/^#CHROM/ {print "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">"} {print}' tmp_q.vcf > tmp_qpr.vcf
check_import tmp_qpr.vcf "$q1|$q2|$q3|$q4"
# (VCF export writes a missing REF as N, so check the .pvar.)
$plink2 --vcf tmp_qpr.vcf --vcf-ref-n-missing --make-pgen --out tmp_qpr > /dev/null
test "$(grep -v '^#' tmp_qpr.pvar | cut -f 3-5 | tr '\n' '|')" = "$(printf 'q1\tA\tG|q2\tC\tT|q3\t.\tA|q4\tA\tC,G,T|')"

# --- GT half-calls. ---
{
    header
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n'
    printf '1\t100\th1\tA\tG\t.\t.\t.\tGT\t0/.\t./1\t1|.\t./.\n'
    printf '1\t200\th2\tC\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t.\n'
} > tmp_h.vcf
check_import_fails tmp_h.vcf
check_import_fails tmp_h.vcf --vcf-half-call e
# $1: mode; $2: expected h1 HOM_REF_CT/HET_REF_ALT_CTS/TWO_ALT_GENO_CTS/
# MISSING_CT.  h2 has no half-calls and must be unaffected.
check_half_call() {
    $plink2 --vcf tmp_h.vcf --vcf-half-call $1 --make-pgen --out tmp_h$1 > /dev/null
    $plink2 --pfile tmp_h$1 --geno-counts --out tmp_h$1 > /dev/null
    test "$(grep -v '^#' tmp_h$1.gcount | cut -f 2,5-7,10 | tr '\n' '|')" = "$(printf "h1\t$2|h2\t1\t1\t1\t1|")"
}
# haploid: 0/. -> 0, ./1 -> 1, 1|. -> 1, which are homozygous on chr1.
check_half_call h '1\t0\t2\t1'
check_half_call haploid '1\t0\t2\t1'
check_half_call m '0\t0\t0\t4'
check_half_call missing '0\t0\t0\t4'
# reference: 0/. -> 0/0, ./1 -> 0/1, 1|. -> 1|0 (phase kept).
check_half_call r '1\t2\t0\t1'
check_half_call reference '1\t2\t0\t1'
$plink2 --pfile tmp_hr --export vcf --out tmp_hr > /dev/null
test "$(grep -v '^#' tmp_hr.vcf | head -n 1 | cut -f 12)" = "1|0"

# --- Polyploid genotypes. ---
{
    header
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
    printf '1\t100\tp1\tA\tG\t.\t.\t.\tGT\t0/0/1\t0/1\n'
    printf '1\t200\tp2\tC\tT\t.\t.\t.\tGT\t1/1\t0|1|1|1\n'
} > tmp_p.vcf
check_import_fails tmp_p.vcf
check_import_fails tmp_p.vcf --polyploid-mode e
check_import tmp_p.vcf 'p1\tA\tG\t./.\t0/1|p2\tC\tT\t1/1\t./.' --polyploid-mode m
check_import tmp_p.vcf 'p1\tA\tG\t./.\t0/1|p2\tC\tT\t1/1\t./.' --polyploid-mode missing

# --- Single-sample VCF without any hom-REF calls. ---
# The sanity check kicks in at 1000 variants with nonmissing GT.  $1: variant
# count; $2: index of the variant to give a 0/0 call (0 for none).
make_nonvar_vcf() {
    {
        header
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        awk -v n=$1 -v homref=$2 'BEGIN {for (i = 1; i <= n; i++) printf("1\t%d\tn%d\tA\tG\t.\t.\t.\tGT\t%s\n", i * 10, i, (i == homref)? "0/0" : ((i % 2)? "1/1" : "0/1"))}'
    } > $3
}
make_nonvar_vcf 999 0 tmp_nv999.vcf
make_nonvar_vcf 1000 0 tmp_nv1000.vcf
make_nonvar_vcf 1000 500 tmp_nv1000b.vcf
$plink2 --vcf tmp_nv999.vcf --make-pgen --out tmp_nv999 > /dev/null
check_import_fails tmp_nv1000.vcf
$plink2 --vcf tmp_nv1000.vcf --vcf-allow-no-nonvar --make-pgen --out tmp_nv1000 > /dev/null
test "$(grep -vc '^#' tmp_nv1000.pvar)" -eq 1000
# A single 0/0 call anywhere suppresses the error.
$plink2 --vcf tmp_nv1000b.vcf --make-pgen --out tmp_nv1000b > /dev/null

# --- Overlong variant IDs (>16000 chars). ---
id16000=$(awk 'BEGIN {for (i = 0; i < 1600; i++) printf("abcdefghij")}')
# $1: middle variant's ID; $2: output filename.
make_overlong_vcf() {
    {
        header
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
        printf '1\t10\tok1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n'
        printf '1\t20\t%s\tC\tT\t.\t.\t.\tGT\t1/1\t0/1\n' "$1"
        printf '1\t30\tok3\tG\tA\t.\t.\t.\tGT\t0/0\t1/1\n'
    } > $2
}
make_overlong_vcf "${id16000}X" tmp_ol.vcf
# Exactly 16000 chars is fine.
make_overlong_vcf "$id16000" tmp_ol16000.vcf
$plink2 --vcf tmp_ol16000.vcf --make-pgen --out tmp_ol16000 > /dev/null
test "$(grep -v '^#' tmp_ol16000.pvar | cut -f 3)" = "$(printf 'ok1\n%s\nok3' "$id16000")"
check_import_fails tmp_ol.vcf
check_import_fails tmp_ol.vcf --import-overlong-var-ids e
# $1: mode; $2: expected IDs, one per line.
check_overlong() {
    local mode=$1
    shift
    $plink2 --vcf tmp_ol.vcf --import-overlong-var-ids $mode "$@" --make-pgen --out tmp_ol_$mode > /dev/null
    grep -v '^#' tmp_ol_$mode.pvar | cut -f 3 > tmp_ol_$mode.ids
}
check_overlong m
test "$(cat tmp_ol_m.ids)" = "$(printf 'ok1\n.\nok3')"
check_overlong missing --missing-var-code NOID
test "$(cat tmp_ol_missing.ids)" = "$(printf 'ok1\nNOID\nok3')"
check_overlong s
test "$(cat tmp_ol_s.ids)" = "$(printf 'ok1\nok3')"
test "$(grep -v '^#' tmp_ol_s.pvar | cut -f 2 | tr '\n' ' ')" = "10 30 "
check_overlong t
test "$(cat tmp_ol_t.ids)" = "$(printf 'ok1\n%s\nok3' "$id16000")"
# Genotypes still line up after a skip.
$plink2 --pfile tmp_ol_s --export vcf --out tmp_ol_s > /dev/null
test "$(grep -v '^#' tmp_ol_s.vcf | cut -f 10-11 | tr '\n' '|')" = "$(printf '0/1\t0/0|0/0\t1/1|')"
