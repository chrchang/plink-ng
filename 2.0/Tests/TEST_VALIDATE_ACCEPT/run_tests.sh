#!/bin/bash

# .pgen encodings that the spec (or, for header modes 8-14, pgenlib) permits,
# but --validate used to reject:
# 1. an LD-compressed record whose base record is all-hom-ref (type 5).  The
#    sparse and raw readers (--sample-counts, --het, --score, --make-pgen)
#    also mishandled this.
# 2. fixed-width storage with a .pvar multiallelic variant that has no
#    ALT2+ calls.
# 5. fused vrtype/length header modes 8-12 and 14 (stale header-size math in
#    the trailing-bit check).
# Each must validate and read exactly like an equivalent common encoding.
#
# Conversely, --validate is allowed to be stricter than the reader, and must
# keep rejecting these (the spec is to be updated to prohibit them):
# 3. multiallelic hardcall track mode byte 0xff.
# 4. phased-dosage existence bitarray with no bits set.

set -exo pipefail

python3 make_pgens.py

# 1.
$1/plink2 $2 $3 --pfile tmp_ldb5 --validate --out tmp_ldb5_v
awk 'NR > 1 && NR <= 601 {print $1}' tmp_ldb0.psam > tmp_keep.txt
printf "v1\tC\t1\n" > tmp_score.txt
for f in tmp_ldb5 tmp_ldb0; do
    $1/plink2 $2 $3 --pfile $f --sample-counts --out ${f}_sc
    $1/plink2 $2 $3 --pfile $f --keep tmp_keep.txt --sample-counts --out ${f}_sck
    $1/plink2 $2 $3 --pfile $f --het --out ${f}_het
    $1/plink2 $2 $3 --pfile $f --score tmp_score.txt 1 2 3 cols=+scoresums --out ${f}_score
    $1/plink2 $2 $3 --pfile $f --keep tmp_keep.txt --score tmp_score.txt 1 2 3 cols=+scoresums --out ${f}_scorek
    $1/plink2 $2 $3 --pfile $f --geno-counts --out ${f}_gc
    $1/plink2 $2 $3 --pfile $f --make-pgen --out ${f}_mp
    $1/plink2 $2 $3 --pfile $f --keep tmp_keep.txt --make-pgen --out ${f}_mpk
    $1/plink2 $2 $3 --pfile ${f}_mp --export A --out ${f}_mp
    $1/plink2 $2 $3 --pfile ${f}_mpk --export A --out ${f}_mpk
done
cmp tmp_ldb5_sc.scount tmp_ldb0_sc.scount
cmp tmp_ldb5_sck.scount tmp_ldb0_sck.scount
cmp tmp_ldb5_het.het tmp_ldb0_het.het
cmp tmp_ldb5_score.sscore tmp_ldb0_score.sscore
cmp tmp_ldb5_scorek.sscore tmp_ldb0_scorek.sscore
cmp tmp_ldb5_gc.gcount tmp_ldb0_gc.gcount
cmp tmp_ldb5_mp.raw tmp_ldb0_mp.raw
cmp tmp_ldb5_mpk.raw tmp_ldb0_mpk.raw

# 2.
$1/plink2 $2 $3 --pfile tmp_fw_multi --validate --out tmp_fw_multi_v

# 3, 4.
if $1/plink2 $2 $3 --pfile tmp_mh_ff --validate --out tmp_mh_ff_v; then
    echo "expected --validate to reject multiallelic mode byte 0xff"
    exit 1
fi
grep -q "Invalid multiallelic hardcall track mode byte (255" tmp_mh_ff_v.log
if $1/plink2 $2 $3 --pfile tmp_dphase_empty --validate --out tmp_dphase_empty_v; then
    echo "expected --validate to reject an empty phased-dosage bitarray"
    exit 1
fi
grep -q "Invalid phased-dosage track" tmp_dphase_empty_v.log

# 5.
for m in 8 9 10 11 12 14; do
    $1/plink2 $2 $3 --pfile tmp_m$m --validate --out tmp_m${m}_v
    $1/plink2 $2 $3 --pfile tmp_m$m --export vcf vcf-dosage=HDS --out tmp_m${m}_e
    $1/plink2 $2 $3 --pfile tmp_s$m --export vcf vcf-dosage=HDS --out tmp_s${m}_e
    diff <(grep -v '^##' tmp_m${m}_e.vcf) <(grep -v '^##' tmp_s${m}_e.vcf)
    if $1/plink2 $2 $3 --pfile tmp_t$m --validate --out tmp_t${m}_v; then
        echo "expected --validate to reject trailing bits in header mode $m"
        exit 1
    fi
    grep -q "Nonzero trailing bits in last vrtype index byte" tmp_t${m}_v.log
done
