#!/bin/bash

# Regression tests for PLINK 1.9 max(T) permutation fixes, plus a
# --linear/--logistic chromosome-skipping fix that these checks uncovered.
# Each check failed before its fix.  Most checks compare --mperm-save-all
# dumps (every permutation statistic) with the .mperm reports and
# .mperm.dump.best files (see check_mperm.awk and check_best.awk).  With
# --threads 1, every run with the same --seed draws the same permutations,
# however they are split into batches, so a dump from one run can be checked
# against the report from another.  That matters because --mperm-save-all
# turns off the precomputed bounds the other runs use.
#
# Usage: ./run_tests.sh
# Requires plink (1.9) to be in the system PATH.

set -uo pipefail

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
# Runs plink with a fixed seed and one thread, discarding the console output.
perm_run() {
    plink "$@" --seed 1 --threads 1 > /dev/null 2>&1
}
# Prints the variant IDs of a .map file.
map_ids() {
    awk '{ print $2 }' "$1"
}

# --assoc/--model max(T): each permutation batch after the first set its
# early-exit threshold from the first batch's maxima instead of its own, so a
# statistic that should have become a permutation's maximum could be skipped.
# Batch sizes depend only on memory here (several thousand permutations even
# with --memory 64, the minimum), and a later permutation's maximum falls
# below the first batch's smallest maximum only about once per batch size, so
# this needs many permutations.  --mperm-save-all disables the threshold, so
# its dump gives the true maxima.
awk -v out=mt -v n=200 -v seed=3 -v mdiff=0 -v sparse=0 \
    -v chrs=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 -f make_data.awk
map_ids tmp_mt.map > tmp_mt_ids.txt
perm_run --file tmp_mt --assoc mperm=60000 perm-count --mperm-save --memory 64 --out tmp_mt_best
perm_run --file tmp_mt --assoc mperm=60000 perm-count --mperm-save-all --memory 64 --out tmp_mt_all
awk -f check_best.awk tmp_mt_all.mperm.dump.all tmp_mt_best.mperm.dump.best
check "--assoc max(T) permutation maxima match --mperm-save-all" $?
awk -f check_mperm.awk tmp_mt_ids.txt tmp_mt_all.mperm.dump.all tmp_mt_best.assoc.mperm
check "--assoc max(T) EMP1/EMP2 match --mperm-save-all" $?
rm -f tmp_mt_all.mperm.dump.all

# --model max(T) with --mperm-save-all: a variant skipped by the trend or
# best-test thread (here, a monomorphic one) didn't write its placeholder
# values, so later variants' permutation statistics landed in the wrong dump
# columns.  Each dump must match its own report, and the report from a run
# without --mperm-save-all.
awk -v out=md -v n=200 -v seed=11 -v mdiff=0 -v sparse=0 \
    -v chrs=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 -f make_data.awk
map_ids tmp_md.map > tmp_md_ids.txt
for mode in trend ""; do
    perm_run --file tmp_md --model mperm=500 $mode perm-count --mperm-save-all --out tmp_md_all
    perm_run --file tmp_md --model mperm=500 $mode perm-count --mperm-save --out tmp_md_best
    report=$(ls tmp_md_all.model.*.mperm)
    md_ok=0
    awk -f check_mperm.awk tmp_md_ids.txt tmp_md_all.mperm.dump.all "$report" || md_ok=1
    awk -f check_mperm.awk tmp_md_ids.txt tmp_md_all.mperm.dump.all "${report/_all/_best}" || md_ok=1
    awk -f check_best.awk tmp_md_all.mperm.dump.all tmp_md_best.mperm.dump.best || md_ok=1
    check "--model ${mode:-best} max(T) report matches --mperm-save-all" $md_ok
    rm -f tmp_md_all.model.*.mperm tmp_md_best.model.*.mperm
done

# --assoc (quantitative) with --pfilter: a variant with fewer than two
# observations skipped its NA in the original-statistics row of the dump.
awk -v out=qt -v n=150 -v seed=5 -v qt=1 -v mdiff=0 -v sparse=4 \
    -v chrs=1,1,1,1,1,1,1,1,1,1,1,1 -f make_data.awk
map_ids tmp_qt.map > tmp_qt_ids.txt
perm_run --file tmp_qt --assoc mperm=500 perm-count --pfilter 0.9 --mperm-save-all --out tmp_qt
awk -f check_mperm.awk tmp_qt_ids.txt tmp_qt.mperm.dump.all tmp_qt.qassoc.mperm
check "--assoc qt --pfilter --mperm-save-all keeps the dump columns aligned" $?

# --linear/--logistic max(T) with variants excluded from analysis (chrX with
# --xchr-model 0, chrY with genotypic coding), here between two autosomes:
#   - The .mperm report listed the excluded variants as if they had been
#     analyzed, which shifted every later variant's results.
#   - Each permutation row of the dump had a value for every variant,
#     including the excluded ones, which were never computed.
#   - With a joint test (genotypic coding), the dump's first row had the ADD
#     statistic, while the permutation rows had the joint test statistic.
# Writes tmp_<$1>.bed/.bim/.fam with chromosome 2 renamed to $2.  ($3: random
# seed, $4: 1 for a quantitative phenotype.)
make_glm_data() {
    awk -v out="$1" -v n=150 -v seed="$3" -v qt="$4" -v mdiff=0 -v sparse=0 \
        -v chrs=1,1,1,1,1,2,2,2,3,3,3,3 -f make_data.awk
    plink --file "tmp_$1" --make-bed --out "tmp_$1" > /dev/null 2>&1
    awk -v c="$2" 'BEGIN { OFS = "\t" } $1 == 2 { $1 = c } { print }' "tmp_$1.bim" > "tmp_$1.bim.new"
    mv "tmp_$1.bim.new" "tmp_$1.bim"
}
# $1: output prefix, $2: TEST column value, $3: linear or logistic
glm_ok() {
    awk -v t="$2" '$5 == t { print $2 }' "$1.assoc.$3" > "$1_ids.txt"
    awk -f check_mperm.awk "$1_ids.txt" "$1.mperm.dump.all" "$1.assoc.$3.mperm"
}
make_glm_data lx X 7 1
perm_run --bfile tmp_lx --linear mperm=200 perm-count --xchr-model 0 --mperm-save-all --out tmp_lx
glm_ok tmp_lx ADD linear
check "--linear --xchr-model 0 max(T) report matches --mperm-save-all" $?
make_glm_data ly Y 9 1
perm_run --bfile tmp_ly --linear genotypic mperm=200 perm-count --mperm-save-all --out tmp_ly
glm_ok tmp_ly GENO_2DF linear
check "--linear genotypic max(T) report matches --mperm-save-all" $?
make_glm_data gx X 7 0
perm_run --bfile tmp_gx --logistic mperm=200 perm-count --xchr-model 0 --mperm-save-all --out tmp_gx
glm_ok tmp_gx ADD logistic
check "--logistic --xchr-model 0 max(T) report matches --mperm-save-all" $?

# --linear/--logistic with --xchr-model 0 (or genotypic, hethom, dominant or
# recessive coding, which force it), when an excluded chromosome comes before
# an analyzed one: the excluded chromosome's variants were analyzed in place
# of the next chromosome's, under their own IDs and the next chromosome's
# code.  The results must match a run without the excluded chromosome.
xy_ok=0
for reg in linear logistic; do
    qt=0
    if [[ $reg == linear ]]; then
        qt=1
    fi
    awk -v out=xy -v n=150 -v seed=21 -v qt=$qt -v mdiff=0 -v sparse=0 \
        -v chrs=1,1,1,X,X,X,XY,XY,XY -f make_data.awk
    plink --file tmp_xy --$reg dominant --out tmp_xy_all > /dev/null 2>&1
    plink --file tmp_xy --chr 1,XY --$reg dominant --out tmp_xy_sub > /dev/null 2>&1
    test -s tmp_xy_sub.assoc.$reg && cmp -s tmp_xy_all.assoc.$reg tmp_xy_sub.assoc.$reg || xy_ok=1
done
check "--linear/--logistic dominant skip chrX before chrXY" $xy_ok

# --test-missing max(T) on chrY (where only males count): the permutation
# counts used the precomputed bounds table of the previous chromosome.
# --mperm-save-all disables the precomputed bounds, so its dump gives
# reference statistics.
awk -v out=tm -v n=300 -v seed=13 -v mdiff=1 -v sparse=0 \
    -v chrs=1,1,1,1,1,1,Y,Y,Y,Y,Y,Y -f make_data.awk
map_ids tmp_tm.map > tmp_tm_ids.txt
perm_run --file tmp_tm --test-missing mperm=1000 perm-count --mperm-save-all --out tmp_tm_all
perm_run --file tmp_tm --test-missing mperm=1000 perm-count --mperm-save --out tmp_tm_best
tm_ok=0
awk -v lo=1 -f check_mperm.awk tmp_tm_ids.txt tmp_tm_all.mperm.dump.all tmp_tm_all.missing.mperm || tm_ok=1
awk -v lo=1 -f check_mperm.awk tmp_tm_ids.txt tmp_tm_all.mperm.dump.all tmp_tm_best.missing.mperm || tm_ok=1
awk -v lo=1 -f check_best.awk tmp_tm_all.mperm.dump.all tmp_tm_best.mperm.dump.best || tm_ok=1
check "--test-missing max(T) with chrY matches --mperm-save-all" $tm_ok

# --test-missing max(T) counts each permutation's case missing calls from
# transposed permutation vectors.  The transposition left half of its 16-byte
# buffer uncleared, and the last word of each variant's missing-call bitfield
# was masked to (sample count mod 32) bits instead of (mod 64).  (300 samples
# exercise both.)  Adaptive permutation counts the same calls directly, and
# with no adaptive check before the last permutation, it draws the same
# permutations and must give the same EMP1 counts.
perm_run --file tmp_tm --test-missing perm perm-count --aperm 5 1000 0.05 1e-4 2000 0 --out tmp_tm_ad
test -s tmp_tm_ad.missing.perm && test "$(awk '{ print $2, $3 }' tmp_tm_ad.missing.perm)" = "$(awk '{ print $2, $3 }' tmp_tm_best.missing.mperm)"
check "--test-missing max(T) EMP1 matches adaptive permutation" $?

exit $failed
