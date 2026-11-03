#!/bin/bash

# --epistasis, checked against PLINK 1.9's .epi.qt report where 1.9 can run at
# all, and against an independent least-squares oracle everywhere else.  1.9
# has no covariate support here, so the covariate cases have no other program
# to compare against.

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

plink --simulate-qt simulate.txt --simulate-n 300 --out tmp_data > /dev/null
# Spread the variants over three chromosomes, so that the report's chromosome
# columns are exercised.
awk 'BEGIN{OFS="\t"} {print int((NR - 1) / 50) + 1, $2, 0, ((NR - 1) % 50) * 400000 + 1, $5, $6}' tmp_data.bim > tmp_relabeled.bim
mv tmp_relabeled.bim tmp_data.bim

# 1. The whole report against PLINK 1.9, with --epi1 1 so that the reporting
#    threshold does not come into it.  (1.9's threshold is on the normal
#    approximation to the t-statistic's p-value, this one is on the t p-value
#    itself, so the two select slightly different pairs in between.)
plink --bfile tmp_data --epistasis --epi1 1 --out plink19
$1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 1 --out plink2
awk -f cmp19.awk plink19.epi.qt plink2.PHENO1.epi.qt
awk -f summary19.awk plink19.epi.qt.summary plink2.PHENO1.epi.qt.summary

# 2. The same with missing calls, which each pair has to drop from its own
#    regression.
plink --simulate-qt simulate.txt --simulate-n 300 --simulate-missing 0.08 --out tmp_miss > /dev/null
plink --bfile tmp_miss --epistasis --epi1 1 --out plink19_m
$1/plink2 $2 $3 --bfile tmp_miss --epistasis --epi1 1 --out plink2_m
awk -f cmp19.awk plink19_m.epi.qt plink2_m.PHENO1.epi.qt
awk -f summary19.awk plink19_m.epi.qt.summary plink2_m.PHENO1.epi.qt.summary

# 3. Against the oracle, which solves each pair's least squares from scratch
#    from scratch.  A smaller fileset, since the oracle is quadratic and slow.
plink --simulate-qt simulate_small.txt --simulate-n 300 --simulate-missing 0.05 --out tmp_small > /dev/null
$1/plink2 $2 $3 --bfile tmp_small --epistasis --epi1 1 --out plink2_s
python3 oracle.py tmp_small > oracle_nocovar.txt
awk -f cmp_oracle.awk oracle_nocovar.txt plink2_s.PHENO1.epi.qt

# 4. Covariates: two quantitative and one binary.  PLINK 1.9 cannot run this,
#    so the oracle is the only check.
awk 'BEGIN {srand(7); print "#FID\tIID\tage\tpc1\tgrp"}
     {printf "%s\t%s\t%.4f\t%.5f\t%d\n", $1, $2, 20 + 60 * rand(), rand() * 2 - 1, (rand() < 0.5)? 0 : 1}' tmp_small.fam > tmp_small.cov
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small.cov --epistasis --epi1 1 --out plink2_c
python3 oracle.py tmp_small tmp_small.cov > oracle_covar.txt
awk -f cmp_oracle.awk oracle_covar.txt plink2_c.PHENO1.epi.qt

# 5. Missing covariate values drop those samples from every pair.
awk 'BEGIN {srand(9)} NR == 1 {print; next} {if (rand() < 0.1) { $4 = "NA" }; print}' OFS='\t' tmp_small.cov > tmp_small_na.cov
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small_na.cov --epistasis --epi1 1 --out plink2_na
python3 oracle.py tmp_small tmp_small_na.cov > oracle_na.txt
awk -f cmp_oracle.awk oracle_na.txt plink2_na.PHENO1.epi.qt

# 6. A constant covariate is dropped with a warning rather than failing every
#    pair, so the run reproduces the one without it.
awk 'NR == 1 {print $0 "\tconst"; next} {print $0 "\t1"}' OFS='\t' tmp_small.cov > tmp_small_const.cov
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small_const.cov --epistasis --epi1 1 --out plink2_k
grep -q "Excluding constant covariate 'const'" plink2_k.log
diff -q plink2_c.PHENO1.epi.qt plink2_k.PHENO1.epi.qt

# 7. --epi1 filters the report without changing the summary's N_TOT.
$1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 0.001 --out plink2_p1
test "$(grep -vc '^#' plink2_p1.PHENO1.epi.qt)" -lt "$(grep -vc '^#' plink2.PHENO1.epi.qt)"
diff <(cut -f2,4 plink2_p1.PHENO1.epi.qt.summary) <(cut -f2,4 plink2.PHENO1.epi.qt.summary)

# 8. --threads 1 reproduces the default run, --parallel chunks concatenate to
#    it, 'zs' round-trips, and 'cols=-p' drops the p-value column.
$1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 1 --threads 1 --out plink2_st
diff -q plink2.PHENO1.epi.qt plink2_st.PHENO1.epi.qt
diff -q plink2.PHENO1.epi.qt.summary plink2_st.PHENO1.epi.qt.summary

# The threads split each row's columns, so the covariate case has to come out
# the same way too.
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small.cov --epistasis --epi1 1 --threads 1 --out plink2_ct
diff -q plink2_c.PHENO1.epi.qt plink2_ct.PHENO1.epi.qt
diff -q plink2_c.PHENO1.epi.qt.summary plink2_ct.PHENO1.epi.qt.summary

for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --epistasis --epi1 1 --parallel $i 3 --out plink2_par$i
done
cat plink2_par1.PHENO1.epi.qt.1 plink2_par2.PHENO1.epi.qt.2 plink2_par3.PHENO1.epi.qt.3 > plink2_par.PHENO1.epi.qt
diff -q plink2.PHENO1.epi.qt plink2_par.PHENO1.epi.qt

$1/plink2 $2 $3 --bfile tmp_data --epistasis zs --epi1 1 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.PHENO1.epi.qt.zst > plink2_zs.PHENO1.epi.qt
diff -q plink2.PHENO1.epi.qt plink2_zs.PHENO1.epi.qt

$1/plink2 $2 $3 --bfile tmp_data --epistasis cols=-p --epi1 1 --out plink2_nop
head -n 1 plink2_nop.PHENO1.epi.qt | grep -qx '#CHROM1	ID1	CHROM2	ID2	BETA_INT	SE	T_STAT'

# 9. Rejected: a case/control phenotype, and the two epistasis commands
#    together.  (The set modifiers are checked in section 14.)
plink --simulate simulate.txt --simulate-ncases 150 --simulate-ncontrols 150 --out tmp_cc > /dev/null
fails $1/plink2 $2 $3 --bfile tmp_cc --epistasis --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis --epistasis-boost --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis --fast-epistasis boost --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis --gap 100 --out plink2_bad

# 10. Every row has as many tab-separated fields as the header, with and
#     without the A1 columns.  (The comparisons above split on runs of
#     whitespace, so they cannot see an empty field.)
for c in "" "cols=+a1" "cols=+pos" "cols=-chrom"
do
$1/plink2 $2 $3 --bfile tmp_data --epistasis $c --epi1 1 --out plink2_fields
awk -F'\t' 'NR == 1 {n = NF} NF != n {exit 1}' plink2_fields.PHENO1.epi.qt
done
head -n 2 plink2.PHENO1.epi.qt | tail -n 1 | awk -F'\t' '{exit ($5 == "")}'

# 11. 'log10' reports -log10(p) under a NEG_LOG10_P header.
$1/plink2 $2 $3 --bfile tmp_data --epistasis log10 --epi1 1 --out plink2_log10
head -n 1 plink2_log10.PHENO1.epi.qt | grep -qx '#CHROM1	ID1	CHROM2	ID2	BETA_INT	SE	T_STAT	NEG_LOG10_P'
paste <(cut -f 8 plink2.PHENO1.epi.qt) <(cut -f 8 plink2_log10.PHENO1.epi.qt) | awk 'NR > 1 {d = -log($1) / log(10) - $2; if (d < 0) {d = -d} if (d > 1e-4 * (1 + $2)) {exit 1}}'

# 12. Sample counts whose bit arrays take an odd number (7 or more) of words:
#     the three per-variant bitvectors were packed at that unaligned stride,
#     and PopcountWords() tripped its alignment assertion.
for n in 400 448
do
$1/plink2 $2 $3 --dummy $n 20 0.1 scalar-pheno --seed 1 --epistasis --epi1 1 --out plink2_n$n
done

# 13. 'ref-based' is accepted, as documented.  With REF as the major allele of
#     every variant, it must not change the report.
$1/plink2 $2 $3 --bfile tmp_data --epistasis ref-based --epi1 1 --out plink2_rb
$1/plink2 $2 $3 --bfile tmp_data --freq --out plink2_rb_freq
if awk 'NR > 1 && $6 > 0.5 {found = 1} END {exit !found}' plink2_rb_freq.afreq; then
    echo "some ALT allele is the major one; skipping the comparison"
else
    diff -q plink2.PHENO1.epi.qt plink2_rb.PHENO1.epi.qt
fi

# 14. The set modes, which PLINK 1.9 also has, so they are checked against it,
#     and against the matching pairs of the all-pairs report.  setA is the
#     first 20 variants (chromosome 1), setB 15 variants on chromosome 2, and
#     setC straddles setA so that two sets can overlap.
{
    printf 'setA\n'
    awk 'NR <= 20 {print $2}' tmp_data.bim
    printf 'END\n\nsetB\n'
    awk 'NR > 50 && NR <= 65 {print $2}' tmp_data.bim
    printf 'END\n\nsetC\n'
    awk 'NR > 10 && NR <= 30 {print $2}' tmp_data.bim
    printf 'END\n'
} > tmp_sets.txt

# 14a. One set: every pair inside it, which is still a triangle.
plink --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis set-by-set --epi1 1 --out plink19_sbs
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis set-by-set --epi1 1 --out plink2_sbs
awk -f cmp19.awk plink19_sbs.epi.qt plink2_sbs.PHENO1.epi.qt
awk -f summary19.awk plink19_sbs.epi.qt.summary plink2_sbs.PHENO1.epi.qt.summary
awk -v expected_n=190 -f cmp_subset.awk plink2.PHENO1.epi.qt plink2_sbs.PHENO1.epi.qt

# 14b. Two sets: every pair across them.
plink --bfile tmp_data --set tmp_sets.txt --set-names setA setB --epistasis set-by-set --epi1 1 --out plink19_2s
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA setB --epistasis set-by-set --epi1 1 --out plink2_2s
awk -f cmp19.awk plink19_2s.epi.qt plink2_2s.PHENO1.epi.qt
# 1.9 only summarizes the row set; this also summarizes the column set, as
# --epistasis-boost does, so compare the row set's variants.
awk 'NR == FNR {k[$2]; next} /^#/ || ($2 in k)' plink19_2s.epi.qt.summary plink2_2s.PHENO1.epi.qt.summary > plink2_2s_rows.summary
awk -f summary19.awk plink19_2s.epi.qt.summary plink2_2s_rows.summary
awk -v expected_n=300 -f cmp_subset.awk plink2.PHENO1.epi.qt plink2_2s.PHENO1.epi.qt

# 14c. set-by-all: the set against every variant, minus the self-pairs.  A pair
#      inside the set is listed once from each side, as in 1.9.
plink --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis set-by-all --epi1 1 --out plink19_sba
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis set-by-all --epi1 1 --out plink2_sba
awk -f cmp19.awk plink19_sba.epi.qt plink2_sba.PHENO1.epi.qt
# (The summaries are not compared here: a variant inside the set is also a
# column for the other rows, which only this program counts.)
awk -v expected_n=2980 -f cmp_subset.awk plink2.PHENO1.epi.qt plink2_sba.PHENO1.epi.qt

# 14d. Overlapping sets do not pair a variant with itself.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA setC --epistasis set-by-set --epi1 1 --out plink2_ov
awk -v expected_n=390 -f cmp_subset.awk plink2.PHENO1.epi.qt plink2_ov.PHENO1.epi.qt
awk 'NR > 1 && $2 == $4 {exit 1}' plink2_ov.PHENO1.epi.qt

# 14e. With covariates and missing calls, against the all-pairs run of the same
#      data.
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small.cov --epistasis --epi1 1 --out plink2_cfull
{ printf 'setA\n'; awk 'NR <= 8 {print $2}' tmp_small.bim; printf 'END\n'; } > tmp_small_sets.txt
n_small="$(awk 'END {print NR}' tmp_small.bim)"
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small.cov --set tmp_small_sets.txt --epistasis set-by-all --epi1 1 --out plink2_csba
awk -v expected_n=$((8 * (n_small - 1))) -f cmp_subset.awk plink2_cfull.PHENO1.epi.qt plink2_csba.PHENO1.epi.qt
$1/plink2 $2 $3 --bfile tmp_small --covar tmp_small.cov --set tmp_small_sets.txt --epistasis set-by-set --epi1 1 --out plink2_csbs
awk -v expected_n=28 -f cmp_subset.awk plink2_cfull.PHENO1.epi.qt plink2_csbs.PHENO1.epi.qt

# 14f. One thread reproduces the default run, and --parallel chunks
#      concatenate to it, in each mode.
for m in "set-by-set" "set-by-all"
do
    t="${m//-/}"
    $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis $m --epi1 1 --out plink2_$t
    $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis $m --epi1 1 --threads 1 --out plink2_${t}_st
    diff -q plink2_$t.PHENO1.epi.qt plink2_${t}_st.PHENO1.epi.qt
    diff -q plink2_$t.PHENO1.epi.qt.summary plink2_${t}_st.PHENO1.epi.qt.summary
    for i in 1 2 3
    do
        $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA --epistasis $m --epi1 1 --parallel $i 3 --out plink2_${t}_p$i
    done
    cat plink2_${t}_p1.PHENO1.epi.qt.1 plink2_${t}_p2.PHENO1.epi.qt.2 plink2_${t}_p3.PHENO1.epi.qt.3 > plink2_${t}_p.PHENO1.epi.qt
    diff -q plink2_$t.PHENO1.epi.qt plink2_${t}_p.PHENO1.epi.qt
done
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA setB --epistasis set-by-set --epi1 1 --threads 1 --out plink2_2s_st
diff -q plink2_2s.PHENO1.epi.qt plink2_2s_st.PHENO1.epi.qt
for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names setA setB --epistasis set-by-set --epi1 1 --parallel $i 3 --out plink2_2s_p$i
done
cat plink2_2s_p1.PHENO1.epi.qt.1 plink2_2s_p2.PHENO1.epi.qt.2 plink2_2s_p3.PHENO1.epi.qt.3 > plink2_2s_p.PHENO1.epi.qt
diff -q plink2_2s.PHENO1.epi.qt plink2_2s_p.PHENO1.epi.qt

# 14g. Rejected: a set mode with no --set, the two modes together, set-by-all
#      with more than one set, and set-by-set with more than two.
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis set-by-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --epistasis set-by-all --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --epistasis set-by-set set-by-all --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --epistasis set-by-all --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --epistasis set-by-set --out plink2_bad

# 15. Several quantitative phenotypes in one run.  Phenotypes with the same
#     missingness pattern share each pair's design-matrix work; the others get
#     their own pass.  Either way, each phenotype's report and summary must be
#     byte-identical to a single-phenotype run that selects it with
#     --pheno-name.  B and C share a missingness pattern (so go through the
#     shared path together), A is complete and D has a pattern of its own (so
#     each is a pass by itself), and CC, a case/control column, is skipped.
#     Covariates with a missing value intersect all four patterns.
awk 'BEGIN{OFS="\t"}
     {
       a = $6
       b = ((NR * 7) % 13) * 0.3 - 1.1 + 0.4 * $6
       c = ((NR * 5) % 11) * 0.25 + 0.7 * $6
       d = ((NR * 3) % 17) * 0.2
       cc = (NR % 3 == 0) ? 2 : 1
       if (NR % 9 == 0) { b = "NA"; c = "NA" }
       if (NR % 14 == 0) { d = "NA" }
       print $1, $2, a, b, c, d, cc
     }' tmp_small.fam > tmp_multi_body.txt
(printf '#FID\tIID\tA\tB\tC\tD\tCC\n'; cat tmp_multi_body.txt) > tmp_multi_pheno.txt
for cov in "" "--covar tmp_small.cov" "--covar tmp_small_na.cov"
do
    $1/plink2 $2 $3 --bfile tmp_small --pheno tmp_multi_pheno.txt $cov --epistasis --epi1 1 --out plink2_multi
    test ! -e plink2_multi.CC.epi.qt
    for p in A B C D
    do
        $1/plink2 $2 $3 --bfile tmp_small --pheno tmp_multi_pheno.txt --pheno-name $p $cov --epistasis --epi1 1 --out plink2_single_$p
        cmp plink2_multi.$p.epi.qt plink2_single_$p.$p.epi.qt
        cmp plink2_multi.$p.epi.qt.summary plink2_single_$p.$p.epi.qt.summary
        test "$(wc -l < plink2_multi.$p.epi.qt)" -gt 100
    done
    if cmp -s plink2_multi.A.epi.qt plink2_multi.B.epi.qt; then exit 1; fi
    # B and C always share a pass.
    grep -q "Processing [0-9]* phenotypes with identical missingness patterns together" plink2_multi.log
    # Thread count changes nothing, in the shared path either.
    $1/plink2 $2 $3 --bfile tmp_small --pheno tmp_multi_pheno.txt $cov --epistasis --epi1 1 --threads 1 --out plink2_multi_t1
    $1/plink2 $2 $3 --bfile tmp_small --pheno tmp_multi_pheno.txt $cov --epistasis --epi1 1 --threads 8 --out plink2_multi_t8
    for p in A B C D
    do
        cmp plink2_multi.$p.epi.qt plink2_multi_t1.$p.epi.qt
        cmp plink2_multi.$p.epi.qt plink2_multi_t8.$p.epi.qt
        cmp plink2_multi.$p.epi.qt.summary plink2_multi_t8.$p.epi.qt.summary
    done
    # --parallel jobs concatenate to the single run, for every phenotype.
    for j in 1 2 3
    do
        $1/plink2 $2 $3 --bfile tmp_small --pheno tmp_multi_pheno.txt $cov --epistasis --epi1 1 --parallel $j 3 --out plink2_multi_par
    done
    for p in A B C D
    do
        cat plink2_multi_par.$p.epi.qt.1 plink2_multi_par.$p.epi.qt.2 plink2_multi_par.$p.epi.qt.3 > tmp_multi_par_cat.txt
        cmp plink2_multi.$p.epi.qt tmp_multi_par_cat.txt
    done
done
# Set modes and compressed output go through the same path.
{
    printf 'setA\n'
    awk 'NR <= 20 {print $2}' tmp_small.bim
    printf 'END\n\nsetB\n'
    awk 'NR > 10 && NR <= 30 {print $2}' tmp_small.bim
    printf 'END\n'
} > tmp_multi_sets.txt
$1/plink2 $2 $3 --bfile tmp_small --set tmp_multi_sets.txt --pheno tmp_multi_pheno.txt --covar tmp_small.cov --epistasis set-by-set zs --epi1 1 --out plink2_multi_set
$1/plink2 $2 $3 --bfile tmp_small --set tmp_multi_sets.txt --pheno tmp_multi_pheno.txt --pheno-name C --covar tmp_small.cov --epistasis set-by-set zs --epi1 1 --out plink2_single_set
$1/plink2 $2 $3 --zst-decompress plink2_multi_set.C.epi.qt.zst > tmp_multi_set_C.txt
$1/plink2 $2 $3 --zst-decompress plink2_single_set.C.epi.qt.zst > tmp_single_set_C.txt
cmp tmp_multi_set_C.txt tmp_single_set_C.txt
$1/plink2 $2 $3 --zst-decompress plink2_multi_set.D.epi.qt.summary.zst > tmp_multi_set_D.txt
$1/plink2 $2 $3 --bfile tmp_small --set tmp_multi_sets.txt --pheno tmp_multi_pheno.txt --pheno-name D --covar tmp_small.cov --epistasis set-by-set --epi1 1 --out plink2_single_set_D
cmp tmp_multi_set_D.txt plink2_single_set_D.D.epi.qt.summary
# Only a case/control phenotype selected: still an error.
fails $1/plink2 $2 $3 --bfile tmp_small --pheno tmp_multi_pheno.txt --pheno-name CC --epistasis --out plink2_bad
