#!/bin/bash

# --glm multinomial: multinomial logistic regression of a categorical
# phenotype, one omnibus test per variant (or per allele).
#
# The fixture (mn.vcf.gz, mn_multi.vcf.gz, pheno.txt, covar.txt, sex.txt)
# comes from make_fixture.py, and the ref_*.txt reference values from
# oracle.py, which fits the same models with statsmodels' MNLogit.  Both need
# numpy, so neither runs here; their outputs are committed.  The fixture has
# 600 samples, a 4-level phenotype, 5 covariates, and 41 variants: ordinary,
# low-frequency, with missing calls, with dosages, one separated by a level
# with no ALT allele, six on chrX (where --glm adds the sex covariate and
# codes male haploid calls 0/2), and four multiallelic ones.  The multiallelic
# variants are in their own file, since VCF import with dosage=DS does not
# accept them.

set -exo pipefail

plink2="$1/plink2 $2 $3"

$plink2 --vcf mn.vcf.gz dosage=DS --update-sex sex.txt --make-pgen --out tmp_data > /dev/null
$plink2 --vcf mn_multi.vcf.gz --update-sex sex.txt --make-pgen --out tmp_multi > /dev/null

# Runs --glm on both filesets and concatenates the reports.
# $1: output name.  Remaining arguments: phenotype/covariate/--glm arguments.
run_both() {
    out=$1
    shift
    $plink2 --pfile tmp_data "$@" --out ${out}_bi > /dev/null
    $plink2 --pfile tmp_multi "$@" --out ${out}_multi > /dev/null
    cat ${out}_bi.*.glm.multinomial > $out.txt
    tail -n +2 ${out}_multi.*.glm.multinomial >> $out.txt
}

# 1. Every model and test, with coefficients, against statsmodels.  In the
#    additive model, a multiallelic variant's ALT alleles are tested jointly;
#    in the others, each ALT allele gets its own row, with the other ALT
#    alleles as nuisance predictors.
for model in add dominant recessive hetonly genotypic hethom; do
    modifier=$model
    if [ $model = add ]; then
        modifier=""
    fi
    for t in lrt score wald; do
        run_both tmp_${model}_$t --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial=$t omit-ref $modifier cols=+beta
        awk -v stat=$(echo $t | tr a-z A-Z) -f compare.awk ref_$model.txt tmp_${model}_$t.txt
    done
done
# The default test is the likelihood ratio test.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref cols=+beta --out tmp_default > /dev/null
cmp tmp_default.CAT.glm.multinomial tmp_add_lrt_bi.CAT.glm.multinomial

# 2. The separated variant must not come out with a p-value.
awk -F '\t' '$3 == "sep1" { if ($NF != "SEPARATION,ALT1") { print "sep1: " $NF; exit 1 }; found = 1 }
             END { if (!found) { print "sep1 missing"; exit 1 } }' tmp_add_lrt.txt

# 3. Without covariates (chrX still gets the sex covariate).
run_both tmp_nocovar --pheno pheno.txt --pheno-name CAT --glm multinomial omit-ref allow-no-covars cols=+beta
awk -v stat=LRT -f compare.awk ref_add_nocovar.txt tmp_nocovar.txt

# 4. A level with a single sample.  On the autosomes, the variants where that
#    sample's dosage is at an extreme of the range are separated, and the rest
#    have ordinary fits.  On chrX, the sex covariate separates that level
#    already (it is empty in one sex), so the covariate-only fit fails: that
#    is an error, or a skipped chromosome with skip-invalid-pheno.
grep -v '^x' ref_single_nocovar.txt > tmp_ref_single_autosomal.txt
run_both tmp_single --not-chr X --pheno pheno.txt --pheno-name SINGLE --glm multinomial omit-ref allow-no-covars cols=+beta
awk -v stat=LRT -f compare.awk tmp_ref_single_autosomal.txt tmp_single.txt
if $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name SINGLE --glm multinomial allow-no-covars --out tmp_single_x > /dev/null 2>&1; then
    echo "chrX covariate-only fit should have failed"
    exit 1
fi
grep -q "covariate-only multinomial logistic regression failed" tmp_single_x.log
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name SINGLE --glm multinomial allow-no-covars skip-invalid-pheno --out tmp_single_x > /dev/null
grep -q "Skipping chrX" tmp_single_x.log
test "$(grep -c '^X' tmp_single_x.SINGLE.glm.multinomial || true)" -eq 0

# 5. With two levels, the model is ordinary logistic regression: coefficients,
#    standard errors, and Wald statistics must match --glm's, for every model.
#    (Multiallelic variants only in the additive model: in the others,
#    logistic regression's standard errors for them differ from statsmodels'
#    by up to 5e-4, while this code's match statsmodels.)
for spec in add:ADD dominant:DOM recessive:REC hetonly:HET genotypic:ADD,DOMDEV hethom:HOM,HET; do
    model=${spec%%:*}
    modifier=$model
    if [ $model = add ]; then
        modifier=""
    fi
    $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name BIN --covar covar.txt --glm no-firth omit-ref hide-covar $modifier cols=+beta --out tmp_logistic > /dev/null
    $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT2 --covar covar.txt --glm multinomial=wald multinomial-ref=lo omit-ref $modifier cols=+beta --out tmp_two_levels > /dev/null
    awk -v terms=${spec#*:} -f logistic_crosscheck.awk tmp_logistic.BIN.glm.logistic tmp_two_levels.CAT2.glm.multinomial
done
$plink2 --pfile tmp_multi --pheno pheno.txt --pheno-name BIN --covar covar.txt --glm no-firth omit-ref hide-covar cols=+beta --out tmp_logistic > /dev/null
$plink2 --pfile tmp_multi --pheno pheno.txt --pheno-name CAT2 --covar covar.txt --glm multinomial=wald multinomial-ref=lo omit-ref cols=+beta --out tmp_two_levels > /dev/null
awk -v terms=ADD -v min_ct=4 -f logistic_crosscheck.awk tmp_logistic.BIN.glm.logistic tmp_two_levels.CAT2.glm.multinomial

# 6. Invariance to the counted allele: swap REF and ALT, keep A1 = ALT.  The
#    statistics stay, the coefficients change sign.
awk '!/^#/ { print $3 "\t" $5 }' tmp_data.pvar > tmp_alt.txt
$plink2 --pfile tmp_data --ref-allele force tmp_alt.txt 2 1 --make-pgen --out tmp_swapped > /dev/null
$plink2 --pfile tmp_swapped --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref cols=+beta --out tmp_swapped > /dev/null
awk -v negate=1 -f same_stats.awk tmp_add_lrt_bi.CAT.glm.multinomial tmp_swapped.CAT.glm.multinomial

# 7. Invariance to the reference level.
run_both tmp_ref_gamma --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref multinomial-ref=gamma
awk -f same_stats.awk tmp_add_lrt.txt tmp_ref_gamma.txt
run_both tmp_ref_gamma_geno --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref genotypic multinomial-ref=gamma
awk -f same_stats.awk tmp_genotypic_lrt.txt tmp_ref_gamma_geno.txt
if $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial multinomial-ref=omega --out tmp_bad_ref > /dev/null 2>&1; then
    echo "unknown multinomial-ref= level should have been rejected"
    exit 1
fi

# 8. Per-level columns: the A1 counts add up, and MIN_EXPECTED is the smallest
#    expected cell of the allele-by-level table.  (Biallelic variants, where
#    each column holds one value.)
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial cols=+a1count,+totallele,+a1countcc,+totallelecc,+a1freqcc --out tmp_cols > /dev/null
awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; ++i) { c[$i] = i }; next }
     {
       a1_sum = 0; obs_sum = 0; min_obs = -1
       split("alpha beta delta gamma", lv, " ")
       for (k = 1; k <= 4; ++k) {
         a1_sum += $c["A1_CT_" lv[k]]
         obs = $c["ALLELE_CT_" lv[k]]
         obs_sum += obs
         if ((min_obs < 0) || (obs < min_obs)) { min_obs = obs }
       }
       if ((a1_sum - $c["A1_CT"]) ^ 2 > 1e-6 || obs_sum != $c["ALLELE_CT"]) { print $3 ": per-level counts do not add up"; exit 1 }
       a1 = $c["A1_CT"]; other = $c["ALLELE_CT"] - a1
       expected = ((a1 < other)? a1 : other) * min_obs / $c["ALLELE_CT"]
       if ((expected - $c["MIN_EXPECTED"]) ^ 2 > (1e-5 * expected) ^ 2 + 1e-12) { print $3 ": MIN_EXPECTED " $c["MIN_EXPECTED"] ", expected " expected; exit 1 }
     }' tmp_cols.CAT.glm.multinomial
# Multiallelic, additive: one comma-separated entry per ALT allele, and
# MIN_EXPECTED over the full allele-by-level table.
$plink2 --pfile tmp_multi --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref cols=+a1count,+totallele,+a1countcc,+totallelecc --out tmp_multi_cols > /dev/null
awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; ++i) { c[$i] = i }; next }
     {
       alt_ct = split($c["A1"], alts, ",")
       if (split($c["A1_CT"], a1, ",") != alt_ct) { print $3 ": A1_CT entry count"; exit 1 }
       split("alpha beta delta gamma", lv, " ")
       min_obs = -1
       for (k = 1; k <= 4; ++k) {
         split($c["A1_CT_" lv[k]], per_level, ",")
         for (j = 1; j <= alt_ct; ++j) { level_sum[j] += per_level[j] }
         obs = $c["ALLELE_CT_" lv[k]]
         if ((min_obs < 0) || (obs < min_obs)) { min_obs = obs }
       }
       min_total = $c["ALLELE_CT"]
       for (j = 1; j <= alt_ct; ++j) {
         if (level_sum[j] != a1[j]) { print $3 ": per-level A1 counts do not add up"; exit 1 }
         if (a1[j] < min_total) { min_total = a1[j] }
         ref_ct -= a1[j]
         level_sum[j] = 0
       }
       ref_ct += $c["ALLELE_CT"]
       if (ref_ct < min_total) { min_total = ref_ct }
       ref_ct = 0
       expected = min_total * min_obs / $c["ALLELE_CT"]
       if ((expected - $c["MIN_EXPECTED"]) ^ 2 > (1e-5 * expected) ^ 2 + 1e-12) { print $3 ": MIN_EXPECTED " $c["MIN_EXPECTED"] ", expected " expected; exit 1 }
     }' tmp_multi_cols.CAT.glm.multinomial

# 9. Results do not depend on the thread count.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref cols=+beta --threads 1 --out tmp_1thread > /dev/null
cmp tmp_1thread.CAT.glm.multinomial tmp_add_lrt_bi.CAT.glm.multinomial
$plink2 --pfile tmp_multi --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial=lrt omit-ref genotypic cols=+beta --threads 1 --out tmp_1thread_multi > /dev/null
cmp tmp_1thread_multi.CAT.glm.multinomial tmp_genotypic_lrt_multi.CAT.glm.multinomial

# 10. --adjust reports every row with a valid statistic, one per ALT allele
#     in the non-additive models.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial --adjust --out tmp_adjust > /dev/null
test "$(tail -n +2 tmp_adjust.CAT.glm.multinomial.adjusted | wc -l)" -eq "$(awk -F '\t' 'NR > 1 && $NF == "."' tmp_adjust.CAT.glm.multinomial | wc -l)"
$plink2 --pfile tmp_multi --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial dominant --adjust --out tmp_adjust_multi > /dev/null
test "$(tail -n +2 tmp_adjust_multi.CAT.glm.multinomial.adjusted | wc -l)" -eq "$(awk -F '\t' 'NR > 1 && $NF == "."' tmp_adjust_multi.CAT.glm.multinomial | wc -l)"

# 11. Without 'multinomial', categorical phenotypes are still skipped, and the
#     modifier leaves other phenotypes' results alone.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name BIN,CAT --covar covar.txt --glm --out tmp_plain > /dev/null
test ! -e tmp_plain.CAT.glm.multinomial
grep -q "Skipping categorical phenotype 'CAT'" tmp_plain.log
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name BIN,CAT --covar covar.txt --glm multinomial --out tmp_with > /dev/null
cmp tmp_plain.BIN.glm.logistic.hybrid tmp_with.BIN.glm.logistic.hybrid
test -e tmp_with.CAT.glm.multinomial

# 12. Modifiers it does not support yet are rejected.
if $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial interaction --out tmp_bad > /dev/null 2>&1; then
    echo "'multinomial interaction' should have been rejected"
    exit 1
fi
