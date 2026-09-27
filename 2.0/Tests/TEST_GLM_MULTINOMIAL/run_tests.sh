#!/bin/bash

# --glm multinomial: multinomial logistic regression of a categorical
# phenotype, one omnibus test per variant.
#
# The fixture (mn.vcf.gz, pheno.txt, covar.txt, sex.txt) comes from
# make_fixture.py, and the ref_*.txt reference values from oracle.py, which
# fits the same models with statsmodels' MNLogit.  Both need numpy, so neither
# runs here; their outputs are committed.  The fixture has 600 samples, a
# 4-level phenotype, 5 covariates, and 37 variants: ordinary, low-frequency,
# with missing calls, with dosages, one separated by a level with no ALT
# allele, and six on chrX (where --glm adds the sex covariate and codes male
# haploid calls 0/2).

set -exo pipefail

plink2="$1/plink2 $2 $3"

$plink2 --vcf mn.vcf.gz dosage=DS --update-sex sex.txt --make-pgen --out tmp_data > /dev/null

# 1. The three tests, with coefficients, against statsmodels.
for t in lrt score wald; do
    $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial=$t omit-ref cols=+beta --out tmp_$t > /dev/null
    awk -v stat=$(echo $t | tr a-z A-Z) -f compare.awk ref_cat.txt tmp_$t.CAT.glm.multinomial
done
# The default test is the likelihood ratio test.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref cols=+beta --out tmp_default > /dev/null
cmp tmp_default.CAT.glm.multinomial tmp_lrt.CAT.glm.multinomial

# 2. The separated variant must not come out with a p-value.
awk -F '\t' '$3 == "sep1" { if ($NF != "SEPARATION,ALT1") { print "sep1: " $NF; exit 1 }; found = 1 }
             END { if (!found) { print "sep1 missing"; exit 1 } }' tmp_lrt.CAT.glm.multinomial

# 3. Without covariates (chrX still gets the sex covariate).
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --glm multinomial omit-ref allow-no-covars cols=+beta --out tmp_nocovar > /dev/null
awk -v stat=LRT -f compare.awk ref_cat_nocovar.txt tmp_nocovar.CAT.glm.multinomial

# 4. A level with a single sample.  On the autosomes, the variants where that
#    sample's dosage is at an extreme of the range are separated, and the rest
#    have ordinary fits.  On chrX, the sex covariate separates that level
#    already (it is empty in one sex), so the covariate-only fit fails: that
#    is an error, or a skipped chromosome with skip-invalid-pheno.
grep -v '^x' ref_single_nocovar.txt > tmp_ref_single_autosomal.txt
$plink2 --pfile tmp_data --not-chr X --pheno pheno.txt --pheno-name SINGLE --glm multinomial omit-ref allow-no-covars cols=+beta --out tmp_single > /dev/null
awk -v stat=LRT -f compare.awk tmp_ref_single_autosomal.txt tmp_single.SINGLE.glm.multinomial
if $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name SINGLE --glm multinomial allow-no-covars --out tmp_single_x > /dev/null 2>&1; then
    echo "chrX covariate-only fit should have failed"
    exit 1
fi
grep -q "covariate-only multinomial logistic regression failed" tmp_single_x.log
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name SINGLE --glm multinomial allow-no-covars skip-invalid-pheno --out tmp_single_x > /dev/null
grep -q "Skipping chrX" tmp_single_x.log
test "$(grep -c '^X' tmp_single_x.SINGLE.glm.multinomial || true)" -eq 0

# 5. Invariance to the counted allele: swap REF and ALT, keep A1 = ALT.  The
#    statistics stay, the coefficients change sign.
awk '!/^#/ { print $3 "\t" $5 }' tmp_data.pvar > tmp_alt.txt
$plink2 --pfile tmp_data --ref-allele force tmp_alt.txt 2 1 --make-pgen --out tmp_swapped > /dev/null
$plink2 --pfile tmp_swapped --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref cols=+beta --out tmp_swapped > /dev/null
awk -v negate=1 -f same_stats.awk tmp_lrt.CAT.glm.multinomial tmp_swapped.CAT.glm.multinomial

# 6. Invariance to the reference level.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref multinomial-ref=gamma --out tmp_ref_gamma > /dev/null
awk -f same_stats.awk tmp_lrt.CAT.glm.multinomial tmp_ref_gamma.CAT.glm.multinomial
head -n 1 tmp_ref_gamma.CAT.glm.multinomial | grep -q "BETA_alpha" && exit 1
if $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial multinomial-ref=omega --out tmp_bad_ref > /dev/null 2>&1; then
    echo "unknown multinomial-ref= level should have been rejected"
    exit 1
fi

# 7. Per-level columns: the A1 counts add up, and MIN_EXPECTED is the smallest
#    expected cell of the allele-by-level table.
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

# 8. Results do not depend on the thread count.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial omit-ref cols=+beta --threads 1 --out tmp_1thread > /dev/null
cmp tmp_1thread.CAT.glm.multinomial tmp_lrt.CAT.glm.multinomial

# 9. --adjust reports every variant with a valid statistic.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial --adjust --out tmp_adjust > /dev/null
test "$(tail -n +2 tmp_adjust.CAT.glm.multinomial.adjusted | wc -l)" -eq "$(awk -F '\t' 'NR > 1 && $NF == "."' tmp_adjust.CAT.glm.multinomial | wc -l)"

# 10. Without 'multinomial', categorical phenotypes are still skipped, and the
#     modifier leaves other phenotypes' results alone.
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name BIN,CAT --covar covar.txt --glm --out tmp_plain > /dev/null
test ! -e tmp_plain.CAT.glm.multinomial
grep -q "Skipping categorical phenotype 'CAT'" tmp_plain.log
$plink2 --pfile tmp_data --pheno pheno.txt --pheno-name BIN,CAT --covar covar.txt --glm multinomial --out tmp_with > /dev/null
cmp tmp_plain.BIN.glm.logistic.hybrid tmp_with.BIN.glm.logistic.hybrid
test -e tmp_with.CAT.glm.multinomial

# 11. Modifiers it does not support yet are rejected.
if $plink2 --pfile tmp_data --pheno pheno.txt --pheno-name CAT --covar covar.txt --glm multinomial genotypic --out tmp_bad > /dev/null 2>&1; then
    echo "'multinomial genotypic' should have been rejected"
    exit 1
fi
