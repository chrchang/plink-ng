#!/bin/bash

# --glm multinomial logistic regression (--mnl-ref).
#
# 1. With two categories the model is ordinary logistic regression, so the
#    per-category rows must match --glm no-firth on the same phenotype coded as
#    case/control.
# 2. Each variant's rows are checked against a run that extracts that variant
#    alone and removes its missing-genotype samples up front (the pattern of
#    TEST_GLM_MISSING), so the null-model refit on each variant's non-missing
#    samples is checked against a fit that never sees the missing calls.
# 3. --threads 1 and --threads 4 must give byte-identical output.
# 4. A fixture whose expected values were computed independently with
#    statsmodels 0.14.6 (MNLogit, Newton-Raphson, full and covariate-only fits
#    per variant), compared at 5 significant digits.
# 5. Multiallelic variants report MULTIALLELIC_UNSUPPORTED, and invalid
#    --mnl-ref usage is rejected.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# The data come from a Park-Miller minimal standard generator, which is exact
# in double-precision arithmetic, so every awk produces the same values and the
# fixture's expected values stay valid.  The phenotype is drawn from integer
# thresholds rather than from exp(), so no libm rounding can move a sample
# between categories.
#
# 400 samples; 20 autosomal variants (every third one carries dosages, which
# are exact multiples of 1/16384 so the import is lossless) and 4 chrX
# variants (outside the pseudoautosomal region).  Variant j blanks 0, 4, 30 or 120 calls, by j % 4.
awk 'function rnd() { seed = (seed * 16807) % 2147483647; return seed / 2147483647 }
     BEGIN {
       seed = 20260923; n = 400; OFS = "\t"
       print "##fileformat=VCFv4.2"
       print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
       print "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage\">"
       print "##contig=<ID=1>"
       print "##contig=<ID=X>"
       line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
       for (s = 1; s <= n; ++s) {
         line = line "\ti" s
         male[s] = (rnd() < 0.5)
       }
       print line
       for (j = 1; j <= 24; ++j) {
         is_x = (j > 20)
         blank_ct = (j % 4 == 0)? 0 : ((j % 4 == 1)? 4 : ((j % 4 == 2)? 30 : 120))
         maf = 0.08 + 0.4 * rnd()
         line = (is_x? "X" : "1") OFS (j * 1000 + is_x * 5000000) OFS (is_x? ("x" (j - 20)) : ("v" j)) OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS "GT:DS"
         for (s = 1; s <= n; ++s) {
           g = (rnd() < maf) + (rnd() < maf)
           if (is_x && male[s]) { g = 2 * (rnd() < maf) }
           cell = ((g == 0)? "0/0" : ((g == 1)? "0/1" : "1/1")) ":" g
           if ((!is_x) && (j % 3 == 1) && (rnd() < 0.4)) {
             k = g * 16384 + int(rnd() * 8001) - 4000
             if (k < 0) { k = 0 }
             if (k > 32768) { k = 32768 }
             cell = ((k < 8192)? "0/0" : ((k < 24576)? "0/1" : "1/1")) ":" sprintf("%.14f", k / 16384)
           }
           if (((s * 37 + j * 11) % n) < blank_ct) { cell = "./.:." }
           line = line OFS cell
           geno[j, s] = g
         }
         print line
       }
       # phenotypes and covariates, written to side files
       print "#IID", "SEX" > "tmp_sex.txt"
       print "#IID", "Q1", "Q2", "GRP" > "tmp_covar.txt"
       print "#IID", "PH4", "B2", "CC", "QT" > "tmp_pheno.txt"
       split("blue gold green red", catname, " ")
       for (s = 1; s <= n; ++s) {
         print "i" s, male[s]? 1 : 2 > "tmp_sex.txt"
         q1 = rnd(); q2 = 3 * rnd() - 1; u = rnd()
         grp = (u < 0.6)? "g1" : ((u < 0.85)? "g2" : "g3")
         print "i" s, sprintf("%.6f", q1), (s % 53 == 0)? "NA" : sprintf("%.6f", q2), grp > "tmp_covar.txt"
         # category thresholds shift with the v2 and v5 genotypes and with Q1
         score = geno[2, s] + geno[5, s] + (q1 > 0.5)
         u = rnd()
         t1 = 0.40 - 0.07 * score; t2 = t1 + 0.20; t3 = t2 + 0.25 - 0.03 * score
         c = (u < t1)? 1 : ((u < t2)? 2 : ((u < t3)? 3 : 4))
         hi = (c == 1) || (c == 4)
         if (rnd() < 0.1) { hi = 1 - hi }
         if (s % 41 == 0) {
           print "i" s, "NONE", "NONE", "NA", "NA" > "tmp_pheno.txt"
         } else {
           print "i" s, catname[c], hi? "hi" : "lo", hi? 2 : 1, sprintf("%.6f", rnd()) > "tmp_pheno.txt"
         }
       }
     }' > tmp_data.vcf
$plink2 --vcf tmp_data.vcf dosage=DS --update-sex tmp_sex.txt --make-pgen --out tmp_data > /dev/null

# 1. K = 2 against logistic regression.
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name CC --covar tmp_covar.txt --glm no-firth --out tmp_logistic > /dev/null
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name B2 --covar tmp_covar.txt --glm --mnl-ref B2=lo --out tmp_k2 > /dev/null
awk -f compare_logistic.awk tmp_logistic.CC.glm.logistic tmp_k2.B2.glm.multinomial

# 2. Missingness self-oracle, on the four-category phenotype.
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --out tmp_main > /dev/null
$plink2 --pfile tmp_data --export A --out tmp_raw > /dev/null
rm -f tmp_miss_*.txt
awk 'NR == 1 { for (j = 7; j <= NF; ++j) { id[j] = $j; sub(/_[^_]*$/, "", id[j]) }; next }
     { for (j = 7; j <= NF; ++j) { if ($j == "NA") { print $2 > ("tmp_miss_" id[j] ".txt") } } }' tmp_raw.raw
grep -v '^#' tmp_data.pvar | cut -f 3 > tmp_variants.txt
rm -f tmp_oracle.PH4.glm.multinomial
for v in $(cat tmp_variants.txt); do
    echo $v > tmp_one.txt
    if [ -f tmp_miss_$v.txt ]; then
        (echo "#IID"; cat tmp_miss_$v.txt) > tmp_remove.txt
        $plink2 --pfile tmp_data --extract tmp_one.txt --remove tmp_remove.txt --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --out tmp_one > /dev/null
    else
        $plink2 --pfile tmp_data --extract tmp_one.txt --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --out tmp_one > /dev/null
    fi
    if [ -f tmp_oracle.PH4.glm.multinomial ]; then
        tail -n +2 tmp_one.PH4.glm.multinomial >> tmp_oracle.PH4.glm.multinomial
    else
        cp tmp_one.PH4.glm.multinomial tmp_oracle.PH4.glm.multinomial
    fi
done
awk -f compare.awk tmp_oracle.PH4.glm.multinomial tmp_main.PH4.glm.multinomial

# 3. Thread count must not change a byte.
$1/plink2 $2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --threads 1 --out tmp_t1 > /dev/null
$1/plink2 $2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --threads 4 --out tmp_t4 > /dev/null
cmp tmp_t1.PH4.glm.multinomial tmp_t4.PH4.glm.multinomial

# 4. statsmodels fixture (autosomes, with the intercept and covariate rows).
$plink2 --pfile tmp_data --chr 1 --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm intercept --mnl-ref PH4=green --out tmp_fixture > /dev/null
awk -f compare_fixture.awk expected.txt tmp_fixture.PH4.glm.multinomial

# 5a. A multiallelic variant gets one row set of NAs.
printf '##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' > tmp_multi.vcf
awk 'NR > 1 { printf "\t%s", $1 }' tmp_data.psam >> tmp_multi.vcf
awk 'BEGIN { printf "\n1\t500\tm1\tA\tC,G\t.\t.\t.\tGT" }
     NR > 1 { split("0/0 0/1 0/2 1/2 1/1", gt, " "); printf "\t%s", gt[(NR * 7) % 5 + 1] }
     END { printf "\n" }' tmp_data.psam >> tmp_multi.vcf
$plink2 --vcf tmp_multi.vcf --make-pgen --out tmp_multi > /dev/null
$plink2 --pfile tmp_multi --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --out tmp_multi > /dev/null
awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; ++i) { if ($i == "ERRCODE") { e = i } }; next }
     { ++n; if ($e != "MULTIALLELIC_UNSUPPORTED") { print "unexpected errcode " $e; exit 1 } }
     END { if (n != 16) { print "expected 16 rows, got " n; exit 1 } }' tmp_multi.PH4.glm.multinomial

# 5b. Invalid usage.  Each of these must fail with the expected message, before
#     any regression output is written.
expect_error() {
    local msg="$1"
    shift
    rm -f tmp_err.*
    if $plink2 --pfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar.txt "$@" --out tmp_err > /dev/null; then
        echo "expected failure: $*"
        exit 1
    fi
    grep -q -- "$msg" tmp_err.log
    if ls tmp_err.*.glm.* > /dev/null 2>&1; then
        echo "output left behind: $*"
        exit 1
    fi
}
expect_error "category 'purple' is not present" --pheno-name PH4 --glm --mnl-ref PH4=purple
# every reference category is checked before the first phenotype is fitted
expect_error "category 'typo' is not present" --pheno-name PH4 B2 --glm --mnl-ref PH4=blue B2=typo
expect_error "'test' column cannot be omitted" --pheno-name PH4 --glm cols=-test --mnl-ref PH4=blue
expect_error "no phenotype named 'NOPE'" --pheno-name PH4 --glm --mnl-ref NOPE=blue
expect_error "phenotype 'QT' is not categorical" --glm --mnl-ref QT=blue
expect_error "phenotype 'PH4' appears more than once" --pheno-name PH4 --glm --mnl-ref PH4=blue PH4=red
expect_error "must be used with --glm" --pheno-name PH4 --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm genotypic --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm interaction --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm firth --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm --parameters 1-2 --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm --tests 1-2 --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm mperm=10 --mnl-ref PH4=blue

# A categorical phenotype without --mnl-ref is still skipped.
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --out tmp_skip > /dev/null
grep -q "Skipping categorical phenotype 'PH4'" tmp_skip.log
test ! -e tmp_skip.PH4.glm.multinomial
