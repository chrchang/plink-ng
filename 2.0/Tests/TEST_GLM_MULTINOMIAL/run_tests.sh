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
# 6. Firth regression, on a second dataset with rare variants, separation and
#    a small category:
#    a. 'firth' against expected values computed independently (direct
#       numerical maximization of the penalized log-likelihood, full and with
#       the genotype coefficients fixed at 0, with scipy), at plink2's print
#       precision;
#    b. the default firth-fallback report: variants flagged FIRTH?=Y are NA
#       under 'no-firth' and identical to the 'firth' rows, the others are
#       identical to the 'no-firth' rows;
#    c. with two categories, the ADD rows match --glm firth logistic;
#    d. --threads 1 and --threads 4 give byte-identical 'firth' and hybrid
#       output;
#    e. missingness self-oracle, as in 2, in 'firth' mode;
#    f. a category emptied by missing calls gives EMPTY_CATEGORY in every
#       mode;
#    g. a covariate that separates a category: 'no-firth' refuses the
#       phenotype, 'firth' and the default use Firth regression throughout.
#
# Parts 1 and 4 use 'no-firth', as they check the unpenalized fit; 2, 3 and 5
# use the default (firth-fallback) mode.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Logistic reference for parts 1 and 6c, fitted one variant at a time with
# --snp and concatenated.  --glm logistic reuses its intercept column from one
# variant to the next, and on master it leaves stale padding in that column
# when a variant with missing calls follows one without (fixed in
# chrchang/plink-ng#526).  Which variants follow which depends on how --threads
# splits the work, so a whole-file run gives different standard errors on
# different machines.  A run of one variant always refills the column.
# Arguments: output file, .pvar file, plink2 arguments (without --out).
logistic_by_variant() {
    out=$1; pvar=$2; shift 2
    rm -f "$out"
    for vid in $(grep -v '^#' "$pvar" | cut -f 3); do
        rm -f tmp_one.*.glm.*
        $plink2 "$@" --snp "$vid" --out tmp_one > /dev/null
        one=$(ls tmp_one.*.glm.*)
        if [ ! -e "$out" ]; then
            head -n 1 "$one" > "$out"
        fi
        tail -n +2 "$one" >> "$out"
    done
}

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
logistic_by_variant tmp_logistic.CC.glm.logistic tmp_data.pvar --pfile tmp_data --pheno tmp_pheno.txt --pheno-name CC --covar tmp_covar.txt --glm no-firth
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name B2 --covar tmp_covar.txt --glm no-firth --mnl-ref B2=lo --out tmp_k2 > /dev/null
awk -f compare_logistic.awk tmp_logistic.CC.glm.logistic tmp_k2.B2.glm.multinomial

# 2. Missingness self-oracle, on the four-category phenotype.
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --out tmp_main > /dev/null
$plink2 --pfile tmp_data --export A --out tmp_raw > /dev/null
rm -f tmp_miss_*.txt
awk 'NR == 1 { for (j = 7; j <= NF; ++j) { id[j] = $j; sub(/_[^_]*$/, "", id[j]) }; next }
     { for (j = 7; j <= NF; ++j) { if ($j == "NA") { print $2 > ("tmp_miss_" id[j] ".txt") } } }' tmp_raw.raw
grep -v '^#' tmp_data.pvar | cut -f 3 > tmp_variants.txt
rm -f tmp_oracle.PH4.glm.multinomial.hybrid
for v in $(cat tmp_variants.txt); do
    echo $v > tmp_one.txt
    if [ -f tmp_miss_$v.txt ]; then
        (echo "#IID"; cat tmp_miss_$v.txt) > tmp_remove.txt
        $plink2 --pfile tmp_data --extract tmp_one.txt --remove tmp_remove.txt --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --out tmp_one > /dev/null
    else
        $plink2 --pfile tmp_data --extract tmp_one.txt --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --out tmp_one > /dev/null
    fi
    if [ -f tmp_oracle.PH4.glm.multinomial.hybrid ]; then
        tail -n +2 tmp_one.PH4.glm.multinomial.hybrid >> tmp_oracle.PH4.glm.multinomial.hybrid
    else
        cp tmp_one.PH4.glm.multinomial.hybrid tmp_oracle.PH4.glm.multinomial.hybrid
    fi
done
awk -f compare.awk tmp_oracle.PH4.glm.multinomial.hybrid tmp_main.PH4.glm.multinomial.hybrid

# 3. Thread count must not change a byte.
$1/plink2 $2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --threads 1 --out tmp_t1 > /dev/null
$1/plink2 $2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --mnl-ref PH4=green --threads 4 --out tmp_t4 > /dev/null
cmp tmp_t1.PH4.glm.multinomial.hybrid tmp_t4.PH4.glm.multinomial.hybrid

# 4. statsmodels fixture (autosomes, with the intercept and covariate rows).
$plink2 --pfile tmp_data --chr 1 --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm intercept no-firth --mnl-ref PH4=green --out tmp_fixture > /dev/null
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
     END { if (n != 16) { print "expected 16 rows, got " n; exit 1 } }' tmp_multi.PH4.glm.multinomial.hybrid

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
expect_error "not supported with multinomial" --pheno-name PH4 --glm --parameters 1-2 --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm --tests 1-2 --mnl-ref PH4=blue
expect_error "not supported with multinomial" --pheno-name PH4 --glm mperm=10 --mnl-ref PH4=blue

# A categorical phenotype without --mnl-ref is still skipped.
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name PH4 --covar tmp_covar.txt --glm --out tmp_skip > /dev/null
grep -q "Skipping categorical phenotype 'PH4'" tmp_skip.log
if ls tmp_skip.PH4.glm.* > /dev/null 2>&1; then
    echo "output written for a skipped categorical phenotype"
    exit 1
fi

# 6. Firth regression.  240 samples (5 with a missing phenotype); PH3 has
# categories a (the reference, ~150 samples), b (~65) and c (~20), shifted by
# Q1 and Q2, and B2 (CC) collapses b and c.  16 variants:
#   r1-r4: MAF ~3%, no A1 allele in category c (separation)
#   r5-r6: 2 or 3 carriers, all in category c (separation)
#   r7-r8: singletons
#   r9-r12: MAF 1-3%, unconstrained
#   r13-r16: MAF 0.15-0.4
# r3, r10 and r14 blank 12 calls (never a carrier's), and r11 carries dosages.
awk 'function rnd() { seed = (seed * 16807) % 2147483647; return seed / 2147483647 }
     BEGIN {
       seed = 918273645; n = 240; OFS = "\t"
       # Phenotype first: category c (small, ~10%), b (~30%), a (the rest),
       # shifted by Q1; B2 collapses {b, c} vs a.
       print "#IID", "PH3", "B2", "CC" > "tmp_firth_pheno.txt"
       print "#IID", "Q1", "Q2" > "tmp_firth_covar.txt"
       for (s = 1; s <= n; ++s) {
         q1 = rnd(); q2 = 2 * rnd() - 1
         u = rnd()
         t1 = 0.06 + 0.08 * q1; t2 = t1 + 0.25 + 0.1 * (q2 > 0)
         cat[s] = (u < t1)? "c" : ((u < t2)? "b" : "a")
         print "f" s, sprintf("%.6f", q1), sprintf("%.6f", q2) > "tmp_firth_covar.txt"
         if (s % 47 == 0) {
           print "f" s, "NONE", "NONE", "NA" > "tmp_firth_pheno.txt"
           cat[s] = "x"
         } else {
           print "f" s, cat[s], (cat[s] == "a")? "lo" : "hi", (cat[s] == "a")? 1 : 2 > "tmp_firth_pheno.txt"
         }
       }
       print "##fileformat=VCFv4.2"
       print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
       print "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage\">"
       print "##contig=<ID=1>"
       line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
       for (s = 1; s <= n; ++s) { line = line "\tf" s }
       print line
       # Variant kinds, by j:
       #   1-4: MAF ~3%, no A1 allele in category c (separation)
       #   5-6: 2 or 3 carriers, all in category c (separation)
       #   7-8: singletons
       #   9-12: MAF 1-3%, unconstrained
       #   13-16: MAF 0.15-0.4
       # Variants 3, 10 and 14 blank 12 calls (never a carrier); variant 11
       # carries dosages (exact multiples of 1/16384).
       for (j = 1; j <= 16; ++j) {
         if (j <= 4) { maf = 0.03 } else if (j <= 8) { maf = 0 } else if (j <= 12) { maf = 0.01 + 0.005 * (j - 9) } else { maf = 0.15 + 0.08 * (j - 13) }
         carriers_left = (j == 5)? 2 : ((j == 6)? 3 : ((j <= 8)? 1 : 0))
         line = "1" OFS (j * 1000) OFS ("r" j) OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS "GT:DS"
         blank_left = ((j == 3) || (j == 10) || (j == 14))? 12 : 0
         for (s = 1; s <= n; ++s) {
           g = (rnd() < maf) + (rnd() < maf)
           if ((j <= 4) && (cat[s] == "c")) { g = 0 }
           if ((j == 5) || (j == 6)) { g = 0; if ((cat[s] == "c") && carriers_left) { g = 1; --carriers_left } }
           if ((j == 7) || (j == 8)) { g = 0; if ((s == 17 * j) && carriers_left) { g = 1; --carriers_left } }
           cell = ((g == 0)? "0/0" : ((g == 1)? "0/1" : "1/1")) ":" g
           if ((j == 11) && g) {
             k = g * 16384 - 1 - int(rnd() * 6000)
             cell = ((k < 24576)? "0/1" : "1/1") ":" sprintf("%.14f", k / 16384)
           }
           if (blank_left && (!g) && (s % 11 == j % 11)) { cell = "./.:."; --blank_left }
           line = line OFS cell
         }
         print line
       }
     }' > tmp_firth.vcf
$plink2 --vcf tmp_firth.vcf dosage=DS --make-pgen --out tmp_firth > /dev/null
firth_args="--pfile tmp_firth --pheno tmp_firth_pheno.txt --covar tmp_firth_covar.txt"

# 6a. Independent fixture: expected_firth.txt was computed by maximizing the
#     penalized log-likelihood l*(b) = l(b) + 0.5 log det I(b) directly
#     (scipy L-BFGS-B, then Newton polishing with a finite-difference Hessian;
#     I built from explicit per-sample Kronecker products and the penalty's
#     gradient from explicit dI/db matrices, so none of plink2's algebra is
#     shared), full and with the genotype coefficients fixed at zero.  The
#     standard errors are those of logistf (inverse information matrix with
#     each sample weighted by 1 + its leverage).  Small likelihood-ratio
#     statistics were recomputed in 40-digit arithmetic.
$plink2 $firth_args --pheno-name PH3 --glm firth intercept --mnl-ref PH3=a --out tmp_ff > /dev/null
awk -f compare_firth_fixture.awk expected_firth.txt tmp_ff.PH3.glm.multinomial.firth

# 6b. Firth-fallback (default) against 'firth' and 'no-firth'.
$plink2 $firth_args --pheno-name PH3 --glm firth --mnl-ref PH3=a --out tmp_fa > /dev/null
$plink2 $firth_args --pheno-name PH3 --glm no-firth --mnl-ref PH3=a --out tmp_fn > /dev/null
$plink2 $firth_args --pheno-name PH3 --glm --mnl-ref PH3=a --out tmp_fh > /dev/null
awk -f compare_hybrid.awk tmp_fa.PH3.glm.multinomial.firth tmp_fn.PH3.glm.multinomial tmp_fh.PH3.glm.multinomial.hybrid
# 'no-firth' reports the separated variants as SEPARATION, and fits the common
# ones
awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; ++i) { if ($i == "ID") { id = i }; if ($i == "ERRCODE") { e = i } }; next }
     ($id ~ /^r[1-8]$/) && ($e !~ /^SEPARATION/) { print "unexpected errcode " $e " on " $id; exit 1 }
     ($id ~ /^r1[3-6]$/) && ($e != ".") { print "unexpected errcode " $e " on " $id; exit 1 }' tmp_fn.PH3.glm.multinomial

# 6c. Two categories: the ADD rows must match --glm firth logistic regression
#     on the same phenotype coded as case/control.  That fit stops once its
#     steps and modified score are below logistf's default tolerances (1e-5),
#     so its coefficients can be about 1e-5 off (relative) from the converged
#     values this one (1e-10) reports: 8.4e-6 on this data, 1.4e-5 and 2.4e-4
#     on other data.  A relative slack of 5e-5 is allowed on top of the
#     rounding allowance: five times that tolerance, well below what a wrong
#     standard-error convention or modified score would produce (percent
#     level).
logistic_by_variant tmp_lf.CC.glm.firth tmp_firth.pvar $firth_args --pheno-name CC --glm firth
$plink2 $firth_args --pheno-name B2 --glm firth --mnl-ref B2=lo --out tmp_k2f > /dev/null
awk -v rel_slack=5e-5 -f compare_logistic.awk tmp_lf.CC.glm.firth tmp_k2f.B2.glm.multinomial.firth

# 6d. Thread count.
for mode in firth hybrid; do
    mod=firth
    if [ $mode = hybrid ]; then
        mod=
    fi
    $1/plink2 $2 $firth_args --pheno-name PH3 --glm $mod --mnl-ref PH3=a --threads 1 --out tmp_ft1 > /dev/null
    $1/plink2 $2 $firth_args --pheno-name PH3 --glm $mod --mnl-ref PH3=a --threads 4 --out tmp_ft4 > /dev/null
    cmp tmp_ft1.PH3.glm.multinomial.$mode tmp_ft4.PH3.glm.multinomial.$mode
done

# 6e. Missingness self-oracle for the variants with missing calls.
$plink2 --pfile tmp_firth --export A --out tmp_fraw > /dev/null
rm -f tmp_fmiss_*.txt tmp_foracle.PH3.glm.multinomial.firth
awk 'NR == 1 { for (j = 7; j <= NF; ++j) { id[j] = $j; sub(/_[^_]*$/, "", id[j]) }; next }
     { for (j = 7; j <= NF; ++j) { if ($j == "NA") { print $2 > ("tmp_fmiss_" id[j] ".txt") } } }' tmp_fraw.raw
printf 'r3\nr10\nr14\n' > tmp_fsubset.txt
for v in r3 r10 r14; do
    echo $v > tmp_one.txt
    (echo "#IID"; cat tmp_fmiss_$v.txt) > tmp_remove.txt
    $plink2 $firth_args --extract tmp_one.txt --remove tmp_remove.txt --pheno-name PH3 --glm firth --mnl-ref PH3=a --out tmp_fone > /dev/null
    if [ -f tmp_foracle.PH3.glm.multinomial.firth ]; then
        tail -n +2 tmp_fone.PH3.glm.multinomial.firth >> tmp_foracle.PH3.glm.multinomial.firth
    else
        cp tmp_fone.PH3.glm.multinomial.firth tmp_foracle.PH3.glm.multinomial.firth
    fi
done
$plink2 $firth_args --extract tmp_fsubset.txt --pheno-name PH3 --glm firth --mnl-ref PH3=a --out tmp_fsub > /dev/null
awk -f compare.awk tmp_foracle.PH3.glm.multinomial.firth tmp_fsub.PH3.glm.multinomial.firth

# 6f. A variant whose missing calls leave category c with no samples is
#     reported as EMPTY_CATEGORY in all three modes, without a fit.
printf '##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' > tmp_empty.vcf
awk 'NR > 1 { printf "\t%s", $1 }' tmp_firth.psam >> tmp_empty.vcf
printf '\n1\t500\te1\tA\tG\t.\t.\t.\tGT' >> tmp_empty.vcf
awk -F '\t' 'NR == FNR { if (FNR > 1) { cat[$1] = $2 }; next }
     FNR > 1 { printf "\t%s", (cat[$1] == "c")? "./." : (((FNR % 5) == 0)? "0/1" : "0/0") }
     END { printf "\n" }' tmp_firth_pheno.txt tmp_firth.psam >> tmp_empty.vcf
$plink2 --vcf tmp_empty.vcf --make-pgen --out tmp_empty > /dev/null
for mode in firth hybrid no-firth; do
    mod=$mode
    if [ $mode = hybrid ]; then
        mod=
    fi
    rm -f tmp_fe.PH3.glm.*
    $plink2 --pfile tmp_empty --pheno tmp_firth_pheno.txt --covar tmp_firth_covar.txt --pheno-name PH3 --glm $mod --mnl-ref PH3=a --out tmp_fe > /dev/null
    awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; ++i) { if ($i == "ERRCODE") { e = i } }; next }
         { ++n; if ($e != "EMPTY_CATEGORY") { print "unexpected errcode " $e; exit 1 } }
         END { if (n != 7) { print "expected 7 rows, got " n; exit 1 } }' tmp_fe.PH3.glm.multinomial*
done

# 6g. A covariate that separates a category: PHS's category c is exactly
#     Q1 > 0.85, so the covariate-only model has no finite maximum-likelihood
#     estimate.  'no-firth' refuses the phenotype; 'firth' and the default fit
#     the covariate-only model with Firth regression instead and use Firth
#     regression for every variant (FIRTH?=Y, rows identical to 'firth').
awk 'BEGIN { OFS = "\t" }
     NR == FNR { if (FNR > 1) { q1[$1] = $2 }; next }
     FNR == 1 { print "#IID", "PHS"; next }
     { v = $2; if (v != "NONE") { v = (q1[$1] > 0.85)? "c" : ((v == "c")? "b" : v) }; print $1, v }' tmp_firth_covar.txt tmp_firth_pheno.txt > tmp_sep_pheno.txt
sep_args="--pfile tmp_firth --pheno tmp_sep_pheno.txt --covar tmp_firth_covar.txt --pheno-name PHS"
rm -f tmp_sep*.PHS.glm.*
if $plink2 $sep_args --glm no-firth --mnl-ref PHS=a --out tmp_sepn > /dev/null; then
    echo "expected 'no-firth' to refuse a covariate-separated phenotype"
    exit 1
fi
grep -q "Firth-fallback was disabled" tmp_sepn.log
$plink2 $sep_args --glm firth --mnl-ref PHS=a --out tmp_sepf > /dev/null
$plink2 $sep_args --glm --mnl-ref PHS=a --out tmp_seph > /dev/null
grep -q "Firth regression will be used for every variant" tmp_seph.log
awk -F '\t' 'NR == FNR { if (FNR > 1) { firth[FNR] = $0 }; next }
     FNR == 1 { for (i = 1; i <= NF; ++i) { if ($i == "FIRTH?") { f = i } }; next }
     { if ($f != "Y") { print "FIRTH?=" $f " on " $3; exit 1 }
       line = ""; for (i = 1; i <= NF; ++i) { if (i != f) { line = line ((line == "")? "" : "\t") $i } }
       if (line != firth[FNR]) { print "hybrid row differs from firth: " $3; exit 1 }
       if ($NF != ".") { print "errcode " $NF " on " $3; exit 1 }
       ++n }
     END { if (n != 112) { print "expected 112 rows, got " n; exit 1 } }' tmp_sepf.PHS.glm.multinomial.firth tmp_seph.PHS.glm.multinomial.hybrid
