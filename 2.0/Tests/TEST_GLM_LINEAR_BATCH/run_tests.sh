#!/bin/bash

# --glm linear regression on several quantitative phenotypes at once, which
# fits them in batches (GlmLinearBatch()), must report exactly what separate
# runs on each phenotype report (GlmLinear()).
#
# In particular, with multiallelic variants and covariate rows shown, the
# batched writer used to skip the first covariate row of each multiallelic
# variant and repeat the last ALT allele's ADD row in its place, since it
# offset the covariate rows' test indices as if 'intercept' were always
# specified.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# 300 samples, 24 variants; every third variant gets a second ALT allele and
# every sixth a third, carried by a deterministic subset of the samples.
$plink2 --dummy 300 24 acgt --seed 5 --out tmp_dense > /dev/null
$plink2 --pfile tmp_dense --export vcf --out tmp_dense > /dev/null
awk 'BEGIN { OFS = "\t" }
     /^#/ { print; next }
     {
       if (vidx % 3 == 0) {
         extra_ct = (vidx % 6 == 0)? 2 : 1
         split("A C G T", bases, " ")
         delete used
         used[$4] = 1
         used[$5] = 1
         alt_ct = 1
         for (b = 1; (b <= 4) && (alt_ct <= extra_ct); ++b) {
           if (!(bases[b] in used)) {
             used[bases[b]] = 1
             $5 = $5 "," bases[b]
             ++alt_ct
             for (s = 10; s <= NF; ++s) {
               # 7 and 300 are coprime, so this touches a varying subset.
               if (((s - 10) * 7 + vidx * 13 + alt_ct) % 300 < 30) {
                 split($s, gt, "/")
                 $s = gt[1] "/" alt_ct
               }
             }
           }
         }
       }
       ++vidx
       print
     }' tmp_dense.vcf > tmp_multi.vcf
$plink2 --vcf tmp_multi.vcf --make-pgen --out tmp_data > /dev/null
test "$(awk '!/^#/ && $5 ~ /,/' tmp_data.pvar | wc -l)" -eq 8

# Three quantitative phenotypes (none missing, so they can be batched) and
# three covariates.
awk 'BEGIN { srand(7); OFS = "\t"; print "#IID", "Q1", "Q2", "Q3" }
     NR > 1 { print $1, sprintf("%.6f", rand()), sprintf("%.6f", 2 * rand() - 1), sprintf("%.6f", 5 * rand()) }' tmp_data.psam > tmp_pheno.txt
awk 'BEGIN { srand(11); OFS = "\t"; print "#IID", "C1", "C2", "C3" }
     NR > 1 { print $1, sprintf("%.6f", rand() - 0.5), sprintf("%.6f", 3 * rand()), sprintf("%.6f", rand() * rand()) }' tmp_data.psam > tmp_covar.txt

# $1: --glm modifiers.
check_batch() {
    $plink2 --pfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar.txt --glm $1 --out tmp_batch > /dev/null
    grep -q "quantitative phenotypes #1-3" tmp_batch.log
    for q in Q1 Q2 Q3; do
        $plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-name $q --covar tmp_covar.txt --glm $1 --out tmp_single > /dev/null
        cmp tmp_batch.$q.glm.linear tmp_single.$q.glm.linear
    done
}

# Additive model (where a multiallelic variant's alleles share one fit), with
# and without 'intercept' and 'hide-covar'; and the dominant and genotypic
# models, which fit each allele separately.
for mods in "" "omit-ref" "intercept" "omit-ref intercept" "hide-covar" "dominant" "genotypic omit-ref"; do
    check_batch "$mods"
done

# Independently of the single-phenotype writer: with covariate rows shown,
# every multiallelic variant reports each covariate exactly once, and each
# non-omitted allele's ADD row exactly once.
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --covar tmp_covar.txt --glm omit-ref --out tmp_batch > /dev/null
awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; ++i) { c[$i] = i }; next }
     $c["ALT"] ~ /,/ {
       id = $c["ID"]
       if (!(id in alt_ct)) { alt_ct[id] = split($c["ALT"], tmp, ",") }
       key = id "\t" (($c["TEST"] == "ADD")? $c["A1"] : $c["TEST"])
       if (++row_ct[key] > 1) { print "duplicate row " key; exit 1 }
       ++per_variant_ct[id]
     }
     END {
       if (length(alt_ct) != 8) { print "expected 8 multiallelic variants"; exit 1 }
       for (id in alt_ct) {
         if (per_variant_ct[id] != alt_ct[id] + 3) { print id ": " per_variant_ct[id] " rows, expected " (alt_ct[id] + 3); exit 1 }
         for (cv = 1; cv <= 3; ++cv) {
           if (!((id "\tC" cv) in row_ct)) { print id ": no C" cv " row"; exit 1 }
         }
       }
     }' tmp_batch.Q1.glm.linear
