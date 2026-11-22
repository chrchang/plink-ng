#!/bin/bash

# --score mean imputation of missing genotypes under the 'dominant' and
# 'recessive' modifiers.
#
# A missing genotype must contribute the expected coded value: 2f for the
# default additive model, 1 - (1-f)^2 for 'dominant' and f^2 for 'recessive',
# where f is the frequency of the scored allele.  Every SCORE1_SUM is checked
# against oracle.awk, which recomputes it from the VCF, for all three models
# with and without 'no-mean-imputation'.  Half the variants score the REF
# allele.  The run is repeated on a copy of the data with fractional dosages,
# since --score takes a different code path when dosages are present.

set -exo pipefail

plink2="$1/plink2 $2 $3"

$plink2 --dummy 300 40 acgt --seed 3 --out tmp_dense > /dev/null
$plink2 --pfile tmp_dense --export vcf --out tmp_dense > /dev/null

# Blank 0, 6, 45 or 120 of the 300 calls of each variant.  tmp_ds.vcf also
# carries a DS field, fractional for every seventh sample.
awk 'BEGIN { OFS = "\t"; split("0 6 45 120", blank_ct, " ") }
     /^##/ { print > "tmp_hard.vcf"; print > "tmp_ds.vcf"; next }
     /^#/ {
       print > "tmp_hard.vcf"
       print "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"dosage\">" > "tmp_ds.vcf"
       print > "tmp_ds.vcf"
       next
     }
     {
       cur_blank_ct = blank_ct[(vidx % 4) + 1]
       line_gt = $1
       line_ds = $1
       for (c = 2; c <= 8; ++c) {
         line_gt = line_gt OFS $c
         line_ds = line_ds OFS $c
       }
       line_gt = line_gt OFS "GT"
       line_ds = line_ds OFS "GT:DS"
       for (s = 10; s <= NF; ++s) {
         gt = $s
         # 37 and 300 are coprime, so this is a permutation of the samples.
         if (((s - 10) * 37 + vidx * 11) % 300 < cur_blank_ct) {
           gt = "./."
           ds = "."
         } else {
           tmp = gt
           ds = gsub(/1/, "1", tmp)
           if ((s - 10) % 7 == 3) {
             ds = (ds == 2)? 1.7 : ds + 0.3
           }
         }
         line_gt = line_gt OFS gt
         line_ds = line_ds OFS gt ":" ds
       }
       print line_gt > "tmp_hard.vcf"
       print line_ds > "tmp_ds.vcf"
       ++vidx
     }' tmp_dense.vcf

$plink2 --vcf tmp_hard.vcf --make-pgen --out tmp_hard > /dev/null
$plink2 --vcf tmp_ds.vcf dosage=DS --make-pgen --out tmp_ds > /dev/null

# Score ALT for even-numbered variants and REF for odd ones.
awk 'BEGIN { OFS = "\t" }
     /^#/ { next }
     { print $3, (vidx % 2)? $4 : $5, ((vidx * 7) % 11 - 5) / 4; ++vidx }' tmp_hard.pvar > tmp_weights.txt

for data in hard ds; do
    $plink2 --pfile tmp_$data --freq --out tmp_$data > /dev/null
    for model in additive dominant recessive; do
        modif=$model
        if [ $model = additive ]; then
            modif=
        fi
        for meanimpute in 1 0; do
            nmi=
            if [ $meanimpute = 0 ]; then
                nmi=no-mean-imputation
            fi
            $plink2 --pfile tmp_$data --read-freq tmp_$data.afreq --score tmp_weights.txt $modif $nmi cols=+scoresums --out tmp_score > /dev/null
            awk -v model=$model -v meanimpute=$meanimpute -f oracle.awk tmp_$data.afreq tmp_weights.txt tmp_$data.vcf > tmp_expected.txt
            awk -v label="$data $model meanimpute=$meanimpute" \
                'function abs(x) { return (x < 0)? -x : x }
                 NR == FNR { expected[$1] = $2; next }
                 FNR == 1 { for (k = 1; k <= NF; ++k) { if ($k == "SCORE1_SUM") { col = k } }; next }
                 {
                   ++n
                   e = expected[$1]
                   tol = 1e-4 * ((abs(e) > 1)? abs(e) : 1)
                   if (abs($col - e) > tol) {
                     failed = 1; print label ": " $1 " SCORE1_SUM " $col ", expected " e
                     exit 1
                   }
                 }
                 END { if ((!failed) && (n != 300)) { print label ": expected 300 samples, got " n + 0; exit 1 } }' tmp_expected.txt tmp_score.sscore
        done
    done
done

# The dosage copy must actually carry dosages, which move the frequencies.
test "$(awk 'NR == FNR { f[$2] = $5; next } FNR > 1 && f[$2] != $5' tmp_hard.afreq tmp_ds.afreq | wc -l)" -gt 30
