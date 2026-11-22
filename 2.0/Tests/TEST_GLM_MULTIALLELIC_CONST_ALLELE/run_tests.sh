#!/bin/bash

# --glm on a multiallelic variant with an allele that's constant (absent)
# among the analyzed samples.
#
# Such an allele is dropped from the regression, so the remaining ALT
# allele's result must match the biallelic version of the variant.
# Linear regression used to compute the t-test df as if the dropped allele
# were still a predictor; logistic regression read the wrong coefficient
# slots.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(31)
n = 80
ids = ["s%d" % i for i in range(n)]
g = [(random.random() < 0.35) + (random.random() < 0.35) for _ in range(n)]
with open("tmp_pheno.txt", "w") as f:
    f.write("#IID\tQT\tCC\tCOV1\n")
    for i in range(n):
        c1 = random.gauss(0, 1)
        qt = 0.5 * g[i] + 0.3 * c1 + random.gauss(0, 1)
        cc = 2 if (0.7 * g[i] + c1 + random.gauss(0, 1)) > 0.8 else 1
        if i >= n - 3:
            # These samples carry ALT2 and have missing phenotypes.
            f.write("%s\tNA\tNA\t%.6f\n" % (ids[i], c1))
        else:
            f.write("%s\t%.6f\t%d\t%.6f\n" % (ids[i], qt, cc, c1))
def write(fn, multi):
    with open(fn, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
        gts = []
        for i in range(n):
            if i >= n - 3:
                gts.append("0/2" if multi else "./.")
            else:
                gts.append(("0/0", "0/1", "1/1")[g[i]])
        f.write("1\t1000\tv1\tA\t%s\t.\t.\t.\tGT\t%s\n" % ("G,T" if multi else "G", "\t".join(gts)))
write("tmp_multi.vcf", True)
write("tmp_bi.vcf", False)
PY
for f in multi bi; do
    $plink2 --vcf tmp_$f.vcf --make-pgen --out tmp_$f > /dev/null
    $plink2 --pfile tmp_$f --pheno tmp_pheno.txt --pheno-name QT,CC --covar tmp_pheno.txt --covar-name COV1 --glm hide-covar no-firth cols=+a1count,+machr2 --out tmp_glm_$f > /dev/null
done
# Compare the A1=G rows (all columns from A1 on, except ALT/MACH_R2-adjacent
# allele lists).
for p in QT.glm.linear CC.glm.logistic; do
    a=$(awk -F '\t' 'NR == 1 {for (i = 1; i <= NF; ++i) {if ($i == "A1") {a1 = i}; if ($i == "TEST") {t = i}}; next} $a1 == "G" {out = ""; for (i = t; i <= NF; ++i) {out = out "\t" $i}; print out}' tmp_glm_multi.$p)
    b=$(awk -F '\t' 'NR == 1 {for (i = 1; i <= NF; ++i) {if ($i == "A1") {a1 = i}; if ($i == "TEST") {t = i}}; next} $a1 == "G" {out = ""; for (i = t; i <= NF; ++i) {out = out "\t" $i}; print out}' tmp_glm_bi.$p)
    if [[ -z "$a" || "$a" != "$b" ]]; then
        echo "$p mismatch:"
        echo "$a"
        echo "$b"
        exit 1
    fi
done
