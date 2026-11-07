#!/bin/bash

# --glm --parameters/--tests on chrY when a covariate is dropped there.
#
# C1 only varies among females, so it's constant (and dropped) on chrY.  The
# chrY parameter and joint-test subsets used to be collapsed against the
# genome-wide covariate set instead of the chrY one, so --tests picked the
# wrong predictors on chrY.  (chrX goes through the same code, with its own
# covariate set.)  Results must match a run where C1 was never specified.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(19)
n = 150
ids = ["s%d" % i for i in range(n)]
sex = [(1, 2, 0)[i % 3] for i in range(n)]
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n##contig=<ID=X,length=100000000>\n##contig=<ID=Y,length=50000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for chrom, pos0 in (("1", 1000), ("X", 10000000), ("Y", 5000000)):
        for v in range(3):
            gts = []
            for i in range(n):
                if chrom == "1" or (chrom == "X" and sex[i] == 2):
                    gts.append(("0/0", "0/1", "1/1")[(random.random() < 0.4) + (random.random() < 0.4)])
                elif chrom == "X" or sex[i] != 2:
                    gts.append(str(int(random.random() < 0.4)))
                else:
                    gts.append(".")
            f.write("%s\t%d\tc%s_%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (chrom, pos0 + v * 100, chrom, v, "\t".join(gts)))
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tQT\tCC\tC1\tC2\tC3\n")
    for i in range(n):
        c1 = random.gauss(0, 1) if sex[i] == 2 else 0.0
        c2 = random.gauss(0, 1)
        c3 = random.gauss(0, 1)
        qt = c2 + 0.5 * c1 + 0.5 * c3 + random.gauss(0, 1)
        cc = 2 if qt + random.gauss(0, 1) > 0 else 1
        f.write("%s\t%s\t%.6f\t%d\t%.6f\t%.6f\t%.6f\n" % (ids[i], "NA" if sex[i] == 0 else sex[i], qt, cc, c1, c2, c3))
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
common="--pfile tmp_data --covar tmp_data.psam"

# --parameters: ADD, C2, C3 (also drops chrX's automatic SEX covariate).
$plink2 $common --pheno-name QT,CC --covar-name C1,C2,C3 --parameters 1,3,4 --glm no-firth --out tmp_param_c1 > /dev/null
$plink2 $common --pheno-name QT,CC --covar-name C2,C3 --parameters 1-3 --glm no-firth --out tmp_param_ref > /dev/null
for p in QT.glm.linear CC.glm.logistic; do
    python3 compare_glm.py tmp_param_c1.$p tmp_param_ref.$p
done

# --tests: joint test of C2 and C3.
$plink2 $common --pheno-name QT --covar-name C1,C2,C3 --tests 3,4 --glm no-firth --out tmp_tests_c1 > /dev/null
$plink2 $common --pheno-name QT --covar-name C2,C3 --tests 2,3 --glm no-firth --out tmp_tests_ref > /dev/null
y_tests() {
    awk -F '\t' '$3 ~ /^cY_/ && $10 ~ /^USER_/ {print $3, $10, $11, $14, $15}' "$1"
}
[[ -n "$(y_tests tmp_tests_c1.QT.glm.linear)" ]] || exit 1
[[ "$(y_tests tmp_tests_c1.QT.glm.linear)" == "$(y_tests tmp_tests_ref.QT.glm.linear)" ]] || exit 1
