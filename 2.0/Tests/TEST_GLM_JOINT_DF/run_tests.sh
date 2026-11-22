#!/bin/bash

# --glm linear joint-test p-values must use the residual degrees of freedom.
#
# The F-test denominator df used to be the sample count, rather than the
# sample count minus the number of predictors.
#
# The expected p-values come from numpy/scipy and are captured in
# ref_joint.txt (see make_ref.py for how to regenerate it), so this test
# does not need numpy or scipy itself.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(21)
n = 40
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("s%d" % i for i in range(n)) + "\n")
    genos = []
    for v in range(3):
        g = [(random.random() < 0.4) + (random.random() < 0.4) for _ in range(n)]
        genos.append(g)
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(("0/0", "0/1", "1/1")[x] for x in g)))
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tQT\tCOV1\tCOV2\n")
    for i in range(n):
        c1 = random.gauss(0, 1)
        c2 = random.gauss(0, 1)
        y = 0.4 * genos[0][i] + 0.3 * (genos[0][i] == 1) + 0.5 * c1 + random.gauss(0, 1)
        f.write("s%d\t%d\t%.6f\t%.6f\t%.6f\n" % (i, 1 + i % 2, y, c1, c2))
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --glm genotypic hide-covar cols=+test --pheno-name QT --covar tmp_data.psam --covar-name COV1,COV2 --covar-variance-standardize --out tmp_glm > /dev/null
python3 check_joint.py tmp_data.vcf tmp_data.psam tmp_glm.QT.glm.linear ref_joint.txt
