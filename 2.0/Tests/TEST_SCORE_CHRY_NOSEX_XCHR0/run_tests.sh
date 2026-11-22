#!/bin/bash

# --score on chrY with unknown-sex samples, and with --xchr-model 0.
#
# 1. The mask used to zero out females' chrY dosages included unknown-sex
#    samples, though those samples' chrY calls count everywhere else.
# 2. --xchr-model 0 removed chrX variants without updating the variant count.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(9)
ids = ["M1", "N1", "F1", "M2", "F2", "N2"]
sex = [1, 0, 2, 1, 2, 0]
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\n")
    for i, s in zip(ids, sex):
        f.write("%s\t%s\n" % (i, "NA" if s == 0 else s))
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##contig=<ID=X,length=100000000>\n##contig=<ID=Y,length=50000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for v in range(5):
        f.write("1\t%d\ta%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 + v, v, "\t".join(random.choice(("0/0", "0/1", "1/1")) for _ in ids)))
    for v in range(3):
        gts = [random.choice(("0", "1")) if s == 1 else random.choice(("0/0", "0/1", "1/1")) for s in sex]
        f.write("X\t%d\tx%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (10000000 + v, v, "\t".join(gts)))
    # M1 and N1 carry the ALT allele.
    f.write("Y\t5000000\ty0\tA\tG\t.\t.\t.\tGT\t1\t1\t.\t0\t.\t0\n")
with open("tmp_weights.txt", "w") as f:
    for v in range(5):
        f.write("a%d\tG\t%g\n" % (v, 0.5 + v))
    for v in range(3):
        f.write("x%d\tG\t%g\n" % (v, 1.5 + v))
    f.write("y0\tG\t2\n")
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

# 1. N1's chrY ALT call scores like M1's.
$plink2 --pfile tmp_data --chr Y --score tmp_weights.txt no-mean-imputation cols=+scoresums --out tmp_y > /dev/null
[[ "$(awk -F '\t' '$1 == "M1" || $1 == "N1" {printf "%s:%s ", $1, $NF}' tmp_y.sscore)" == "M1:2 N1:2 " ]] || exit 1

# 2. --xchr-model 0 must match leaving chrX out entirely.
$plink2 --pfile tmp_data --score tmp_weights.txt no-mean-imputation cols=+scoresums --xchr-model 0 --out tmp_x0 > /dev/null
$plink2 --pfile tmp_data --not-chr X --score tmp_weights.txt no-mean-imputation cols=+scoresums --out tmp_no_x > /dev/null
cmp tmp_x0.sscore tmp_no_x.sscore
