#!/bin/bash

# --pheno-svd output must match numpy's SVD.  The numpy result is captured in
# ref_svd.txt (see make_ref.py for how to regenerate it), so this test does
# not need numpy itself.
#
# dgesvd's U is column-major with leading dimension = sample count, but it was
# read as if it were row-major with leading dimension = SVD rank.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(9)
n = 25
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tQ1\tQ2\tQ3\n")
    for i in range(n):
        base = random.gauss(0, 1)
        f.write("s%d\t%d\t%.6f\t%.6f\t%.6f\n" % (i, 1 + i % 2, base + random.gauss(0, 0.3), random.gauss(0, 1), -base + random.gauss(0, 0.5)))
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("s%d" % i for i in range(n)) + "\n")
    f.write("1\t100\tv1\tA\tG\t.\t.\t.\tGT\t" + "\t".join(["0/1"] * n) + "\n")
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --pheno-svd 2 --out tmp_svd > /dev/null
python3 check_svd.py tmp_data.psam tmp_svd.svd.pheno tmp_svd.svd.pheno_wts ref_svd.txt
