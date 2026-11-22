#!/bin/bash

# --make-king-table counts on sparsely-handled variants.
#
# Variants with very few non-common calls take a sparse path.  For a pair of
# rare calls where one sample was missing and the other was the non-common
# homozygote, a "hom-hom" observation was added anyway, so NSNP came out too
# high.  All counts are checked against a brute-force count, both for the
# full table and with --king-table-require-xor.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(6)
# The sparse path needs <= n/33 non-common calls per variant.
n = 330
ids = ["s%d" % i for i in range(n)]
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for v in range(300):
        if v % 3:
            common, other = ("0/0", "1/1") if v % 2 else ("1/1", "0/0")
            gts = [common] * n
            for i in random.sample(range(n), 9):
                gts[i] = random.choice(("./.", "./.", "0/1", other))
        else:
            gts = []
            for i in range(n):
                if random.random() < 0.03:
                    gts.append("./.")
                else:
                    gts.append(("0/0", "0/1", "1/1")[(random.random() < 0.4) + (random.random() < 0.4)])
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(gts)))
with open("tmp_require.txt", "w") as f:
    for i in range(0, n, 11):
        f.write("s%d\n" % i)
PY
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
cols="counts cols=+nsnp,+hethet,+ibs0,+ibs1"
$plink2 --pfile tmp_data --make-king-table $cols --out tmp_all > /dev/null
python3 check_king.py tmp_data.vcf tmp_all.kin0
$plink2 --pfile tmp_data --make-king-table $cols --king-table-require-xor tmp_require.txt --out tmp_xor > /dev/null
python3 check_king.py tmp_data.vcf tmp_xor.kin0 tmp_require.txt
