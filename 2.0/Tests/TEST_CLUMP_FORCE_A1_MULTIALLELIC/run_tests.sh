#!/bin/bash

# --clump --clump-force-a1 with a multiallelic variant after a biallelic one.
#
# The "A1 is ALT" flag set by a biallelic line wasn't reset for a following
# multiallelic line, so its clump-member entry named the next allele.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(3)
n = 100
ids = ["s%d" % i for i in range(n)]
v1, v2 = [], []
for i in range(n):
    h = [random.random() < 0.4 for _ in range(2)]
    v1.append("%d/%d" % (h[0], h[1]))
    # v2's ALT1 (C) tracks v1's ALT, with an occasional ALT2 (T).
    a = [1 if x else (2 if random.random() < 0.05 else 0) for x in h]
    v2.append("%d/%d" % tuple(sorted(a)))
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    f.write("1\t1000\tv1\tA\tG\t.\t.\t.\tGT\t" + "\t".join(v1) + "\n")
    f.write("1\t2000\tv2\tA\tC,T\t.\t.\t.\tGT\t" + "\t".join(v2) + "\n")
PY
printf '#CHROM\tPOS\tID\tA1\tP\n1\t1000\tv1\tG\t1e-10\n1\t2000\tv2\tC\t1e-6\n' > tmp_assoc.txt
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --clump cols=+sp2 tmp_assoc.txt --clump-force-a1 --out tmp_clump > /dev/null
[[ "$(awk -F '\t' '$3 == "v1" {print $4, $NF}' tmp_clump.clumps)" == "G v2(C)" ]] || exit 1
