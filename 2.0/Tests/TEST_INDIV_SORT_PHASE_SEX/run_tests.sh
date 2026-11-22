#!/bin/bash

# --indiv-sort with phased data and with --set-invalid-haploid-missing.
#
# 1. Permuting phase information across a word boundary shifted a word by
#    64 bits (undefined; in practice, by 0), corrupting phase.
# 2. The male/female masks used by --set-invalid-haploid-missing were in
#    input order but applied after the sample permutation.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(17)
n = 150
ids = ["s%03d" % i for i in range(n)]
sexes = [1 if random.random() < 0.5 else 2 for _ in range(n)]
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\n")
    for i in range(n):
        f.write("%s\t%d\n" % (ids[i], sexes[i]))
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=10000000>\n##contig=<ID=X,length=100000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for v in range(6):
        gts = ["%d|%d" % (random.random() < 0.5, random.random() < 0.5) for _ in range(n)]
        f.write("1\t%d\ta%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(gts)))
    for v in range(4):
        gts = []
        for i in range(n):
            if sexes[i] == 1 and random.random() < 0.7:
                gts.append(str(int(random.random() < 0.5)))
            else:
                gts.append("%d|%d" % (random.random() < 0.5, random.random() < 0.5))
        f.write("X\t%d\tx%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (10000000 + 1000 * v, v, "\t".join(gts)))
# Reverse order.
with open("tmp_order.txt", "w") as f:
    for i in reversed(range(n)):
        f.write("%s\n" % ids[i])
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

$plink2 --pfile tmp_data --indiv-sort f tmp_order.txt --make-pgen --out tmp_sorted > /dev/null
$plink2 --pfile tmp_sorted --export vcf --out tmp_sorted > /dev/null
$plink2 --pfile tmp_data --indiv-sort f tmp_order.txt --set-invalid-haploid-missing --make-pgen --out tmp_sorted_shm > /dev/null
$plink2 --pfile tmp_sorted_shm --export vcf --out tmp_sorted_shm > /dev/null
python3 check_sorted.py tmp_data.vcf tmp_data.psam tmp_sorted.vcf tmp_sorted_shm.vcf
