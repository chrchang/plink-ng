#!/usr/bin/env python3

# Writes tmp_data.vcf and tmp_data.psam: 60 samples, 6 biallelic variants.
# QT depends on v1 and C1; CC depends on v1 and C4.  C2 is nearly collinear
# with C1, and C5 is nearly collinear with v2's genotypes.

import random

random.seed(5)
n = 60
nv = 6
ids = [f"S{i + 1}" for i in range(n)]
geno = [[random.choice([0, 0, 1, 1, 2]) for _ in range(n)] for _ in range(nv)]
gt_strs = ["0/0", "0/1", "1/1"]
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for v in range(nv):
        f.write(f"1\t{1000 * (v + 1)}\tv{v + 1}\tA\tG\t.\t.\t.\tGT\t" + "\t".join(gt_strs[g] for g in geno[v]) + "\n")
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tQT\tCC\tC1\tC2\tC4\tC5\n")
    for i in range(n):
        c1 = random.gauss(0, 1)
        c2 = c1 + random.gauss(0, 0.1)
        c4 = random.gauss(0, 1)
        c5 = geno[1][i] + random.gauss(0, 0.07)
        qt = 0.5 * geno[0][i] + 0.3 * c1 + random.gauss(0, 1)
        cc = 2 if (geno[0][i] + c4 + random.gauss(0, 1)) > 1 else 1
        f.write(f"{ids[i]}\t{1 + (i % 2)}\t{qt:.4f}\t{cc}\t{c1:.4f}\t{c2:.4f}\t{c4:.4f}\t{c5:.4f}\n")
