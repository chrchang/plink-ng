#!/usr/bin/env python3

# Writes tmp_data.vcf: 200 samples, 300 phased biallelic variants 1 kb apart.
# Every 10th variant (and some others) is independent; the rest copy a
# nearby variant with 8% noise, so LD decays with distance.  Also writes
# tmp_swap.vcf, which adds v5swap right after v5: the same genotypes, but
# with every heterozygous call's phase flipped.  That pair has genotype r^2 1
# but low haplotype r^2.

import random
random.seed(11)
nsample = 200
nvar = 300
haps = []  # per variant: list of 2*nsample alleles
for v in range(nvar):
    if v % 10 == 0 or random.random() < 0.15:
        p = random.uniform(0.15, 0.5)
        h = [1 if random.random() < p else 0 for _ in range(2 * nsample)]
    else:
        # copy a recent variant with some noise, so LD decays with distance
        src = haps[v - random.randint(1, min(v % 10, 4) or 1)]
        h = [(a if random.random() > 0.08 else (1 - a)) for a in src]
    haps.append(h)
with open("tmp_data.vcf", "w") as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=10000000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("S%d" % i for i in range(nsample)) + "\n")
    for v, h in enumerate(haps):
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t" % (1000 * (v + 1), v) + "\t".join("%d|%d" % (h[2*s], h[2*s+1]) for s in range(nsample)) + "\n")
with open("tmp_data.vcf") as f:
    lines = f.read().splitlines()
with open("tmp_swap.vcf", "w") as f:
    for line in lines:
        f.write(line + "\n")
        fields = line.split("\t")
        if len(fields) > 2 and fields[2] == "v5":
            fields[1] = str(int(fields[1]) + 500)
            fields[2] = "v5swap"
            fields[9:] = [gt[::-1] for gt in fields[9:]]
            f.write("\t".join(fields) + "\n")
