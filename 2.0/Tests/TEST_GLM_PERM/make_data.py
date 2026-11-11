# Writes tmp_data.vcf and tmp_data.psam for the --glm permutation tests.
#   chr1: 8 biallelic variants
#   chr2: 4 triallelic variants (two ALT alleles each)
#   chrX: 4 biallelic variants (males haploid)
#   chrY: 3 biallelic variants (males only)
# Phenotypes:
#   QT    quantitative, with effects from several variants so that adaptive
#         permutation stops at a range of permutation counts
#   CC    case/control
#   TWO   quantitative with only two values, so that permutations often
#         reproduce the original statistic exactly
#   CL    one category per sample, for identity "permutations" with
#         --permute-within
# Covariates: COV1, COV2.
import random

random.seed(7)
n = 300
sexes = [1 + (i % 2) for i in range(n)]
variants = []  # (chrom, id, alts, freqs, effects)
for j in range(8):
    variants.append(("1", "a%d" % j, ["G"], [0.15 + 0.04 * j], [(0.35, 0.0, 0.2, 0.0, 0.0, 0.15, 0.0, 0.1)[j]]))
for j in range(4):
    variants.append(("2", "m%d" % j, ["G", "T"], [(0.3, 0.2), (0.25, 0.25), (0.4, 0.1), (0.2, 0.3)][j], [(0.3, 0.0), (0.0, 0.25), (0.2, 0.2), (0.0, 0.3)][j]))
for j in range(4):
    variants.append(("X", "x%d" % j, ["G"], [0.2 + 0.05 * j], [(0.25, 0.0, 0.1, 0.0)[j]]))
for j in range(3):
    variants.append(("Y", "y%d" % j, ["G"], [0.3 + 0.1 * j], [(0.2, 0.0, 0.0)[j]]))

y = [random.gauss(0, 1) for _ in range(n)]
lines = []
for chrom, vid, alts, freqs, effects in variants:
    gts = []
    for i in range(n):
        male = sexes[i] == 1
        ploidy = 1 if (chrom in ("X", "Y") and male) else 2
        if chrom == "Y" and not male:
            gts.append(".")
            continue
        alleles = []
        for _ in range(ploidy):
            u = random.random()
            a = 0
            cum = 0.0
            for k, f in enumerate(freqs):
                cum += f
                if u < cum:
                    a = k + 1
                    break
            alleles.append(a)
        for a in alleles:
            if a:
                y[i] += effects[a - 1] * (2 // ploidy)
        gts.append("/".join(str(a) for a in alleles))
    # outside the chrX pseudoautosomal regions
    pos = 10000000 + 1000 * (len(lines) + 1)
    lines.append("%s\t%d\t%s\tA\t%s\t.\t.\t.\tGT\t%s\n" % (chrom, pos, vid, ",".join(alts), "\t".join(gts)))

with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c in ("1", "2", "X", "Y"):
        f.write("##contig=<ID=%s,length=100000000>\n" % c)
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("s%d" % i for i in range(n)) + "\n")
    f.writelines(lines)

median = sorted(y)[n // 2]
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tQT\tCC\tTWO\tCL\tCOV1\tCOV2\n")
    for i in range(n):
        cc = 2 if y[i] > median else 1
        two = 10 + random.randint(0, 1)
        f.write("s%d\t%d\t%.6f\t%d\t%d\tk%d\t%.6f\t%.6f\n" % (i, sexes[i], y[i], cc, two, i, random.random(), 2 * random.random()))
