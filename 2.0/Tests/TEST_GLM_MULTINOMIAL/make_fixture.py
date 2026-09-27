#!/usr/bin/env python3
# Generates the committed fixture for TEST_GLM_MULTINOMIAL (needs numpy; not
# run by the test itself):
#   mn.vcf.gz  37 variants: 31 on chr1, 6 on chrX (haploid male calls)
#   sex.txt    sample sexes
#   pheno.txt  CAT: 4 levels (alpha, beta, delta, gamma)
#              SINGLE: CAT with one sample moved to a level of its own (zeta)
#              BIN: case/control, for the check that other phenotypes are
#              unaffected by the 'multinomial' modifier
#   covar.txt  3 quantitative covariates, 2 binary ones
# Variant roles on chr1:
#   v1-v21  ordinary; about a third have an effect on some levels
#   v22-v24 low frequency (1-3%)
#   v25-v27 missing calls
#   v28-v30 dosages
#   sep1    no ALT allele in level delta (quasi-complete separation)
import gzip
import numpy as np

rng = np.random.default_rng(518)
n = 600
levels = ["alpha", "beta", "delta", "gamma"]
K = len(levels)
sex = rng.integers(1, 3, n)  # 1 = male, 2 = female
C = rng.normal(size=(n, 3)) * [1.0, 2.0, 0.5] + [0.0, 1.0, 0.0]
B = rng.binomial(1, [0.4, 0.2], size=(n, 2))
Z = np.column_stack([C, B])
Zs = (Z - Z.mean(0)) / Z.std(0)
coef_z = rng.normal(scale=0.3, size=(5, K - 1))
alpha = np.array([0.4, -0.3, 0.2])

def draw_geno(maf, male_haploid=False):
    g = rng.binomial(2, maf, n).astype(float)
    if male_haploid:
        male = sex == 1
        g[male] = 2 * rng.binomial(1, maf, male.sum())
    return g

var_ids = []
genos = []
chroms = []
for v in range(1, 31):
    maf = rng.uniform(0.01, 0.03) if 22 <= v <= 24 else rng.uniform(0.08, 0.5)
    var_ids.append("v%d" % v); genos.append(draw_geno(maf)); chroms.append("1")
for v in range(1, 7):
    var_ids.append("x%d" % v); genos.append(draw_geno(rng.uniform(0.1, 0.5), True)); chroms.append("X")
G = np.array(genos)
effects = np.zeros((len(var_ids), K - 1))
for v in range(21):
    if v % 3 == 0:
        effects[v] = rng.normal(scale=0.5, size=K - 1)
effects[var_ids.index("x2")] = [0.4, -0.3, 0.0]
eta = np.column_stack([np.zeros(n), alpha + Zs @ coef_z + ((G - G.mean(1, keepdims=True)).T @ effects)])
P = np.exp(eta - eta.max(1, keepdims=True))
P /= P.sum(1, keepdims=True)
y = np.array([rng.choice(K, p=p) for p in P])

# sep1: ALT alleles only outside level delta
sep = np.where(y == 2, 0.0, rng.binomial(2, 0.2, n).astype(float))
var_ids.insert(30, "sep1")
G = np.vstack([G[:30], sep, G[30:]])
chroms.insert(30, "1")

missing = {v: set(rng.choice(n, size=s, replace=False)) for v, s in (("v25", 5), ("v26", 20), ("v27", 40))}
dosage_vars = {"v28", "v29", "v30"}

gt = {0.0: "0/0", 1.0: "0/1", 2.0: "1/1"}
hap = {0.0: "0", 2.0: "1"}
vcf_lines = []
vcf_lines.append("##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=X>\n")
vcf_lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
vcf_lines.append('##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">\n')
vcf_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("i%d" % i for i in range(n)) + "\n")
for vidx, vid in enumerate(var_ids):
    chrom = chroms[vidx]
    is_ds = vid in dosage_vars
    cells = []
    for i in range(n):
        gv = G[vidx, i]
        if i in missing.get(vid, ()):
            cells.append("./.")
        elif chrom == "X" and sex[i] == 1:
            cells.append(hap[gv])
        elif is_ds:
            ds = min(2.0, max(0.0, gv + rng.normal(scale=0.2)))
            cells.append("%s:%.3f" % (gt[gv], ds))
        else:
            cells.append(gt[gv])
    # chrX positions lie past PAR1
    pos = 1000 * (vidx + 1) + (5000000 if chrom == "X" else 0)
    vcf_lines.append("%s\t%d\t%s\tA\tG\t.\t.\t.\t%s\t%s\n" % (chrom, pos, vid, "GT:DS" if is_ds else "GT", "\t".join(cells)))
# mtime=0 keeps the committed file byte-identical across regenerations
with gzip.GzipFile("mn.vcf.gz", "wb", compresslevel=9, mtime=0) as f:
    f.write("".join(vcf_lines).encode())

single = [levels[k] for k in y]
single[int(np.flatnonzero(y == 0)[0])] = "zeta"
bin_pheno = 1 + (y >= 2)
with open("pheno.txt", "w") as f:
    f.write("#IID\tCAT\tSINGLE\tBIN\n")
    for i in range(n):
        f.write("i%d\t%s\t%s\t%d\n" % (i, levels[y[i]], single[i], bin_pheno[i]))
with open("covar.txt", "w") as f:
    f.write("#IID\tC1\tC2\tC3\tB1\tB2\n")
    for i in range(n):
        f.write("i%d\t%.4f\t%.4f\t%.4f\t%d\t%d\n" % (i, C[i, 0], C[i, 1], C[i, 2], B[i, 0], B[i, 1]))
with open("sex.txt", "w") as f:
    f.write("#IID\tSEX\n")
    for i in range(n):
        f.write("i%d\t%d\n" % (i, sex[i]))
print("level counts", np.bincount(y, minlength=K))
