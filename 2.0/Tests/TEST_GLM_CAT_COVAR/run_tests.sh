#!/bin/bash

# --glm logistic regression with categorical covariates.
#
# Double-precision logistic regression exploits the structure of categorical
# covariates' dummy columns (each sample has at most one dummy "on") when
# computing the Hessian, gradient, and linear predictor.  Results must match
# a run where the same covariates are pre-expanded into 0/1 numeric columns,
# which takes the generic dense path.  Covers several categorical covariates
# (including a 2-level one), a sample with a missing categorical value,
# missing genotype calls, a multiallelic variant, chrX, the genotypic and
# interaction models, 'intercept' (scale-only covariate rescaling), Firth
# fallback, and multiple threads.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import math
import random
random.seed(31)
n = 700
ids = ["s%d" % i for i in range(n)]
sex = [1 + (i % 2) for i in range(n)]
# Distinct category frequencies, so the omitted (most common) category is
# unambiguous.
def draw(levels, weights):
    return random.choices(levels, weights)[0]
batch_levels = ["b%d" % j for j in range(12)]
batch_w = [30 - 2 * j for j in range(12)]
center_levels = ["c%d" % j for j in range(5)]
center_w = [9, 7, 5, 3, 1]
flag_levels = ["yes", "no"]
flag_w = [3, 2]
rows = []
for i in range(n):
    rows.append({
        "PC1": round(random.gauss(0, 1), 5),
        "PC2": round(random.gauss(0, 1), 5),
        "BATCH": draw(batch_levels, batch_w),
        "CENTER": draw(center_levels, center_w),
        "FLAG": draw(flag_levels, flag_w),
    })
rows[17]["BATCH"] = "NONE"  # plink2's default missing-category name
genos = []
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n##contig=<ID=X,length=100000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for chrom, pos0, vct in (("1", 1000, 7), ("X", 10000000, 2)):
        for v in range(vct):
            maf = 0.1 + 0.06 * v
            multi = (chrom == "1" and v == 6)
            gts = []
            for i in range(n):
                if random.random() < 0.03:
                    gts.append("./." if chrom == "1" or sex[i] == 2 else ".")
                    continue
                if chrom == "1" or sex[i] == 2:
                    a1 = int(random.random() < maf)
                    a2 = int(random.random() < maf)
                    if multi:
                        a1 = a1 * (1 + (random.random() < 0.4))
                        a2 = a2 * (1 + (random.random() < 0.4))
                    gts.append("%d/%d" % (min(a1, a2), max(a1, a2)))
                else:
                    gts.append(str(int(random.random() < maf)))
            f.write("%s\t%d\tv%s_%d\tA\t%s\t.\t.\t.\tGT\t%s\n" % (chrom, pos0 + 100 * v, chrom, v, "G,T" if multi else "G", "\t".join(gts)))
batch_eff = {lv: random.gauss(0, 0.5) for lv in batch_levels}
center_eff = {lv: random.gauss(0, 0.5) for lv in center_levels}
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tCC\n")
    for i in range(n):
        r = rows[i]
        z = 0.4 * r["PC1"] - 0.3 * r["PC2"] + center_eff[r["CENTER"]] + (0.3 if r["FLAG"] == "yes" else 0.0)
        if r["BATCH"] != "NONE":
            z += batch_eff[r["BATCH"]]
        cc = 2 if random.random() < 1.0 / (1.0 + math.exp(-z)) else 1
        f.write("%s\t%d\t%d\n" % (ids[i], sex[i], cc))
cats = ("BATCH", "CENTER", "FLAG")
with open("tmp_covar_cat.txt", "w") as f:
    f.write("#IID\tPC1\tBATCH\tPC2\tCENTER\tFLAG\n")
    for i in range(n):
        r = rows[i]
        f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ids[i], r["PC1"], r["BATCH"], r["PC2"], r["CENTER"], r["FLAG"]))
# Pre-expanded version: same column order, one 0/1 column per non-omitted
# category (omitting the most common one among samples plink2 keeps).
kept = [i for i in range(n) if rows[i]["BATCH"] != "NONE"]
dummies = {}
for c, levels in (("BATCH", batch_levels), ("CENTER", center_levels), ("FLAG", flag_levels)):
    counts = {lv: sum(1 for i in kept if rows[i][c] == lv) for lv in levels}
    largest = max(levels, key=lambda lv: counts[lv])
    assert sorted(counts.values()).count(counts[largest]) == 1
    dummies[c] = [lv for lv in levels if lv != largest and counts[lv]]
order = ["PC1", "BATCH", "PC2", "CENTER", "FLAG"]
with open("tmp_covar_num.txt", "w") as f:
    hdr = []
    for c in order:
        hdr.extend([c] if c not in cats else ["%s_%s" % (c, lv) for lv in dummies[c]])
    f.write("#IID\t" + "\t".join(hdr) + "\n")
    for i in range(n):
        r = rows[i]
        vals = []
        for c in order:
            if c not in cats:
                vals.append(str(r[c]))
            else:
                for lv in dummies[c]:
                    vals.append("NA" if r[c] == "NONE" else ("1" if r[c] == lv else "0"))
        f.write(ids[i] + "\t" + "\t".join(vals) + "\n")
with open("tmp_num_names.txt", "w") as f:
    f.write(",".join(hdr) + "\n")
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
num_names=$(cat tmp_num_names.txt)

# run_pair <name> <args...>
run_pair() {
    name=$1
    shift
    $plink2 --pfile tmp_data --covar tmp_covar_cat.txt "$@" --out tmp_cat_$name > /dev/null
    $plink2 --pfile tmp_data --covar tmp_covar_num.txt --covar-name $num_names "$@" --out tmp_num_$name > /dev/null
    found=0
    for f in tmp_num_$name.*.glm.*; do
        suffix=${f#tmp_num_$name.}
        python3 compare.py tmp_cat_$name.$suffix $f
        found=1
    done
    [ $found -eq 1 ]
}

run_pair default --glm cols=+beta
run_pair nofirth --glm no-firth --threads 3
run_pair genotypic --glm genotypic no-firth
run_pair interaction --glm interaction no-firth
run_pair intercept --glm intercept no-firth hide-covar
