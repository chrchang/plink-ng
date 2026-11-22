#!/bin/bash

# --glm logistic regression with categorical covariates and groups of
# mutually exclusive indicator covariates.
#
# Double-precision logistic regression exploits the structure of categorical
# covariates' dummy columns, and of runs of two-valued covariates whose "on"
# (larger) values mark disjoint sample sets (e.g. from --split-cat-pheno),
# when computing the Hessian, gradient, and linear predictor.  Results must
# match a run where the same indicators are negated (coded 0/-1; equivalent
# after internal rescaling, up to sign), which takes the generic dense path.  Covers several categorical covariates
# (including a 2-level one), a sample with a missing categorical value,
# missing genotype calls, a multiallelic variant, chrX and chrY (which have
# their own sample sets), 0/1 and 1/2 indicator coding, --parameters with
# the indicator groups (dropping some
# indicators of a group, whole groups, or the genotype term), the genotypic
# and interaction models, 'intercept' (scale-only covariate
# rescaling), Firth fallback, and multiple threads.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import math
import random
random.seed(31)
n = 700
ids = ["s%d" % i for i in range(n)]
sex = [1 + (i % 2) for i in range(n)]
# Skewed category frequencies, so the omitted (most common) category is
# unambiguous, and the same among males (the chrY sample set) as overall.
def draw(levels, weights):
    return random.choices(levels, weights)[0]
batch_levels = ["b%d" % j for j in range(12)]
batch_w = [80] + [24 - 2 * j for j in range(11)]
center_levels = ["c%d" % j for j in range(5)]
center_w = [20, 7, 5, 3, 2]
flag_levels = ["yes", "no"]
flag_w = [4, 2]
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
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n##contig=<ID=X,length=100000000>\n##contig=<ID=Y,length=50000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for chrom, pos0, vct in (("1", 1000, 7), ("X", 10000000, 2), ("Y", 5000000, 2)):
        for v in range(vct):
            maf = 0.1 + 0.06 * v
            multi = (chrom == "1" and v == 6)
            gts = []
            for i in range(n):
                if chrom == "Y" and sex[i] == 2:
                    gts.append(".")
                    continue
                if random.random() < 0.03:
                    gts.append("./." if chrom == "1" or (chrom == "X" and sex[i] == 2) else ".")
                    continue
                if chrom == "1" or (chrom == "X" and sex[i] == 2):
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
# Pre-expanded versions: same column order, one indicator column per
# non-omitted category (omitting the most common one among samples plink2
# keeps).  tmp_covar_num01.txt and tmp_covar_num12.txt code them 0/1 and 1/2
# (--split-cat-pheno's default), so --glm recognizes each covariate's columns
# as a group of mutually exclusive indicators; tmp_covar_numneg.txt codes
# them 0/-1 (-1 = in category), whose larger value marks overlapping sample
# sets, which forces the dense path.
kept = [i for i in range(n) if rows[i]["BATCH"] != "NONE"]
dummies = {}
for c, levels in (("BATCH", batch_levels), ("CENTER", center_levels), ("FLAG", flag_levels)):
    counts = {lv: sum(1 for i in kept if rows[i][c] == lv) for lv in levels}
    largest = max(levels, key=lambda lv: counts[lv])
    assert sorted(counts.values()).count(counts[largest]) == 1
    male_counts = {lv: sum(1 for i in kept if rows[i][c] == lv and sex[i] == 1) for lv in levels}
    assert max(levels, key=lambda lv: male_counts[lv]) == largest
    assert min(male_counts.values()) > 0
    dummies[c] = [lv for lv in levels if lv != largest and counts[lv]]
order = ["PC1", "BATCH", "PC2", "CENTER", "FLAG"]
hdr = []
for c in order:
    hdr.extend([c] if c not in cats else ["%s_%s" % (c, lv) for lv in dummies[c]])
for off_val, on_val, fn in (("0", "1", "tmp_covar_num01.txt"), ("1", "2", "tmp_covar_num12.txt"), ("0", "-1", "tmp_covar_numneg.txt")):
    with open(fn, "w") as f:
        f.write("#IID\t" + "\t".join(hdr) + "\n")
        for i in range(n):
            r = rows[i]
            vals = []
            for c in order:
                if c not in cats:
                    vals.append(str(r[c]))
                else:
                    for lv in dummies[c]:
                        vals.append("NA" if r[c] == "NONE" else (on_val if r[c] == lv else off_val))
            f.write(ids[i] + "\t" + "\t".join(vals) + "\n")
with open("tmp_num_names.txt", "w") as f:
    f.write(",".join(hdr) + "\n")
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
num_names=$(cat tmp_num_names.txt)

# run_cmp <name> <cat | num01 | num12> <args...>: runs --glm with the given
# covariate file and with tmp_covar_numneg.txt, and compares every output
# file.
run_cmp() {
    name=$1
    kind=$2
    shift 2
    if [ $kind = cat ]; then
        $plink2 --pfile tmp_data --covar tmp_covar_cat.txt "$@" --out tmp_${kind}_$name > /dev/null
    else
        $plink2 --pfile tmp_data --covar tmp_covar_$kind.txt --covar-name $num_names "$@" --out tmp_${kind}_$name > /dev/null
    fi
    $plink2 --pfile tmp_data --covar tmp_covar_numneg.txt --covar-name $num_names "$@" --out tmp_numneg_$name > /dev/null
    found=0
    for f in tmp_numneg_$name.*.glm.*; do
        suffix=${f#tmp_numneg_$name.}
        python3 compare.py tmp_${kind}_$name.$suffix $f
        found=1
    done
    [ $found -eq 1 ]
}

for kind in cat num01 num12; do
    run_cmp default $kind --glm cols=+beta
    run_cmp nofirth $kind --glm no-firth --threads 3
    run_cmp genotypic $kind --glm genotypic no-firth
done
# 'intercept' and 'interaction' prevent centering, so the 1/2 coding isn't
# equivalent there.
for kind in cat num01; do
    run_cmp interaction $kind --glm interaction no-firth
    run_cmp intercept $kind --glm intercept no-firth hide-covar
done

# --parameters can't be used with categorical covariates directly, but can
# with the (detected) indicator groups.  1 = ADD, 2 = PC1, 3-13 = BATCH
# indicators, 14 = PC2, 15-18 = CENTER indicators, 19 = FLAG; with
# 'interaction', 20 = ADDxPC1, 21-31 = ADDxBATCH_*, 32 = ADDxPC2, etc.
# Covariates with no remaining parameter are dropped before regression, so
# the first two cases shrink/remove groups up front.
run_cmp params_partial num01 --parameters 1-2,4,6-9,14-16,19 --glm no-firth
run_cmp params_dropfam num01 --parameters 1-2,14,19 --glm no-firth
run_cmp params_partial12 num12 --parameters 1-2,4,6-9,14-16,19 --glm no-firth
# Here BATCH_b11 is kept for ADDxBATCH_b11 (31) while its main effect (13) is
# excluded, so the regression sees 10 of the BATCH group's 11 indicators; the
# main genotype term is also omitted.
run_cmp params_remap num01 --parameters 2-12,15-20,31 --glm interaction no-firth --threads 2
