#!/bin/bash

# --glm local covariates, checked against the same values supplied as
# ordinary covariates.
#
# 1. Interaction terms for the ordinary covariates were read from the wrong
#    column (offset by the local covariate count).
# 2. With 1-based 'local-cats=', the last category wasn't omitted, so every
#    category got an indicator column and the model became collinear.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(13)
n = 120
nv = 6
ids = ["s%d" % i for i in range(n)]
g = [[(random.random() < 0.4) + (random.random() < 0.4) for _ in range(n)] for _ in range(nv)]
lc = [random.gauss(0, 1) for _ in range(n)]
cat = [random.choice((1, 2, 3)) for _ in range(n)]
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for v in range(nv):
        # Some missing calls, so each variant has its own sample subset.
        gts = ["./." if (i * 7 + v) % 13 == 0 else ("0/0", "0/1", "1/1")[x] for i, x in enumerate(g[v])]
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(gts)))
# D1/D2 are the category-1/category-2 indicators, shifted by 10 so they
# aren't read as case/control columns (the shift only changes the intercept).
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tQT\tCC\tCOVL\tCOV1\tD1\tD2\n")
    for i in range(n):
        c1 = random.gauss(0, 1)
        qt = 0.4 * g[0][i] + 0.5 * lc[i] + 0.3 * c1 + 0.4 * (cat[i] == 1) + random.gauss(0, 1)
        cc = 2 if (0.5 * g[0][i] + lc[i] + 0.5 * c1 + random.gauss(0, 1)) > 0.7 else 1
        f.write("%s\t%d\t%.6f\t%d\t%.6f\t%.6f\t%d\t%d\n" % (ids[i], 1 + i % 2, qt, cc, lc[i], c1, 10 + (cat[i] == 1), 10 + (cat[i] == 2)))
with open("tmp_local.psam", "w") as f:
    f.write("#IID\tSEX\n")
    for i in range(n):
        f.write("%s\t%d\n" % (ids[i], 1 + i % 2))
with open("tmp_local_real.txt", "w") as f:
    for v in range(nv):
        f.write("\t".join("%.6f" % x for x in lc) + "\n")
with open("tmp_local_cat.txt", "w") as f:
    for v in range(nv):
        f.write("\t".join(str(x) for x in cat) + "\n")
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
common="--pfile tmp_data --pheno-name QT,CC --covar tmp_data.psam"
local_files="local-psam=tmp_local.psam local-pvar=tmp_data.pvar"

# 1. Real-valued local covariate, with interaction terms.
$plink2 $common --covar-name COV1 --glm interaction no-firth local-covar=tmp_local_real.txt $local_files --out tmp_local_real > /dev/null
$plink2 $common --covar-name COVL,COV1 --glm interaction no-firth --out tmp_ordinary_real > /dev/null
for p in QT.glm.linear CC.glm.logistic; do
    python3 compare_glm.py tmp_local_real.$p tmp_ordinary_real.$p LOCAL1=COVL
done

# 2. local-cats=3: categories 1 and 2 get indicators, category 3 is omitted.
$plink2 $common --covar-name COV1 --glm no-firth local-covar=tmp_local_cat.txt $local_files local-cats=3 --out tmp_local_cat > /dev/null
$plink2 $common --covar-name D1,D2,COV1 --glm no-firth --out tmp_ordinary_cat > /dev/null
for p in QT.glm.linear CC.glm.logistic; do
    python3 compare_glm.py tmp_local_cat.$p tmp_ordinary_cat.$p LOCAL1=D1 LOCAL2=D2
done
