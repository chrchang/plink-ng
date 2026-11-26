#!/bin/bash

# --glm internal covariate rescaling.
#
# YOB (a year-of-birth-like covariate, mean ~1975, sd ~14) and MID (mean 1e4,
# sd 1) used to make --glm error out with "genotype/covariate scales vary too
# widely".  --glm now always rescales the covariates internally (centering
# them, unless 'intercept' or 'interaction' is specified) and translates
# covariate coefficients back to original units.  Results must match runs with the
# covariates pre-shifted/scaled to a safe range:
#   YOB_REF = (YOB - 1975) / 100, MID_REF = MID - 10000
# or, for 'intercept'/'interaction' runs (where shifting would change the
# intercept and the ADD main effect):
#   YOB_REF = YOB / 100, W_REF = W / 100.
# W (mean 0, sd 500) is the covariate used to check interaction terms, since
# an uncentered YOB makes ADD and ADDxYOB too collinear to report.
# (Genotype terms, P-values, and joint tests are unchanged; YOB coefficients
# are 1/100 of the YOB_REF ones; MID coefficients are identical.)

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(23)
n = 400
ids = ["s%d" % i for i in range(n)]
sex = [1 + (i % 2) for i in range(n)]
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n##contig=<ID=X,length=100000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for chrom, pos0 in (("1", 1000), ("X", 10000000)):
        for v in range(4):
            gts = []
            maf = 0.15 + 0.1 * v
            for i in range(n):
                if random.random() < 0.03:
                    gts.append("./." if chrom == "1" or sex[i] == 2 else ".")
                elif chrom == "1" or sex[i] == 2:
                    gts.append(("0/0", "0/1", "1/1")[(random.random() < maf) + (random.random() < maf)])
                else:
                    gts.append(str(int(random.random() < maf)))
            alt = "G,T" if (chrom == "1" and v == 3) else "G"
            if alt == "G,T":
                gts = [g.replace("1/1", "1/2") if random.random() < 0.3 else g for g in gts]
            f.write("%s\t%d\tv%s_%d\tA\t%s\t.\t.\t.\tGT\t%s\n" % (chrom, pos0 + v * 100, chrom, v, alt, "\t".join(gts)))
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tQT1\tQT2\tCC\tMN\n")
    yobs = []
    mids = []
    cats = []
    ws = []
    for i in range(n):
        yob = random.randint(1950, 2000)
        mid = 10000 + round(random.gauss(0, 1), 4)
        cat = random.choice("abc")
        w = round(random.gauss(0, 500), 1)
        ws.append(w)
        yobs.append(yob)
        mids.append(mid)
        cats.append(cat)
        lin = 0.03 * (yob - 1975) + 0.4 * (mid - 10000) + (0.3 if cat == "b" else 0.0) + 0.001 * w
        qt1 = lin + random.gauss(0, 1)
        qt2 = -lin + random.gauss(0, 1)
        cc = 2 if lin + random.gauss(0, 1) > 0 else 1
        mn = "xyz"[min(2, max(0, int(lin + random.gauss(0, 1) + 1)))]
        f.write("%s\t%d\t%.6f\t%.6f\t%d\t%s\n" % (ids[i], sex[i], qt1, qt2, cc, mn))
with open("tmp_covar.txt", "w") as f:
    f.write("#IID\tYOB\tMID\tCAT\tW\n")
    for i in range(n):
        f.write("%s\t%d\t%.4f\t%s\t%.1f\n" % (ids[i], yobs[i], mids[i], cats[i], ws[i]))
with open("tmp_covar_ref.txt", "w") as f:
    f.write("#IID\tYOB\tMID\tCAT\n")
    for i in range(n):
        f.write("%s\t%.2f\t%.4f\t%s\n" % (ids[i], (yobs[i] - 1975) / 100.0, mids[i] - 10000, cats[i]))
with open("tmp_covar_ref_nc.txt", "w") as f:
    f.write("#IID\tYOB\tCAT\tW\n")
    for i in range(n):
        f.write("%s\t%.2f\t%s\t%.3f\n" % (ids[i], yobs[i] / 100.0, cats[i], ws[i] / 100.0))
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
common="--pfile tmp_data"
all3="--covar-name YOB,MID,CAT"
mults="YOB=0.01 W=0.01"

# Sanity check: the original covariates used to be rejected.  (master would
# error out here.)
$plink2 $common --covar tmp_covar.txt $all3 --pheno-name QT1 --glm --out tmp_auto_chk > /dev/null

# run_pair <name> <reference covariate file> <args...>: runs --glm with the original and reference
# covariates, then compares every output file.
run_pair() {
    name=$1
    ref=$2
    shift 2
    $plink2 $common --covar tmp_covar.txt "$@" --out tmp_auto_$name > /dev/null
    $plink2 $common --covar $ref "$@" --out tmp_ref_$name > /dev/null
    found=0
    for f in tmp_ref_$name.*.glm.*; do
        suffix=${f#tmp_ref_$name.}
        python3 compare_rescaled.py tmp_auto_$name.$suffix $f $mults
        found=1
    done
    [ $found -eq 1 ]
}

# Centering mode: linear (single phenotype, and the QT batch path), logistic,
# Firth, and multinomial, with covariate rows reported.  Includes chrX (with
# its automatic SEX covariate), a multiallelic variant, and missing calls.
run_pair lin tmp_covar_ref.txt $all3 --pheno-name QT1 --glm
run_pair batch tmp_covar_ref.txt $all3 --pheno-name QT1,QT2 --glm
run_pair logistic tmp_covar_ref.txt $all3 --pheno-name CC --glm no-firth
run_pair firth tmp_covar_ref.txt $all3 --pheno-name CC --glm firth
run_pair mn tmp_covar_ref.txt $all3 --pheno-name MN --glm multinomial
# --tests joint test of covariates (no categorical covariates allowed there),
# and hide-covar
run_pair tests tmp_covar_ref.txt --covar-name YOB,MID --pheno-name QT1,CC --tests 2,3 --glm no-firth
run_pair hide tmp_covar_ref.txt $all3 --pheno-name QT1,CC --glm hide-covar no-firth

# Scale-only mode ('intercept' and/or 'interaction').  YOB and W can be fixed
# by scaling alone, so these must work, with INTERCEPT and the ADD main effect
# unchanged.
run_pair intercept tmp_covar_ref_nc.txt --pheno-name QT1,CC --covar-name YOB,CAT,W --glm intercept no-firth
run_pair interact tmp_covar_ref_nc.txt --pheno-name QT1,CC --covar-name W,CAT --glm interaction intercept no-firth
run_pair interact_firth tmp_covar_ref_nc.txt --pheno-name CC --covar-name W --glm interaction firth
# make sure the interaction comparisons weren't vacuous
[ $(awk -F '\t' '$10 == "ADDxW" && $NF == "." {n++} END {print n+0}' tmp_auto_interact.QT1.glm.linear) -ge 6 ]
[ $(awk -F '\t' '$10 == "ADDxW" && $NF == "." {n++} END {print n+0}' tmp_auto_interact_firth.CC.glm.firth) -ge 6 ]

# MID can't be fixed without centering, so 'interaction' must still error out.
if $plink2 $common --covar tmp_covar.txt --covar-name YOB,MID --pheno-name QT1 --glm interaction --out tmp_fail > /dev/null 2>&1; then
    exit 1
fi
grep -q 'automatic covariate rescaling' tmp_fail.log
