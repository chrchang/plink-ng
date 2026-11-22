#!/bin/bash

# --mind and --geno thresholds under --y-nosex-missing-stats.
#
# With that flag, unknown-sex samples have chrY missingness counted, so chrY
# must also be in their --mind denominator, and they must be in the chrY
# --geno denominator.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# M1, M2 male; N unknown sex; F1 female.  10 autosomal + 10 chrY variants.
# N is missing 3 autosomal calls and no chrY calls: 3/20 = 0.15 when chrY
# counts, 3/10 = 0.3 when it doesn't.
python3 - <<'PY'
samples = ["M1", "M2", "N", "F1"]
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##contig=<ID=Y,length=100000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    for i in range(10):
        gts = ["0/1", "0/0", "./." if i < 3 else "0/1", "1/1"]
        f.write("1\t%d\ta%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 + i, i, "\t".join(gts)))
    for i in range(10):
        # y0: M2 missing, so 1/3 missing with N in the denominator, 1/2 without.
        gts = ["1", "." if i == 0 else "0", "1", "."]
        f.write("Y\t%d\ty%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 + i, i, "\t".join(gts)))
PY
printf '#IID\tSEX\nM1\t1\nM2\t1\nN\tNA\nF1\t2\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

$plink2 --pfile tmp_data --y-nosex-missing-stats --mind 0.2 --make-just-psam --out tmp_mind > /dev/null
got=$(awk 'NR > 1 {printf "%s ", $1}' tmp_mind.psam)
[[ "$got" == "M1 M2 N F1 " ]] || exit 1
# Without the flag, chrY doesn't count for N, so N is removed.
$plink2 --pfile tmp_data --mind 0.2 --make-just-psam --out tmp_mind_default > /dev/null
got=$(awk 'NR > 1 {printf "%s ", $1}' tmp_mind_default.psam)
[[ "$got" == "M1 M2 F1 " ]] || exit 1

$plink2 --pfile tmp_data --y-nosex-missing-stats --geno 0.4 --make-just-pvar --out tmp_geno > /dev/null
grep -q $'\ty0\t' tmp_geno.pvar
$plink2 --pfile tmp_data --geno 0.4 --make-just-pvar --out tmp_geno_default > /dev/null
! grep -q $'\ty0\t' tmp_geno_default.pvar || exit 1
