#!/bin/bash

# --fst block-jackknife SE when the block size divides the variant count.
#
# The final block was already saved (and the running sums reset) in that
# case, but was then dropped from the jackknife.  The SE must match a
# weighted block jackknife (Busing et al. 1999) computed from the
# per-variant numerators/denominators.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(4)
n = 60
ids = ["S%d" % i for i in range(n)]
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tPOP\n")
    for i, s in enumerate(ids):
        f.write("%s\t%d\t%s\n" % (s, 1 + i % 2, "A" if i < 30 else "B"))
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for v in range(120):
        pa = 0.1 + 0.8 * random.random()
        pb = min(0.95, max(0.05, pa + random.gauss(0, 0.15)))
        gts = []
        for i in range(n):
            p = pa if i < 30 else pb
            gts.append(("0/0", "0/1", "1/1")[(random.random() < p) + (random.random() < p)])
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(gts)))
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

# 120 variants: blocksize=20 divides evenly, blocksize=25 doesn't.
for bs in 20 25; do
    $plink2 --pfile tmp_data --fst POP blocksize=$bs report-variants vcols=+fstfrac --out tmp_bs$bs > /dev/null
    expected=$(python3 wjack.py tmp_bs$bs.A.B.fst.var $bs)
    got=$(awk 'NR == 2 {print $3, $4}' tmp_bs$bs.fst.summary)
    python3 -c "
import sys
e = [float(x) for x in '$expected'.split()]
g = [float(x) for x in '$got'.split()]
assert all(abs(x - y) <= 1e-5 * abs(x) for x, y in zip(e, g)), (e, g)
"
done
