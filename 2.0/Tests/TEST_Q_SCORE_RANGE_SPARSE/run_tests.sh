#!/bin/bash

# --q-score-range with several ranges, on sparsely-stored variants.
#
# The sparse-variant correction was added to ddosage_sums as if it were
# sample-major, when it's range-major, so ranges got each other's
# corrections.  Each range's sums must match a plain --score run restricted
# to that range's variants.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(3)
n = 50
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("s%d" % i for i in range(n)) + "\n")
    for v in range(60):
        # Mostly hom-ALT, so the common genotype is nonzero.
        p = 0.95 if v % 3 else 0.5
        gts = [("0/0", "0/1", "1/1")[(random.random() < p) + (random.random() < p)] for _ in range(n)]
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(gts)))
with open("tmp_score.txt", "w") as f:
    for v in range(60):
        f.write("v%d\tG\t%g\n" % (v, 0.5 + (v % 7) * 0.25))
with open("tmp_pvals.txt", "w") as f:
    for v in range(60):
        f.write("v%d\t%g\n" % (v, (v % 4) / 4.0 + 0.1))
with open("tmp_ranges.txt", "w") as f:
    f.write("r1 0 0.3\nr2 0.3 0.6\nr3 0.6 1\n")
for name, lo, hi in (("r1", 0, 0.3), ("r2", 0.3, 0.6), ("r3", 0.6, 1)):
    with open("tmp_extract_%s.txt" % name, "w") as f:
        for v in range(60):
            p = (v % 4) / 4.0 + 0.1
            if lo <= p <= hi:
                f.write("v%d\n" % v)
PY
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

$plink2 --pfile tmp_data --score tmp_score.txt cols=+scoresums --q-score-range tmp_ranges.txt tmp_pvals.txt --out tmp_q > /dev/null
for r in r1 r2 r3; do
    $plink2 --pfile tmp_data --extract tmp_extract_$r.txt --score tmp_score.txt cols=+scoresums --out tmp_single_$r > /dev/null
    cmp tmp_q.$r.sscore tmp_single_$r.sscore
done
