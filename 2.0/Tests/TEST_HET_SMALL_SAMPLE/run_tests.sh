#!/bin/bash

# --het small-sample with nonfounders, on sparsely-stored variants.
#
# 'small-sample' makes --het estimate allele frequencies from founders
# itself.  On variants stored as a difflist, the common genotype was
# credited with (founder count - all difflist entries), nonfounder entries
# included.  E(HOM) is checked against a direct computation.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(7)
n_fam = 60
samples = [("f%d" % i, "0", "0", 1 + (i % 2)) for i in range(n_fam)]
samples += [("c%d" % i, "f%d" % (2 * (i // 2)), "f%d" % (2 * (i // 2) + 1), 1 + (i % 2)) for i in range(n_fam)]
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tPAT\tMAT\tSEX\n")
    for s in samples:
        f.write("%s\t%s\t%s\t%d\n" % s)
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=10000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(s[0] for s in samples) + "\n")
    for v in range(200):
        # Mostly rare variants (stored as difflists), plus some common ones.
        maf = 0.03 if v % 10 else 0.3
        gts = [("0/0", "0/1", "1/1")[(random.random() < maf) + (random.random() < maf)] for _ in samples]
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(gts)))
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --het small-sample --out tmp_het > /dev/null
python3 check_het.py tmp_data.vcf tmp_het.het
