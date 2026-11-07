#!/bin/bash

# --ld <variant> <variant> on female-only chrX.
#
# Hardcall phase was discarded on chrX whenever there were no males, so the
# result didn't match the same calls on an autosome.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
samples = ["s%d" % i for i in range(160)]
with open("tmp_female.psam", "w") as f:
    f.write("#IID\tSEX\n")
    for s in samples:
        f.write("%s\t2\n" % s)
def hap_pair():
    # Two loci in partial LD.
    a = random.random() < 0.4
    b = a if random.random() < 0.8 else (random.random() < 0.4)
    return int(a), int(b)
def write_vcf(fn, chrom, phased):
    sep = "|" if phased else "/"
    with open(fn, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=%s,length=10000000>\n" % chrom)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        rows = [[], []]
        for _ in samples:
            h1 = hap_pair()
            h2 = hap_pair()
            for v in range(2):
                rows[v].append("%d%s%d" % (h1[v], sep, h2[v]))
        for v in range(2):
            f.write("%s\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (chrom, 5000000 + 1000 * v, v, "\t".join(rows[v])))
random.seed(12)
write_vcf("tmp_x.vcf", "X", True)
random.seed(12)
write_vcf("tmp_x_as_auto.vcf", "1", True)
PY

ld_line() {
    grep -E '^ *r\^2 = ' "$1".log
}

# Female-only chrX must give the same (phased) result as the same calls on an
# autosome.
$plink2 --vcf tmp_x.vcf --psam tmp_female.psam --make-pgen --out tmp_x > /dev/null
$plink2 --vcf tmp_x_as_auto.vcf --psam tmp_female.psam --make-pgen --out tmp_x_as_auto > /dev/null
$plink2 --pfile tmp_x --ld v0 v1 --out tmp_ld_x > /dev/null
$plink2 --pfile tmp_x_as_auto --ld v0 v1 --out tmp_ld_x_as_auto > /dev/null
[[ "$(ld_line tmp_ld_x)" == "$(ld_line tmp_ld_x_as_auto)" ]]
