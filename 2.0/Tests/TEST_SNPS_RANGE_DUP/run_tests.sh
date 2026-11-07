#!/bin/bash

# --snps range parsing.
#
# The duplicate check compared a range's end with earlier IDs using the
# range start's length, so "rs115,rs1-rs11" was rejected as a duplicate.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
    ids = ["rs%d" % i for i in range(1, 13)] + ["rs115"]
    for i, vid in enumerate(ids):
        f.write("1\t%d\t%s\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n" % (100 * (i + 1), vid))
PY
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

$plink2 --pfile tmp_data --snps rs115,rs1-rs11 --make-just-pvar --out tmp_snps > /dev/null
[[ $(grep -vc '^#' tmp_snps.pvar) -eq 12 ]] || exit 1
# A genuine duplicate is still an error.
if $plink2 --pfile tmp_data --snps rs3,rs1-rs3 --make-just-pvar --out tmp_dup > /dev/null 2>&1; then
    echo "duplicate range end accepted"
    exit 1
fi
