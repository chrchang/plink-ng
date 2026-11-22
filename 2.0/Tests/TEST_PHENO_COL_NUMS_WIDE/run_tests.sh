#!/bin/bash

# --pheno-col-nums on a wide headerless phenotype file.
#
# The selected-column count was popcounted over the file's full column count
# instead of the bitarray's size, reading past the end of the bitarray.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
import random
random.seed(1)
ids = ["S%d" % i for i in range(6)]
with open("tmp_pheno.txt", "w") as f:
    for s in ids:
        f.write(s + " " + s + " " + " ".join("%.2f" % random.random() for _ in range(1024)) + "\n")
with open("tmp_data.vcf", "w") as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    f.write("1\t1\tv1\tA\tG\t.\t.\t.\tGT\t" + "\t".join(["0/1"] * 6) + "\n")
PY
$plink2 --vcf tmp_data.vcf --double-id --make-pgen --out tmp_data > /dev/null
# Columns 3-5 are the first three phenotypes; 1026 is the last one.
$plink2 --pfile tmp_data --pheno tmp_pheno.txt --pheno-col-nums 3-5,1026 --make-just-psam --out tmp_out > /dev/null
[[ "$(head -n 1 tmp_out.psam)" == "$(printf '#FID\tIID\tSEX\tPHENO1\tPHENO2\tPHENO3\tPHENO4')" ]] || exit 1
[[ "$(awk 'NR == 2 {print $4, $5, $6, $7}' tmp_out.psam)" == "$(awk 'NR == 1 {print $3, $4, $5, $1026}' tmp_pheno.txt)" ]] || exit 1
