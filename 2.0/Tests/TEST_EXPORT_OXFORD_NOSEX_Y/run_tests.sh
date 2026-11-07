#!/bin/bash

# .sample missing-rate column for unknown-sex samples when chrY is present.
#
# Unknown-sex samples have chrY calls exported, so their missing rate must
# use the male (chrY-inclusive) denominator, as males' does.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# 4 autosomal + 4 chrY variants.  N (unknown sex) and M (male) each miss two
# calls, so both have missing rate 2/8; F (female) misses one autosomal call:
# 1/4.
python3 - <<'PY'
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##contig=<ID=Y,length=100000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM\tN\tF\n")
    for i in range(4):
        f.write("1\t%d\ta%d\tA\tG\t.\t.\t.\tGT\t%s\t%s\t%s\n" % (1000 + i, i, "./." if i == 0 else "0/1", "./." if i == 0 else "0/1", "./." if i == 0 else "0/0"))
    for i in range(4):
        f.write("Y\t%d\ty%d\tA\tG\t.\t.\t.\tGT\t%s\t%s\t.\n" % (1000 + i, i, "." if i == 0 else "1", "." if i == 0 else "0"))
PY
printf '#IID\tSEX\nM\t1\nN\tNA\nF\t2\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --export oxford --out tmp_ox > /dev/null
$plink2 --pfile tmp_data --export bgen-1.2 --out tmp_bgen > /dev/null
for f in tmp_ox.sample tmp_bgen.sample; do
    [[ "$(awk 'NR > 2 {printf "%s:%s ", $2, $3}' $f)" == "M:0.25 N:0.25 F:0.25 " ]] || exit 1
done
