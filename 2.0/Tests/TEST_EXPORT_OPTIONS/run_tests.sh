#!/bin/bash

set -exo pipefail

# Covers some less-frequently-used export options.

$1/plink2 $2 $3 --dummy 21 333 0.05 --seed 1 --threads 1 --memory 640 --out tmp_data

# There should be some variants with ALT freq 0, and some with ALT freq 1.
# Confirm that trim-alts affects variants with ALT freq 0, but not those with
# any other ALT freq.
$1/plink2 $2 $3 --pfile tmp_data --make-pgen trim-alts --out trimmed_alts
$1/plink2 $2 $3 --pfile trimmed_alts --freq
if [[ $(cat plink2.afreq | awk '{if (($4 == ".") && ($5 == "0")) print $0}' | wc -c) -eq 0 ]]; then
    exit 1
fi
if [[ $(cat plink2.afreq | awk '{if (($4 == ".") && ($5 != "0")) print $0}' | wc -c) -ne 0 ]]; then
    exit 1
fi
# Confirm that --output-missing-genotype works properly with .bim and .pvar
# output.
$1/plink2 $2 $3 --pfile trimmed_alts --output-missing-genotype 0 --make-bed
if [[ $(cat plink2.bim | awk '{if ($5 == "0") print $0}' | wc -c) -eq 0 ]]; then
    exit 1
fi
if [[ $(cat plink2.bim | awk '{if ($6 == "0") print $0}' | wc -c) -ne 0 ]]; then
    exit 1
fi
$1/plink2 $2 $3 --pfile trimmed_alts --output-missing-genotype 0 --make-pgen
if [[ $(cat plink2.pvar | awk '{if ($5 == "0") print $0}' | wc -c) -eq 0 ]]; then
    exit 1
fi
if [[ $(cat plink2.pvar | awk '{if ($4 == "0") print $0}' | wc -c) -ne 0 ]]; then
    exit 1
fi

# A single-'#' comment line that merely starts with "#CHROM" must not be taken
# for the header line when the INFO column is reloaded for .pvar output.
$1/plink2 $2 $3 --dummy 4 2 --make-pgen --out tmp_ri
printf '##fileformat=VCFv4.3\n#CHROMOSOME names below are GRCh38\n##INFO=<ID=AF,Number=A,Type=Float,Description="af">\n#CHROM\tPOS\tID\tREF\tALT\tINFO\n1\t100\tv1\tA\tG\tAF=0.1\n1\t200\tv2\tC\tT\tAF=0.2\n' > tmp_ri.pvar
$1/plink2 $2 $3 --pfile tmp_ri --make-just-pvar --out plink2_ri
test "$(grep -v '^#' plink2_ri.pvar | cut -f 6 | tr '\n' ' ')" = "AF=0.1 AF=0.2 "

# --export bgz=<level> must change only the compressed size.
$1/plink2 $2 $3 --pfile tmp_data --export vcf bgz --out plink2_bgz
$1/plink2 $2 $3 --pfile tmp_data --export vcf bgz=1 --out plink2_bgz1
$1/plink2 $2 $3 --pfile tmp_data --export vcf bgz=12 --out plink2_bgz12
gzip -dc plink2_bgz.vcf.gz > plink2_bgz.vcf
gzip -dc plink2_bgz1.vcf.gz | diff -q plink2_bgz.vcf -
gzip -dc plink2_bgz12.vcf.gz | diff -q plink2_bgz.vcf -
test $(wc -c < plink2_bgz1.vcf.gz) -gt $(wc -c < plink2_bgz12.vcf.gz)
$1/plink2 $2 $3 --pfile tmp_data --export bcf bgz=1 --out plink2_bgz1b
$1/plink2 $2 $3 --bcf plink2_bgz1b.bcf --export vcf --out plink2_bgz1b
if $1/plink2 $2 $3 --pfile tmp_data --export vcf bgz=0 --out plink2_bgz0; then
    exit 1
fi
if $1/plink2 $2 $3 --pfile tmp_data --export vcf bgz=13 --out plink2_bgz13; then
    exit 1
fi
if $1/plink2 $2 $3 --pfile tmp_data --export A bgz=1 --out plink2_bgza; then
    exit 1
fi
