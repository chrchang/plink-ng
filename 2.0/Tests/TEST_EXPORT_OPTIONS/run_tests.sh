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
