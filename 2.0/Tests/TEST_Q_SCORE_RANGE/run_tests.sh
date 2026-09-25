#!/bin/bash

# --q-score-range range names become part of the output filename, so a name
# that doesn't fit must be rejected.  The length check used to compose its
# error message without acting on it, and the name was then copied into the
# fixed-size output filename buffer.

set -exo pipefail

$1/plink2 $2 $3 --dummy 60 30 --seed 1 --make-pgen --out tmp_data
grep -v '^#' tmp_data.pvar | awk '{print $3 "\t" $5 "\t" (NR % 7) / 10}' > tmp_score.txt
grep -v '^#' tmp_data.pvar | awk '{print $3 "\t" (NR % 10) / 10}' > tmp_pvals.txt

# Sanity check: ordinary range names work.
printf 'r1\t0\t0.3\nr2\t0.2\t0.9\n' > tmp_ranges_ok.txt
$1/plink2 $2 $3 --pfile tmp_data --score tmp_score.txt --q-score-range tmp_ranges_ok.txt tmp_pvals.txt --out plink2_ok
test -e plink2_ok.r1.sscore
test -e plink2_ok.r2.sscore

python3 -c "
open('tmp_ranges_long.txt', 'w').write('r' * 5000 + '\t0\t0.3\n')
"
if $1/plink2 $2 $3 --pfile tmp_data --score tmp_score.txt --q-score-range tmp_ranges_long.txt tmp_pvals.txt --out plink2_long 2> tmp_err.txt; then
    echo "expected --q-score-range to reject the overlong range name"
    exit 1
fi
grep -q "Name too long" tmp_err.txt
