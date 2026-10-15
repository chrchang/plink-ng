#!/bin/bash

# --gen import parses genotype probabilities on worker threads.  The output
# must not depend on the thread count, and a malformed line must be reported
# with its own line number even when a later line is also malformed.

set -exo pipefail

$1/plink2 $2 $3 --dummy 300 3000 0.02 dosage-freq=0.2 --seed 1 --export oxford --out tmp_data

for t in 1 8; do
    $1/plink2 $2 $3 --gen tmp_data.gen ref-first --sample tmp_data.sample --threads $t --make-pgen --out plink2_t$t
done
cmp plink2_t1.pgen plink2_t8.pgen
cmp plink2_t1.pvar plink2_t8.pvar

# Line 2000 loses its last token; line 1200 gets an invalid dosage.  The
# earlier one must be reported.
awk 'NR == 2000 { $NF = "" } NR == 1200 { $20 = "0.x5" } { print }' tmp_data.gen > tmp_bad.gen
for t in 1 8; do
    if $1/plink2 $2 $3 --gen tmp_bad.gen ref-first --sample tmp_data.sample --threads $t --make-pgen --out plink2_bad 2> /dev/null; then
        echo "expected --gen to reject tmp_bad.gen"
        exit 1
    fi
    grep -q 'Line 1200 of .gen file has an invalid dosage value' plink2_bad.log
done

awk 'NR == 2000 { $NF = "" } { print }' tmp_data.gen > tmp_bad2.gen
for t in 1 8; do
    if $1/plink2 $2 $3 --gen tmp_bad2.gen ref-first --sample tmp_data.sample --threads $t --make-pgen --out plink2_bad2 2> /dev/null; then
        echo "expected --gen to reject tmp_bad2.gen"
        exit 1
    fi
    grep -q 'Line 2000 of .gen file has fewer tokens than expected' plink2_bad2.log
done
