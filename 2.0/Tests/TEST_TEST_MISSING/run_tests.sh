#!/bin/bash

# --test-missing: differential missingness between cases and controls, tested
# against PLINK 1.9 on autosomes, chrX and chrY.
#
# The chrX and chrY parts matter because PLINK 1.x counts a male heterozygous
# call on a haploid chromosome as missing, and restricts chrY to males.

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.05 --simulate-ncases 400 --simulate-ncontrols 400 --out tmp_data

# Relabel the three 100-variant blocks as chr1, chrX and chrY, and give half
# the samples male sex so that the haploid special cases are exercised.
awk 'BEGIN{OFS="\t"} {if (NR <= 100) $1 = 1; else if (NR <= 200) $1 = 23; else $1 = 24; print}' tmp_data.bim > tmp_relabeled.bim
mv tmp_relabeled.bim tmp_data.bim
awk 'BEGIN{OFS=" "} {$5 = (NR % 2) + 1; print}' tmp_data.fam > tmp_relabeled.fam
mv tmp_relabeled.fam tmp_data.fam

compare() {
    # $1: PLINK 1.9 .missing report, $2: plink2 .missing report.
    awk '
        function abs(x) { return (x < 0)? -x : x }
        function close_enough(a, b) { return abs(a - b) <= 1e-6 + 1e-3 * abs(a) }
        FNR == NR {
            if (FNR > 1) {
                fa[$2] = $3; fu[$2] = $4; p[$2] = $5; ++n1;
            }
            next
        }
        /^#/ { next }
        {
            ++n2;
            if (!($3 in fa)) { print "extra row in plink2 output: " $3; failed = 1; exit 1 }
            if (!close_enough(fa[$3], $4) || !close_enough(fu[$3], $5) || !close_enough(p[$3], $6)) {
                print "mismatch on " $3 ": " fa[$3] "/" fu[$3] "/" p[$3] " vs " $4 "/" $5 "/" $6;
                failed = 1; exit 1
            }
            delete fa[$3];
        }
        END {
            if (failed) { exit 1 }
            if (n1 != n2) { print "row count mismatch: " n1 " vs " n2; exit 1 }
            if (n1 == 0) { print "no rows compared"; exit 1 }
            print n1 " rows matched"
        }
    ' "$1" "$2"
}

# 1. Fisher's exact test.
plink --bfile tmp_data --test-missing --out plink19
$1/plink2 $2 $3 --bfile tmp_data --test-missing --out plink2
compare plink19.missing plink2.missing

# 2. mid-p adjustment.
plink --bfile tmp_data --test-missing midp --out plink19_midp
$1/plink2 $2 $3 --bfile tmp_data --test-missing midp --out plink2_midp
compare plink19_midp.missing plink2_midp.missing

# 3. Same results after a round trip through .pgen, which stores the male
#    heterozygous chrX/chrY calls the same way.
$1/plink2 $2 $3 --bfile tmp_data --make-pgen --out tmp_pgen
$1/plink2 $2 $3 --pfile tmp_pgen --test-missing --out plink2_pgen
diff -q plink2.missing plink2_pgen.missing

# 4. 'zs' output is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --test-missing zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.missing.zst > plink2_zs.missing
diff -q plink2.missing plink2_zs.missing
