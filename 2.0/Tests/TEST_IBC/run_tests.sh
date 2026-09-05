#!/bin/bash

# --ibc, checked against PLINK 1.9.
#
# The two disagree by design on monomorphic variants: PLINK 1.9 keeps them in
# the denominator, plink2 skips them, since the estimators are undefined
# there.  The comparison therefore runs on a --mac 1 filtered fileset, and the
# difference itself is checked separately.

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.02 --simulate-ncases 100 --simulate-ncontrols 100 --out tmp_raw > /dev/null
plink --bfile tmp_raw --mac 1 --make-bed --out tmp_data > /dev/null

compare() {
    awk '
        function abs(x) { return (x < 0)? -x : x }
        function close_enough(a, b) { return abs(a - b) <= 1e-6 + 1e-6 * abs(a) }
        FNR == NR {
            if (FNR > 1) { nobs[$2] = $3; f1[$2] = $4; f2[$2] = $5; f3[$2] = $6; ++n1 }
            next
        }
        /^#/ { next }
        {
            ++n2;
            if (!($2 in nobs)) { print "sample missing from the PLINK 1.9 report: " $2; failed = 1; exit 1 }
            if (nobs[$2] != $3) { print "OBS_CT mismatch on " $2 ": " nobs[$2] " vs " $3; failed = 1; exit 1 }
            if (!close_enough(f1[$2], $4) || !close_enough(f2[$2], $5) || !close_enough(f3[$2], $6)) {
                print "Fhat mismatch on " $2 ": " f1[$2] "/" f2[$2] "/" f3[$2] " vs " $4 "/" $5 "/" $6;
                failed = 1; exit 1
            }
        }
        END {
            if (failed) { exit 1 }
            if (n1 != n2) { print "sample count mismatch: " n1 " vs " n2; exit 1 }
            if (n1 == 0) { print "no samples compared"; exit 1 }
            print n1 " samples matched"
        }
    ' "$1" "$2"
}

# 1. All three estimators, against PLINK 1.9.
plink --bfile tmp_data --ibc --out plink19
$1/plink2 $2 $3 --bfile tmp_data --ibc --out plink2
compare plink19.ibc plink2.ibc

# 2. Same on a sample subset, which changes the imputed allele frequencies.
head -n 120 tmp_data.fam | cut -f1,2 -d' ' > tmp_keep.txt
plink --bfile tmp_data --keep tmp_keep.txt --ibc --out plink19_sub
$1/plink2 $2 $3 --bfile tmp_data --keep tmp_keep.txt --ibc --out plink2_sub
compare plink19_sub.ibc plink2_sub.ibc

# 3. --read-freq: the estimators are frequency-based, so this has to be
#    honored the same way.
plink --bfile tmp_data --freq --out plink19_freq > /dev/null
plink --bfile tmp_data --read-freq plink19_freq.frq --ibc --out plink19_rf
$1/plink2 $2 $3 --bfile tmp_data --read-freq plink19_freq.frq --ibc --out plink2_rf
compare plink19_rf.ibc plink2_rf.ibc

# 4. Monomorphic variants: PLINK 1.9 and GCTA keep them in OBS_CT with a zero
#    contribution, and so does plink2, so the unfiltered fileset has to match
#    too.  simulate.txt's rare block is what produces them.
plink --bfile tmp_raw --freq --out tmp_freq > /dev/null
test "$(awk 'NR > 1 && $5 == 0' tmp_freq.frq | wc -l)" -gt 0
plink --bfile tmp_raw --ibc --out plink19_raw
$1/plink2 $2 $3 --bfile tmp_raw --ibc --out plink2_raw
compare plink19_raw.ibc plink2_raw.ibc

# 5. cols= drops the columns it doesn't name.
$1/plink2 $2 $3 --bfile tmp_data --ibc cols=fhat2 --out plink2_cols
head -n 1 plink2_cols.ibc | grep -qx '#IID	FHAT2'

# 6. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --ibc zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.ibc.zst > plink2_zs.ibc
diff -q plink2.ibc plink2_zs.ibc
