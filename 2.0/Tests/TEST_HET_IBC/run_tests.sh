#!/bin/bash

# --het's fhat1/fhat2/fhat3 columns, checked against PLINK 1.9's --ibc.
#
# The two disagree by design on monomorphic variants: PLINK 1.9 keeps them in
# the denominator, while these columns share --het's OBS_CT, which excludes
# them.  The comparison therefore runs on a --mac 1 filtered fileset, and the
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
$1/plink2 $2 $3 --bfile tmp_data --het cols=maybefid,maybesid,nobs,fhat1,fhat2,fhat3 --out plink2
compare plink19.ibc plink2.het

# 2. Same on a sample subset, which changes the imputed allele frequencies.
#    The subset is materialized first, since variants that are monomorphic
#    within it would otherwise trip the OBS_CT difference from section 4.
head -n 120 tmp_data.fam | cut -f1,2 -d' ' > tmp_keep.txt
plink --bfile tmp_data --keep tmp_keep.txt --mac 1 --make-bed --out tmp_sub > /dev/null
plink --bfile tmp_sub --ibc --out plink19_sub
$1/plink2 $2 $3 --bfile tmp_sub --het cols=maybefid,maybesid,nobs,fhat1,fhat2,fhat3 --out plink2_sub
compare plink19_sub.ibc plink2_sub.het

# 3. --read-freq: the estimators are frequency-based, so this has to be
#    honored the same way.
plink --bfile tmp_data --freq --out plink19_freq > /dev/null
plink --bfile tmp_data --read-freq plink19_freq.frq --ibc --out plink19_rf
$1/plink2 $2 $3 --bfile tmp_data --read-freq plink19_freq.frq --het cols=maybefid,maybesid,nobs,fhat1,fhat2,fhat3 --out plink2_rf
compare plink19_rf.ibc plink2_rf.het

# 4. Monomorphic variants: PLINK 1.9 and GCTA keep them in OBS_CT with a zero
#    contribution; these columns share --het's OBS_CT, which excludes them.
#    The unfiltered fileset therefore has to differ, and OBS_CT is where.
plink --bfile tmp_raw --freq --out tmp_freq > /dev/null
test "$(awk 'NR > 1 && $5 == 0' tmp_freq.frq | wc -l)" -gt 0
plink --bfile tmp_raw --ibc --out plink19_raw
$1/plink2 $2 $3 --bfile tmp_raw --het cols=maybefid,maybesid,nobs,fhat1,fhat2,fhat3 --out plink2_raw
if compare plink19_raw.ibc plink2_raw.het > /dev/null 2>&1; then
    echo "monomorphic variants no longer change the result"
    exit 1
fi
# Every sample loses the monomorphic variants it had a call at, so plink2's
# OBS_CT is strictly smaller for all of them.  The exact drop varies by
# sample, since PLINK 1.9 only counts a monomorphic variant for samples that
# are nonmissing there.
awk '
    FNR == NR { if (FNR > 1) {nobs[$2] = $3}; next }
    /^#/ { next }
    {
        if (nobs[$2] - $3 <= 0) {
            print "OBS_CT did not drop on " $2 ": " nobs[$2] " vs " $3; exit 1
        }
        ++n
    }
    END { if (n == 0) { print "no samples compared"; exit 1 } }
' plink19_raw.ibc plink2_raw.het

# 5. cols= drops the columns it doesn't name, and the default --het report is
#    unchanged by any of this.
$1/plink2 $2 $3 --bfile tmp_data --het cols=fhat2 --out plink2_cols
head -n 1 plink2_cols.het | grep -qx '#IID	FHAT2'
$1/plink2 $2 $3 --bfile tmp_data --het --out plink2_default
head -n 1 plink2_default.het | grep -qx '#FID	IID	O(HOM)	E(HOM)	OBS_CT	F'
# The f column is not the same statistic as fhat2, so they must differ.
$1/plink2 $2 $3 --bfile tmp_data --het cols=maybefid,maybesid,f,fhat2 --out plink2_both
awk 'NR > 1 && $3 == $4 {++same} END {if (same == NR - 1) {print "f and fhat2 are identical"; exit 1}}' plink2_both.het

# 6. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --het zs cols=maybefid,maybesid,nobs,fhat1,fhat2,fhat3 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.het.zst > plink2_zs.het
diff -q plink2.het plink2_zs.het

# 7. Multiallelic variants, which PLINK 1.9's --ibc could not handle at all.
#    Checked against the closed forms directly, since there is nothing to
#    compare against.
$1/plink2 $2 $3 --vcf multi.vcf --make-pgen --out tmp_multi
$1/plink2 $2 $3 --pfile tmp_multi --freq --out tmp_multi_freq
$1/plink2 $2 $3 --pfile tmp_multi --het cols=nobs,fhat1,fhat2,fhat3 --out plink2_multi
python3 check_multi.py multi.vcf tmp_multi_freq.afreq plink2_multi.het

# 8. --ibc itself is retired, and says where the statistics went.
if $1/plink2 $2 $3 --bfile tmp_data --ibc --out plink2_bad > /dev/null 2> tmp_err.txt; then
    echo "--ibc unexpectedly succeeded"
    exit 1
fi
tr '\n' ' ' < tmp_err.txt | tr -s ' ' | grep -q -- '--het cols=+fhat1,+fhat2,+fhat3'
