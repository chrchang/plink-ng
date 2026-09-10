#!/bin/bash

# --neighbour, checked against PLINK 1.9 where 1.9 is right, and against
# plink2's own output where it is not.

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.03 --simulate-ncases 60 --simulate-ncontrols 60 --out tmp_data > /dev/null

compare() {
    awk -f compare.awk "$1" "$2"
}

# 1. The single nearest neighbour.
plink --bfile tmp_data --neighbour 1 1 --out plink19_1
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 1 --out plink2_1
compare plink19_1.nearest plink2_1.nearest

# 2. A range starting at 1.
plink --bfile tmp_data --neighbour 1 5 --out plink19_15
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 5 --out plink2_15
compare plink19_15.nearest plink2_15.nearest

# 3. A range that does not start at 1, and a single rank above 1.  These are
#    not compared against PLINK 1.9: released 1.9 mislabels the rows whenever
#    the first rank is above 1 (it prints rank 1 upward under the labels n1
#    upward), so the reference would be wrong.  Step 4 checks them against
#    plink2's own --neighbour 1 <n2> instead, which is the same requirement
#    without depending on a fixed 1.9.
$1/plink2 $2 $3 --bfile tmp_data --neighbour 3 7 --out plink2_37
$1/plink2 $2 $3 --bfile tmp_data --neighbour 4 4 --out plink2_44

# 4. --neighbour <n1> <n2> is a slice of --neighbour 1 <n2>.
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 7 --out plink2_17
awk '
    FNR == NR { if (FNR > 1) { full[$1 "\t" $2 "\t" $3] = $4 "\t" $5 "\t" $6 "\t" $7 }; next }
    /^#/ { next }
    {
        k = $1 "\t" $2 "\t" $3;
        v = $4 "\t" $5 "\t" $6 "\t" $7;
        if (!(k in full)) { print "row " k " is not in the 1..7 report"; exit 1 }
        if (full[k] != v) { print "row " k " disagrees with the 1..7 report"; exit 1 }
        ++n
    }
    END { print n " rows agree with the 1..7 report" }
' plink2_17.nearest plink2_37.nearest

# 5. IBS values match "--distance ibs flat-missing".
$1/plink2 $2 $3 --bfile tmp_data --distance ibs flat-missing square --out plink2_mibs
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 3 --out plink2_m3
awk '
    FILENAME ~ /mibs\.id$/ { if (FNR > 1) { id[$2] = FNR - 1 }; next }
    FILENAME ~ /mibs$/ { ++r; for (i = 1; i <= NF; ++i) { v[r "," i] = $i }; next }
    /^#/ { next }
    {
        a = id[$2]; b = id[$7];
        if (a == "" || b == "") { print "unknown sample on " $0; exit 1 }
        d = v[a "," b];
        if (d == "") { print "no matrix cell for " $2 "/" $7; exit 1 }
        diff = d - $4; if (diff < 0) { diff = -diff }
        if (diff > 1e-9) { print "IBS " $4 " does not match the .mibs cell " d " for " $2 "/" $7; exit 1 }
        ++n
    }
    END { print n " IBS values match the .mibs matrix" }
' plink2_mibs.mibs.id plink2_mibs.mibs plink2_m3.nearest

# 6. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 3 zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.nearest.zst > plink2_zs.nearest
diff -q plink2_m3.nearest plink2_zs.nearest

# 7. cols= trims the report, and only the named columns appear.
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 3 cols=id,nn,ibs --out plink2_cols
head -1 plink2_cols.nearest | grep -qx '#IID	NN	IBS'
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 3 cols=maybefid,id,nn,ibs,z,id2 --out plink2_cols2
head -1 plink2_cols2.nearest | grep -qx '#FID	IID	NN	IBS	Z	FID2	IID2'
diff -q plink2_m3.nearest plink2_cols2.nearest

# 8. --parallel is rejected.
fails() {
    "$@" && exit 1 || true
}
fails $1/plink2 $2 $3 --bfile tmp_data --neighbour 1 3 --parallel 1 2 --out plink2_bad

# 9. --neighbour with a --distance that wants the default rescaling is rejected.
fails $1/plink2 $2 $3 --bfile tmp_data --neighbour 1 3 --distance --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_data --neighbour 1 3 --distance ibs flat-missing --out plink2_both
diff -q plink2_m3.nearest plink2_both.nearest
