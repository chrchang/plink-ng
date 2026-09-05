#!/bin/bash

# --distance, --distance-matrix and --ibs-matrix, checked against PLINK 1.9.
#
# The two programs agree exactly on 'flat-missing', which is the mode where
# neither weights variants.  The default frequency-weighted mode agrees to
# within PLINK 1.9's weight quantization (it rounds the per-variant weights to
# integers summing to just under 2^32; plink2 keeps them as doubles), so those
# comparisons carry a tolerance.  See section 6 for the one case where the two
# disagree on purpose.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

plink --simulate simulate.txt --simulate-missing 0.06 --simulate-ncases 60 --simulate-ncontrols 60 --out tmp_raw > /dev/null
# PLINK 1.9 gives monomorphic variants a weight of 1 rather than 0 (see 6), so
# the comparison fileset excludes them.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_raw --mac 1 --make-bed --out tmp_data

# Compare two whitespace-delimited matrices cell by cell.  $3 is the maximum
# relative difference allowed.
close() {
    awk -v tol="$3" '
        function abs(x) { return (x < 0)? -x : x }
        FNR == NR { for (i = 1; i <= NF; ++i) {a[FNR, i] = $i}; rows = FNR; cols[FNR] = NF; next }
        {
            if (FNR > rows) { print "extra row " FNR; exit 1 }
            if (NF != cols[FNR]) { print "row " FNR ": " cols[FNR] " vs " NF " columns"; exit 1 }
            for (i = 1; i <= NF; ++i) {
                d = abs(a[FNR, i] - $i)
                if (d > tol * abs(a[FNR, i]) + 1e-12) {
                    print "row " FNR " col " i ": " a[FNR, i] " vs " $i; exit 1
                }
                ++n
            }
        }
        END {
            if (FNR != rows) { print "row count " rows " vs " FNR; exit 1 }
            if (n == 0) { print "no cells compared"; exit 1 }
        }
    ' "$1" "$2"
}

# 1. 'flat-missing' is byte-for-byte identical across every shape, since
#    nothing is weighted and both programs group the arithmetic the same way.
for shape in "" square square0 triangle; do
    suffix=${shape:-default}
    plink --bfile tmp_data --distance flat-missing $shape --out plink19_f_$suffix
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing $shape --out plink2_f_$suffix
    cmp plink19_f_$suffix.dist plink2_f_$suffix.dist
done

# ...including the double-precision binary matrix.
plink --bfile tmp_data --distance flat-missing bin --out plink19_fbin
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing bin --out plink2_fbin
cmp plink19_fbin.dist.bin plink2_fbin.dist.bin
plink --bfile tmp_data --distance flat-missing bin square0 --out plink19_fbin0
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing bin square0 --out plink2_fbin0
cmp plink19_fbin0.dist.bin plink2_fbin0.dist.bin

# 2. The IBS and 1-IBS reports, and the deprecated commands that wrap them.
plink --bfile tmp_data --distance flat-missing ibs 1-ibs allele-ct --out plink19_i
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing ibs 1-ibs allele-ct --out plink2_i
cmp plink19_i.dist plink2_i.dist
cmp plink19_i.mibs plink2_i.mibs
cmp plink19_i.mdist plink2_i.mdist
# The .mibs triangle carries its diagonal (self-IBS is 1); .dist and .mdist do
# not, so they have one row fewer.
test "$(wc -l < plink2_i.mibs)" -eq "$(wc -l < tmp_data.fam)"
test "$(wc -l < plink2_i.mdist)" -eq "$(($(wc -l < tmp_data.fam) - 1))"

plink --bfile tmp_data --distance-matrix --out plink19_dm
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance-matrix --out plink2_dm
cmp plink19_dm.mdist plink2_dm.mdist
plink --bfile tmp_data --ibs-matrix --out plink19_im
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --ibs-matrix --out plink2_im
cmp plink19_im.mibs plink2_im.mibs

# 3. The default, frequency-weighted missingness correction.  1.9's weights
#    are quantized, so this is a tolerance comparison.
plink --bfile tmp_data --distance --out plink19_w
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance --out plink2_w
close plink19_w.dist plink2_w.dist 1e-6
# The weighting has to actually do something, or the test above proves nothing.
if cmp -s plink19_f_default.dist plink19_w.dist; then
    echo "fixture has no missing calls; the weighted path is untested"
    exit 1
fi

# 4. Sample IDs, compression, threads and --parallel.
test "$(wc -l < plink2_w.dist.id)" -eq "$(($(wc -l < tmp_data.fam) + 1))"
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing zs --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.dist.zst > plink2_zs.dist
cmp plink2_f_default.dist plink2_zs.dist

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing --threads 1 --out plink2_t1
cmp plink2_f_default.dist plink2_t1.dist

for piece in 1 2 3; do
    # The loop's stdout is the concatenated matrix, so plink2's own console
    # output has to go somewhere else.
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing --parallel $piece 3 --out plink2_p > /dev/null
    cat plink2_p.dist.$piece
done > plink2_parallel.dist
cmp plink2_f_default.dist plink2_parallel.dist
# ...and each piece matches the one 1.9 writes.
plink --bfile tmp_data --distance flat-missing --parallel 2 3 --out plink19_p
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing --parallel 2 3 --out plink2_p2
cmp plink19_p.dist.2 plink2_p2.dist.2

# 5. Bad modifier combinations are rejected.
fails() {
    if "$@" > /dev/null 2>&1; then
        echo "expected failure: $*"
        exit 1
    fi
}
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance square triangle --out plink2_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance bin bin4 --out plink2_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance zs bin --out plink2_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance nonsense --out plink2_bad

# 6. Two deliberate divergences from PLINK 1.9, both on the monomorphic
#    variants excluded above and on 'bin4'.
#
# 6a. PLINK 1.9's weight loop leaves the frequency in place instead of
#     computing a weight when it is exactly 0 or 1, so a monomorphic variant
#     whose observed allele is A2 gets weight 1 -- more than five times the
#     largest weight any polymorphic variant can get (0.1875, at MAF 0.5).
#     plink2 gives it weight 0, which is what 1.9's own comment says it wants.
#     The two agree again once the variant is dropped.
cat > tmp_mono.ped << 'EOF'
F1 s1 0 0 1 1 A A  A A  A A  A G
F1 s2 0 0 1 1 A G  A A  G G  A G
F1 s3 0 0 1 1 G G  A A  A A  G G
F1 s4 0 0 1 1 0 0  A A  A G  A A
EOF
cat > tmp_mono.map << 'EOF'
1 v1 0 1
1 v2 0 2
1 v3 0 3
1 v4 0 4
EOF
plink --file tmp_mono --make-bed --out tmp_mono
plink --bfile tmp_mono --distance --out plink19_mono
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_mono --distance --out plink2_mono
if cmp -s plink19_mono.dist plink2_mono.dist; then
    echo "PLINK 1.9 no longer weights monomorphic variants at 1; drop this case"
    exit 1
fi
echo v2 > tmp_mono_drop.txt
plink --bfile tmp_mono --exclude tmp_mono_drop.txt --distance --out plink19_mono2
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_mono --exclude tmp_mono_drop.txt --distance --out plink2_mono2
cmp plink19_mono2.dist plink2_mono2.dist

# 6b. PLINK 1.9's single-precision square writer reuses the variable holding
#     the last lower-triangle value for the diagonal, so its diagonal is that
#     value instead of 0.  plink2 writes 0.  Everything off the diagonal
#     matches.
plink --bfile tmp_data --distance flat-missing bin4 --out plink19_b4
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --distance flat-missing bin4 --out plink2_b4
python3 - << 'EOF'
import struct
def rd(path):
    with open(path, 'rb') as f:
        data = f.read()
    return struct.unpack('<%df' % (len(data) // 4), data)
a = rd('plink19_b4.dist.bin')
b = rd('plink2_b4.dist.bin')
assert len(a) == len(b), (len(a), len(b))
n = int(round(len(a) ** 0.5))
assert n * n == len(a)
off_diag_diffs = 0
diag_zero = 0
for i, (x, y) in enumerate(zip(a, b)):
    if i // n == i % n:
        assert y == 0.0, 'plink2 diagonal at %d is %r' % (i, y)
        diag_zero += 1
    elif x != y:
        off_diag_diffs += 1
assert off_diag_diffs == 0, '%d off-diagonal differences' % off_diag_diffs
assert diag_zero == n
EOF
