#!/bin/bash

# --lfile/--lgen, checked by round trip and against PLINK 1.9.
#
# The 200-sample fixture is deliberately one whose packed and aligned word
# counts differ: the per-variant genotype slices were spaced by the packed
# count, so every other slice was unaligned and GenovecInvertUnsafe() aborted.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

plink --simulate simulate.txt --simulate-missing 0.03 --simulate-ncases 100 --simulate-ncontrols 100 --out tmp_data > /dev/null
test "$(wc -l < tmp_data.fam)" -eq 200

# 1. Round trip: export to long format, read it back, and the genotypes have
#    to survive.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export lgen --out tmp_lg
$BUILD/plink2 $EXTRA1 $EXTRA2 --lfile tmp_lg --make-bed --out plink2_rt
cmp tmp_data.bed plink2_rt.bed
diff -q <(cut -f1,2,4 tmp_data.bim) <(cut -f1,2,4 plink2_rt.bim)

# 2. PLINK 1.9 reading the same .lgen has to agree.
plink --lfile tmp_lg --make-bed --out plink19_rt
cmp plink19_rt.bed plink2_rt.bed
diff -q plink19_rt.bim plink2_rt.bim

# 3. --lgen with the .map and .fam named separately.
$BUILD/plink2 $EXTRA1 $EXTRA2 --lgen tmp_lg.lgen --map tmp_lg.map --fam tmp_lg.fam --make-bed --out plink2_split
cmp plink2_rt.bed plink2_split.bed

# 4. --reference: calls absent from the .lgen become homozygous for the named
#    allele instead of missing.  Dropping every call of one variant and naming
#    its A2 as the reference has to reproduce a fully homozygous variant.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export lgen-ref --out tmp_lgref
test -s tmp_lgref.ref
$BUILD/plink2 $EXTRA1 $EXTRA2 --lfile tmp_lgref --reference tmp_lgref.ref --make-bed --out plink2_ref
plink --lfile tmp_lgref --reference tmp_lgref.ref --make-bed --out plink19_ref
cmp plink19_ref.bed plink2_ref.bed

# 5. The 'lgen-ref' export is smaller than the plain one, since the reference
#    calls are omitted.
test "$(wc -c < tmp_lgref.lgen)" -lt "$(wc -c < tmp_lg.lgen)"

# 6. A missing .map is an error rather than a silent empty result.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --lgen tmp_lg.lgen --fam tmp_lg.fam --make-bed --out plink2_nomap 2> tmp_nomap_err.txt; then
    echo "--lgen unexpectedly succeeded with no .map"
    exit 1
fi


# 7. Duplicate variant IDs are rejected rather than resolved to an arbitrary
#    variant.  (The duplicate-tolerant hash tables flag a duplicate in the high
#    bit of the stored index, which this lookup path would use as an index.)
awk 'NR == 2 { $2 = "dupid" } NR == 3 { $2 = "dupid" } { print }' tmp_lg.map > tmp_dup.map
cp tmp_lg.fam tmp_dup.fam
cp tmp_lg.lgen tmp_dup.lgen
if $BUILD/plink2 $EXTRA1 $EXTRA2 --lfile tmp_dup --make-bed --out plink2_dup 2> tmp_dup_err.txt; then
    echo "--lgen unexpectedly succeeded with duplicate variant IDs"
    exit 1
fi
grep -q 'duplicate variant IDs' tmp_dup_err.txt

# 8. When the genotype matrix does not fit in memory, the first block of
#    variants is filled in directly and the rest of the calls are spilled to a
#    temporary file that is reread once per block.  --memory cannot go low
#    enough to force that on a fixture this small, so
#    '--debug lgen-block-size=7' caps the block at 7 variants instead.  The
#    result has to match the in-memory import, and the temporary file has to
#    be gone afterwards.  Plain --debug must not cap it.
$BUILD/plink2 $EXTRA1 $EXTRA2 --debug lgen-block-size=7 --lfile tmp_lg --make-bed --out plink2_spill > plink2_spill.stdout
grep -q 'spilled to' plink2_spill.log
cmp plink2_rt.bed plink2_spill.bed
diff -q plink2_rt.bim plink2_spill.bim
test ! -e plink2_spill-temporary.lgen.tmp
$BUILD/plink2 $EXTRA1 $EXTRA2 --debug --lfile tmp_lg --make-bed --out plink2_nospill > /dev/null
if grep -q 'spilled to' plink2_nospill.log; then
    echo "plain --debug should not cap the --lgen block size"
    exit 1
fi
cmp plink2_rt.bed plink2_nospill.bed
$BUILD/plink2 $EXTRA1 $EXTRA2 --debug lgen-block-size=7 --lfile tmp_lgref --reference tmp_lgref.ref --make-bed --out plink2_spill_ref > /dev/null
cmp plink2_ref.bed plink2_spill_ref.bed

# 9. A later line for the same (sample, variant) pair overrides an earlier
#    one, and the spilled records have to keep that order: append a missing
#    call for every 97th line, then compare the two paths.
{
    cat tmp_lg.lgen
    awk 'NR % 97 == 0 { print $1, $2, $3, "0", "0" }' tmp_lg.lgen
} > tmp_ovr.lgen
cp tmp_lg.map tmp_ovr.map
cp tmp_lg.fam tmp_ovr.fam
$BUILD/plink2 $EXTRA1 $EXTRA2 --lfile tmp_ovr --make-bed --out plink2_ovr
$BUILD/plink2 $EXTRA1 $EXTRA2 --debug lgen-block-size=7 --lfile tmp_ovr --make-bed --out plink2_ovr_spill > /dev/null
cmp plink2_ovr.bed plink2_ovr_spill.bed
if cmp -s plink2_rt.bed plink2_ovr.bed; then
    echo "the overriding missing calls had no effect"
    exit 1
fi
