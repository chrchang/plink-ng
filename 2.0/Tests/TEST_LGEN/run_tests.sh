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
