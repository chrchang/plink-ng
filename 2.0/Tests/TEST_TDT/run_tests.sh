#!/bin/bash

# --tdt, checked against PLINK 1.9.
#
# The fixture is 150 trios with Mendelian transmission, plus 1% missing calls
# and 1% deliberate Mendel errors, which are the two cases a trio gets skipped
# for.

set -exo pipefail

awk -f make_fixture.awk
$1/plink2 $2 $3 --vcf tmp_trios.vcf --id-delim _ --make-bed --out tmp_data
cp tmp_ped.txt tmp_data.fam

# 1. Chi-square test.
plink --bfile tmp_data --tdt --out plink19
$1/plink2 $2 $3 --bfile tmp_data --tdt --out plink2
awk -v p19chisq=1 -v p19col=10 -v p2col=12 -f compare.awk plink19.tdt plink2.tdt

# 2. Exact binomial test, and its mid-p variant.  PLINK 1.9 keeps its CHISQ
#    column there and plink2 drops it, so the p-value moves.
plink --bfile tmp_data --tdt exact --out plink19_exact
$1/plink2 $2 $3 --bfile tmp_data --tdt exact --out plink2_exact
awk -v p19chisq=0 -v p19col=9 -v p2col=11 -f compare.awk plink19_exact.tdt plink2_exact.tdt

plink --bfile tmp_data --tdt exact-midp --out plink19_midp
$1/plink2 $2 $3 --bfile tmp_data --tdt exact-midp --out plink2_midp
awk -v p19chisq=0 -v p19col=9 -v p2col=11 -f compare.awk plink19_midp.tdt plink2_midp.tdt

# 3. With the affected children dropped there is no trio left, which has to be
#    an error rather than an empty report.
awk '$2 != "kid"' tmp_ped.txt | cut -d' ' -f1,2 > tmp_parents.txt
if $1/plink2 $2 $3 --bfile tmp_data --keep tmp_parents.txt --tdt --out plink2_noaff 2> tmp_noaff_err.txt; then
    echo "--tdt unexpectedly succeeded with no trios"
    exit 1
fi
grep -q "requires at least one trio" tmp_noaff_err.txt

# 4. cols= drops the columns it doesn't name.
$1/plink2 $2 $3 --bfile tmp_data --tdt cols=chrom,t,u,p --out plink2_cols
head -n 1 plink2_cols.tdt | grep -qx '#CHROM	ID	T	U	P'

# 5. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --tdt zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.tdt.zst > plink2_zs.tdt
diff -q plink2.tdt plink2_zs.tdt
