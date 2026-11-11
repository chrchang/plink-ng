#!/bin/bash

# Regression tests for the PLINK 1.9 --dfam and --mds-plot eigendecomp fixes.
# Each check failed before its fix.  (Both were verified against PLINK 1.07;
# plink2 has neither command, so these checks compare PLINK 1.9 against
# itself on equivalent inputs.)
#
# Usage: ./run_tests.sh
# Requires plink (1.9) to be in the system PATH.

set -uo pipefail

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
# The columns of a .dfam report other than OBS and EXP.
dfam_stats() {
    awk '{ print $1, $2, $3, $4, $7, $8 }' "$1"
}

# --dfam: in a family with both parents genotyped, the variance denominator
# counted unaffected children too; only affected children belong there (as
# in PLINK 1.07's dfam.cpp).  So adding an unaffected child to each such
# family must not change CHISQ or P.  (OBS and EXP do both shift, in PLINK
# 1.07 as well.)  There are no missing calls: when a parent's genotype is
# missing, the family falls back to the sibship calculation, where
# unaffected children legitimately count.
awk -v extra=0 -f make_dfam.awk
mv tmp_df.ped tmp_df0.ped
mv tmp_df.map tmp_df0.map
awk -v extra=1 -f make_dfam.awk
plink --file tmp_df0 --dfam --out tmp_df0 > /dev/null
plink --file tmp_df --dfam --out tmp_df1 > /dev/null
test -s tmp_df0.dfam && test "$(dfam_stats tmp_df0.dfam)" = "$(dfam_stats tmp_df1.dfam)"
check "--dfam ignores unaffected children of genotyped parents" $?

# --dfam with --remove (or other sample filters) and --within clusters of
# unrelated samples: copy_bitarr_subset_excl() was given the filtered sample
# count as the bit array length, so samples near the end of the file lost
# their case/control status.  Removing samples on the fly must give the same
# result as removing them beforehand.
plink --file tmp_df --make-bed --out tmp_df > /dev/null
plink --bfile tmp_df --within tmp_df_within.txt --remove tmp_df_remove.txt --dfam --out tmp_df_rm > /dev/null
plink --bfile tmp_df --remove tmp_df_remove.txt --make-bed --out tmp_df_pre > /dev/null
plink --bfile tmp_df_pre --within tmp_df_within.txt --dfam --out tmp_df_pre > /dev/null
test -s tmp_df_pre.dfam && cmp -s tmp_df_rm.dfam tmp_df_pre.dfam
check "--dfam with --remove matches a pre-filtered fileset" $?

# --mds-plot eigendecomp: the double-centering loop started at 1, leaving
# B[0][0] uncentered.  It must give the same coordinates as the default (SVD)
# path, up to the sign of each dimension, both per sample and by cluster.
awk -f make_mds.awk
plink --file tmp_mds --make-bed --out tmp_mds > /dev/null
mds_ok=0
for mode in "" "by-cluster"; do
    plink --bfile tmp_mds --cluster --K 12 --mds-plot 4 $mode --out tmp_mds_svd > /dev/null
    plink --bfile tmp_mds --cluster --K 12 --mds-plot 4 $mode eigendecomp --out tmp_mds_eig > /dev/null
    awk -f compare_mds.awk tmp_mds_svd.mds tmp_mds_eig.mds || mds_ok=1
done
check "--mds-plot eigendecomp matches the default MDS" $mds_ok

exit $failed
