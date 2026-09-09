#!/bin/bash

# "--make-king-table cols=+rt,+pkin": the two pedigree-derived columns.
#
# These describe the pedigree rather than the genotypes, so they are checked
# against a hand-built pedigree whose expected kinship coefficients are known
# in closed form.  The genotypes are simulated and mutually independent, which
# is deliberate: it keeps the observed KINSHIP column irrelevant to what is
# being tested here.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

# 1. Build the fileset.  Only the .fam matters; the genotypes just have to
#    exist.
printf '1000 null 0.05 0.95 1 1\n' > tmp_sim.txt
plink --simulate tmp_sim.txt --make-bed --out tmp_base --seed 42 > /dev/null

# The pedigree, all in one family:
#   s1 x s2  -> s3, s4        (s3/s4 full siblings)
#   s1 x s6  -> s5            (s5 half sibling of s3 and s4)
#   s3 x s7  -> s8            (s8 grandchild of s1/s2, nephew of s4,
#                              half-nephew of s5)
cat > tmp_ped.fam << 'PEDEOF'
F1 s1 0 0 1 1
F1 s2 0 0 2 1
F1 s3 s1 s2 1 1
F1 s4 s1 s2 2 1
F1 s6 0 0 2 1
F1 s5 s1 s6 1 1
F1 s7 0 0 2 1
F1 s8 s3 s7 1 1
PEDEOF

head -8 tmp_base.fam | awk '{print $1, $2}' > tmp_keep.txt
plink --bfile tmp_base --keep tmp_keep.txt --make-bed --out tmp_sub > /dev/null
# (via a temporary: redirecting straight onto tmp_sub.fam would truncate it
# before awk reads it)
paste -d' ' <(cut -d' ' -f1-4 tmp_ped.fam) <(awk '{print $5, $6}' tmp_sub.fam) > tmp_sub_fam.new
mv tmp_sub_fam.new tmp_sub.fam

# 2. Every pair is reported, including the unrelated ones.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_sub --make-king-table cols=+rt,+pkin --king-table-filter -9 --out tmp_kt

if [[ $(grep -cv '^#' tmp_kt.kin0) != 28 ]]; then
    echo "--make-king-table did not report all 28 pairs."
    exit 1
fi

# 3. Each labelled relationship, with the kinship coefficient the pedigree
#    implies.  Column order is FID1 IID1 FID2 IID2 ... RT PEDIGREE_KINSHIP.
check_pair() {
    local iid1=$1
    local iid2=$2
    local expected_rt=$3
    local expected_pkin=$4
    local observed
    observed=$(awk -v a="$iid1" -v b="$iid2" '($2 == a && $4 == b) || ($2 == b && $4 == a) {print $(NF-1), $NF}' tmp_kt.kin0)
    if [[ "$observed" != "$expected_rt $expected_pkin" ]]; then
        echo "$iid1/$iid2: expected '$expected_rt $expected_pkin', got '$observed'."
        exit 1
    fi
}

check_pair s1 s3 PO 0.25     # parent/offspring
check_pair s2 s4 PO 0.25
check_pair s3 s8 PO 0.25
check_pair s3 s4 FS 0.25     # full siblings
check_pair s3 s5 HS 0.125    # half siblings, shared father
check_pair s4 s5 HS 0.125
check_pair s1 s8 GG 0.125    # grandparent/grandchild
check_pair s2 s8 GG 0.125
check_pair s4 s8 AV 0.125    # aunt/nephew
check_pair s5 s8 REL 0.0625  # half-aunt/nephew: named by no simpler label
check_pair s1 s2 UN 0        # spouses, unrelated to each other
check_pair s2 s5 UN 0        # s5 is the father's child by another mother
check_pair s6 s3 UN 0

# 4. Neither column is emitted unless asked for, and 'rt' and 'pkin' are
#    independent of one another.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_sub --make-king-table --out tmp_kt_default
if grep -q 'RT' tmp_kt_default.kin0; then
    echo "RT column present without cols=+rt."
    exit 1
fi

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_sub --make-king-table cols=+pkin --out tmp_kt_pkin
if grep -q 'RT' tmp_kt_pkin.kin0; then
    echo "RT column present with cols=+pkin alone."
    exit 1
fi
if ! grep -q 'PEDIGREE_KINSHIP' tmp_kt_pkin.kin0; then
    echo "PEDIGREE_KINSHIP column missing with cols=+pkin."
    exit 1
fi

# 5. A fileset with no parental IDs at all: every pair is unrelated, and the
#    pedigree machinery must not misreport that as anything else.
awk '{print $1, $2, 0, 0, $5, $6}' tmp_sub.fam > tmp_nofam.fam
cp tmp_sub.bed tmp_nofam.bed
cp tmp_sub.bim tmp_nofam.bim
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_nofam --make-king-table cols=+rt,+pkin --king-table-filter -9 --out tmp_kt_nofam
if [[ $(awk '!/^#/ && ($(NF-1) != "UN" || $NF != 0)' tmp_kt_nofam.kin0 | wc -l | tr -d '[:space:]') != 0 ]]; then
    echo "Pedigree-free fileset reported a relationship."
    exit 1
fi

echo "--make-king-table rt/pkin tests passed."
