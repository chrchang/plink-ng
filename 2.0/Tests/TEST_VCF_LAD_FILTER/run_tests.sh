#!/bin/bash

# --vcf-min-lad / --vcf-max-lad: genotype filter on the sum of FORMAT/LAD.
# Checks FORMAT order (GQ before/after LAD, no GQ), missing and dropped LAD
# (count as 0), haploid ceil(min/2) at both short and >10-allele sites, LAD
# values needing int16 in BCF, BCF/VCF agreement, and that --vcf-min-dp is
# unaffected.  The phased record (v6) is VCF-only.

set -exo pipefail

python3 make_inputs.py

$1/plink2 $2 $3 --vcf lad.vcf --export vcf --out tmp_none
$1/plink2 $2 $3 --vcf lad.vcf --vcf-min-gq 20 --vcf-min-lad 10 --export vcf --out tmp_vcf
$1/plink2 $2 $3 --bcf lad.bcf --vcf-min-gq 20 --vcf-min-lad 10 --export vcf --out tmp_bcf
$1/plink2 $2 $3 --vcf lad.vcf --vcf-max-lad 250 --export vcf --out tmp_vcf_max
$1/plink2 $2 $3 --bcf lad.bcf --vcf-max-lad 250 --export vcf --out tmp_bcf_max
$1/plink2 $2 $3 --vcf dp.vcf --export vcf --out tmp_dp_none
$1/plink2 $2 $3 --vcf dp.vcf --vcf-min-dp 10 --export vcf --out tmp_dp

python3 - << 'PYEOF'
def gts(fname):
    rows = {}
    for line in open(fname):
        if not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            rows[f[2]] = f[9:]
    return rows

def check(filtered, unfiltered, expect_missing, skip=()):
    got, base = gts(filtered), gts(unfiltered)
    base = {k: v for k, v in base.items() if k not in skip}
    assert got.keys() == base.keys(), (filtered, got.keys())
    # compare calls only; separator display for homozygous/haploid calls varies
    norm = lambda g: "./." if g in (".", "./.", ".|.") else g.replace("|", "/")
    for vid, row in base.items():
        for s, g in enumerate(row):
            want = "./." if s in expect_missing.get(vid, ()) else norm(g)
            assert norm(got[vid][s]) == want, (filtered, vid, "abcdefgh"[s], got[vid][s], want)

# sample index: a=0 ... h=7
mins = {"v1": {1, 2, 3, 4, 7}, "v2": {1, 2, 3, 6}, "v3": {1, 3, 4, 5}, "v4": {1, 3, 6, 7}, "v5": {3, 5, 7}, "v6": {1, 2, 4, 6}}
check("tmp_vcf.vcf", "tmp_none.vcf", mins)
check("tmp_bcf.vcf", "tmp_none.vcf", mins, skip={"v6"})
maxes = {"v5": {4}}
check("tmp_vcf_max.vcf", "tmp_none.vcf", maxes)
check("tmp_bcf_max.vcf", "tmp_none.vcf", maxes, skip={"v6"})
# --vcf-min-dp still reads FORMAT/DP (LAD ignored; missing DP passes)
check("tmp_dp.vcf", "tmp_dp_none.vcf", {"d1": {0, 4, 7}})
PYEOF

# DP and LAD flags cannot be combined.
set +e
$1/plink2 $2 $3 --vcf lad.vcf --vcf-min-dp 5 --vcf-min-lad 5 --export vcf --out tmp_conflict > /dev/null 2>&1
rc=$?
set -e
[[ $rc -ne 0 && $rc -lt 128 ]] || exit 1
grep -q 'cannot be used with' tmp_conflict.log
