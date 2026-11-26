#!/bin/bash

# --vcf-min-ad-sum / --vcf-max-ad-sum: genotype filter on the sum of
# FORMAT/LAD, or FORMAT/AD when LAD is absent.
# Checks FORMAT order (GQ before/after LAD, no GQ), missing and dropped LAD
# (count as 0), haploid ceil(min/2) at both short and >10-allele sites, LAD
# values needing int16 in BCF, BCF/VCF agreement, and that --vcf-min-dp is
# unaffected.  The phased record (v6) is VCF-only.
# AD: used when a line has no LAD; LAD takes precedence whichever comes first
# (for both min and max); AD is used without a LAD header line, and an
# undeclared LAD is then ignored (VCF only); with neither header line, the
# flags are ignored with a warning.

set -exo pipefail

python3 make_inputs.py

$1/plink2 $2 $3 --vcf lad.vcf --export vcf --out tmp_none
$1/plink2 $2 $3 --vcf lad.vcf --vcf-min-gq 20 --vcf-min-ad-sum 10 --export vcf --out tmp_vcf
$1/plink2 $2 $3 --bcf lad.bcf --vcf-min-gq 20 --vcf-min-ad-sum 10 --export vcf --out tmp_bcf
$1/plink2 $2 $3 --vcf lad.vcf --vcf-max-ad-sum 250 --export vcf --out tmp_vcf_max
$1/plink2 $2 $3 --bcf lad.bcf --vcf-max-ad-sum 250 --export vcf --out tmp_bcf_max
$1/plink2 $2 $3 --vcf ad.vcf --export vcf --out tmp_ad_none
$1/plink2 $2 $3 --vcf ad.vcf --vcf-min-gq 20 --vcf-min-ad-sum 10 --export vcf --out tmp_ad_vcf
$1/plink2 $2 $3 --bcf ad.bcf --vcf-min-gq 20 --vcf-min-ad-sum 10 --export vcf --out tmp_ad_bcf
$1/plink2 $2 $3 --vcf ad.vcf --vcf-max-ad-sum 250 --export vcf --out tmp_ad_vcf_max
$1/plink2 $2 $3 --bcf ad.bcf --vcf-max-ad-sum 250 --export vcf --out tmp_ad_bcf_max
$1/plink2 $2 $3 --vcf adonly.vcf --export vcf --out tmp_adonly_none
$1/plink2 $2 $3 --vcf adonly.vcf --vcf-min-ad-sum 10 --export vcf --out tmp_adonly
$1/plink2 $2 $3 --vcf nodepth.vcf --vcf-min-ad-sum 10 --export vcf --out tmp_nodepth
grep -q 'No FORMAT/LAD or FORMAT/AD key found' tmp_nodepth.log
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
ad_mins = {"a1": {1, 3, 4, 5}, "a2": {0, 2, 3, 7}, "a3": {1, 3, 4, 6}, "a4": {1, 2, 3, 4}}
check("tmp_ad_vcf.vcf", "tmp_ad_none.vcf", ad_mins)
check("tmp_ad_bcf.vcf", "tmp_ad_none.vcf", ad_mins)
ad_maxes = {"a3": {5}, "a4": {7}}
check("tmp_ad_vcf_max.vcf", "tmp_ad_none.vcf", ad_maxes)
check("tmp_ad_bcf_max.vcf", "tmp_ad_none.vcf", ad_maxes)
check("tmp_adonly.vcf", "tmp_adonly_none.vcf", {"z1": {0, 2, 5}})
# --vcf-min-dp still reads FORMAT/DP (LAD ignored; missing DP passes)
check("tmp_dp.vcf", "tmp_dp_none.vcf", {"d1": {0, 4, 7}})
PYEOF

# DP and AD-sum flags cannot be combined.
set +e
$1/plink2 $2 $3 --vcf lad.vcf --vcf-min-dp 5 --vcf-min-ad-sum 5 --export vcf --out tmp_conflict > /dev/null 2>&1
rc=$?
set -e
[[ $rc -ne 0 && $rc -lt 128 ]] || exit 1
grep -q 'cannot be used with' tmp_conflict.log
