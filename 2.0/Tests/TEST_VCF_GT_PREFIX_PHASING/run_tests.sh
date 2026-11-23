#!/bin/bash

# VCFv4.4+ GT prefix phasing ("|0/1", "/0|1", "|1").
#
# plink2 skips the optional leading phasing indicator; for a diploid call,
# phase comes from the indicator between the two alleles.  This matches the
# BCF importer, which only looks at the second allele's phase bit.
#
# For each family and sample count, make_vcf.py writes a prefixed VCF and a
# canonical copy with the leading indicators removed, and
# Tests/utils/vcf_to_bcf.py converts the prefixed VCF to BCF the way htslib
# 1.22 does: the first allele's phase bit is set, and every GT vector in a
# record gets an extra END_OF_VECTOR slot when any of its samples has a
# prefix (so all-diploid records have width-3 GT vectors, which used to be
# rejected as polyploid).  Prefixed VCF and BCF must both import exactly like
# the canonical VCF.  Sample counts straddle the 32-sample genovec-word
# boundary, so the GT-only family also runs into the SWAR fast path.

set -eo pipefail

PLINK2="$1/plink2 $2 $3"
SAMPLE_CTS="1 2 31 32 33 34 65 66"

case_ct=0

# $1: family, $2: sample count, $3: --vcf/--bcf modifiers, $4: other flags,
# $5: expected outcome ("ok" or "fail")
check() {
    local fam=$1 n=$2 mods=$3 flags=$4 want=$5
    local rc_c=0 rc_v=0 rc_b=0
    $PLINK2 --vcf tmp_${fam}_${n}_canon.vcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_c > /dev/null 2>&1 || rc_c=$?
    $PLINK2 --vcf tmp_${fam}_$n.vcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_v > /dev/null 2>&1 || rc_v=$?
    $PLINK2 --bcf tmp_${fam}_$n.bcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_b > /dev/null 2>&1 || rc_b=$?
    case_ct=$((case_ct + 1))
    local desc="family=$fam samples=$n mods=[$mods] flags=[$flags]"
    if [ $rc_c -ne $rc_v ] || [ $rc_c -ne $rc_b ]; then
        echo "exit code mismatch (canonical $rc_c, prefixed vcf $rc_v, bcf $rc_b): $desc"
        exit 1
    fi
    if [ $want = ok ] && [ $rc_c -ne 0 ]; then
        echo "unexpected failure (exit $rc_c): $desc"
        exit 1
    fi
    if [ $want = fail ] && [ $rc_c -eq 0 ]; then
        echo "unexpected success: $desc"
        exit 1
    fi
    if [ $rc_c -eq 0 ]; then
        for out in tmp_v tmp_b; do
            if ! diff -q <(grep -v '^##' tmp_c.vcf) <(grep -v '^##' $out.vcf) > /dev/null; then
                echo "output mismatch ($out vs canonical): $desc"
                diff <(grep -v '^##' tmp_c.vcf) <(grep -v '^##' $out.vcf) | cut -c 1-300 | head -n 20
                exit 1
            fi
        done
    fi
}

# (--import-dosage-certainty 0.72 below: generated DS/HDS/GP values are multiples
# of 0.05, so none sit exactly on a threshold, where BCF's float32 and VCF's
# double parsing could legitimately round to opposite sides.)
for n in $SAMPLE_CTS; do
    for fam in geno gt ds hds gp poly; do
        python3 make_vcf.py $fam $n tmp_${fam}_$n.vcf tmp_${fam}_${n}_canon.vcf
        python3 ../utils/vcf_to_bcf.py tmp_${fam}_$n.vcf tmp_${fam}_$n.bcf
    done
    for flags in "--vcf-half-call h" "--vcf-half-call r" "--vcf-half-call m" "--vcf-half-call m --polyploid-mode missing"; do
        check geno $n "" "$flags" ok
        check gt $n "" "$flags" ok
    done
    check geno $n "" "--vcf-half-call r --vcf-min-gq 20 --vcf-min-dp 10 --vcf-max-dp 40" ok
    check geno $n "" "" fail  # half-calls error out by default
    for fam in ds hds; do
        check $fam $n "dosage=DS" "--vcf-half-call h" ok
        check $fam $n "dosage=HDS" "--vcf-half-call h" ok
        check $fam $n "dosage=HDS" "--vcf-half-call h --import-dosage-certainty 0.72" ok
    done
    check gp $n "dosage=GP-force" "--vcf-half-call h" ok
    check gp $n "dosage=GP-force" "--vcf-half-call h --import-dosage-certainty 0.72" ok
    check poly $n "" "--vcf-half-call m" fail
    check poly $n "" "--vcf-half-call m --polyploid-mode missing" ok
done

# Spot-check the rule itself, so the canonical comparison isn't the only
# thing pinning it down.  Record v9 is all 0/0 except a final "|0|1"; v10
# ends in "|1/0"; v7 (haploid only) starts with "|.", "1", "0", "|1", "/0".
$PLINK2 --vcf tmp_gt_33.vcf --vcf-half-call m --export vcf --out tmp_spot > /dev/null
python3 - << 'PYEOF'
rows = {}
for line in open("tmp_spot.vcf"):
    if not line.startswith("#"):
        f = line.rstrip("\n").split("\t")
        rows[f[2]] = f[9:]
# hets: "/0|1", "|0|1" and "/1|0" are phased; "|1/0" and "|0/1" aren't
assert rows["v1"][1:6] == ["0|1", "0|1", "0/1", "1|0", "0/1"], rows["v1"][1:6]
assert rows["v9"][-1] == "0|1", rows["v9"][-1]
assert rows["v10"][-1] == "0/1", rows["v10"][-1]
# haploid calls: compare alleles only (export separator display varies)
assert [g[0] for g in rows["v7"][:5]] == [".", "1", "0", "1", "0"], rows["v7"][:5]
PYEOF

# Only one leading indicator is allowed.
printf '##fileformat=VCFv4.4\n##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n1\t100\tv1\tA\tC\t.\t.\t.\tGT\t0/1\t||0/1\n' > tmp_double.vcf
rc=0
$PLINK2 --vcf tmp_double.vcf --export vcf --out tmp_double > /dev/null 2>&1 || rc=$?
[ $rc -ne 0 ] && [ $rc -lt 128 ]
grep -q 'invalid GT field' tmp_double.log

echo "$case_ct parity cases passed."
