#!/bin/bash

# VCFv4.4+ GT prefix phasing ("|0/1", "/0|1", "|1").
#
# Each phasing indicator applies to the allele after it, and plink2 treats a
# diploid call as phased only if both alleles are phased: "0|1" and "|0|1"
# are phased; "/0|1", "|0/1", and "/0/1" aren't.
#
# For each family and sample count, make_vcf.py writes a prefixed VCF, a
# canonical copy with that rule applied (and no prefixes), and a
# "BCF-expected" copy.  Tests/utils/vcf_to_bcf.py converts the prefixed VCF to
# BCF the way HTSlib 1.22+ does: the first allele's phase bit is set, and
# every GT vector in a record gets an extra END_OF_VECTOR slot when any of its
# samples has a prefix (so all-diploid records have width-3 GT vectors, which
# used to be rejected as polyploid).  Since BCF writers before HTSlib 1.22
# never set the first-allele bit, plink2 only uses it in records where some
# first value has it set; the BCF-expected file encodes that.  Checks:
#   * prefixed VCF imports like the canonical VCF;
#   * the BCF imports like the BCF-expected VCF;
#   * the canonical VCF, converted to BCF with legacy (pre-1.22) phase
#     encoding, imports like the canonical VCF, i.e. legacy phased calls stay
#     phased.
# Sample counts straddle the 32-sample genovec-word boundary, so the GT-only
# family also runs into the SWAR fast path.

set -eo pipefail

PLINK2="$1/plink2 $2 $3"
# The weekly workflow (auto_checkcommits_tests_linuxmac.yml) sets
# PLINK2_FULL_TESTS=1 to sweep every sample count and --vcf-half-call mode.
# Otherwise only 33 samples (just past a 32-sample genovec word; the spot
# checks below also use it) and two half-call modes are run.
if [ -n "$PLINK2_FULL_TESTS" ]; then
    SAMPLE_CTS="1 2 31 32 33 34 65 66"
    HALF_CALL_FLAGS=("--vcf-half-call h" "--vcf-half-call r" "--vcf-half-call m" "--vcf-half-call m --polyploid-mode missing")
else
    SAMPLE_CTS="33"
    HALF_CALL_FLAGS=("--vcf-half-call h" "--vcf-half-call m --polyploid-mode missing")
fi

case_ct=0

# $1: family, $2: sample count, $3: --vcf/--bcf modifiers, $4: other flags,
# $5: expected outcome ("ok" or "fail")
check() {
    local fam=$1 n=$2 mods=$3 flags=$4 want=$5
    local rc_c=0 rc_v=0 rc_e=0 rc_b=0 rc_l=0
    $PLINK2 --vcf tmp_${fam}_${n}_canon.vcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_c > /dev/null 2>&1 || rc_c=$?
    $PLINK2 --vcf tmp_${fam}_$n.vcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_v > /dev/null 2>&1 || rc_v=$?
    $PLINK2 --vcf tmp_${fam}_${n}_bcfexp.vcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_e > /dev/null 2>&1 || rc_e=$?
    $PLINK2 --bcf tmp_${fam}_$n.bcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_b > /dev/null 2>&1 || rc_b=$?
    $PLINK2 --bcf tmp_${fam}_${n}_legacy.bcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_l > /dev/null 2>&1 || rc_l=$?
    case_ct=$((case_ct + 1))
    local desc="family=$fam samples=$n mods=[$mods] flags=[$flags]"
    if [ $rc_c -ne $rc_v ] || [ $rc_c -ne $rc_e ] || [ $rc_c -ne $rc_b ] || [ $rc_c -ne $rc_l ]; then
        echo "exit code mismatch (canonical $rc_c, prefixed vcf $rc_v, bcf-expected $rc_e, bcf $rc_b, legacy bcf $rc_l): $desc"
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
        for pair in tmp_c:tmp_v tmp_e:tmp_b tmp_c:tmp_l; do
            local want_f=${pair%%:*} got_f=${pair##*:}
            if ! diff -q <(grep -v '^##' $want_f.vcf) <(grep -v '^##' $got_f.vcf) > /dev/null; then
                echo "output mismatch ($got_f vs $want_f): $desc"
                diff <(grep -v '^##' $want_f.vcf) <(grep -v '^##' $got_f.vcf) | cut -c 1-300 | head -n 20
                exit 1
            fi
        done
    fi
}

# (--import-dosage-certainty 0.72 below: generated DS/HDS/GP values are multiples
# of 0.05, so none sit exactly on a threshold, where BCF's float32 and VCF's
# double parsing could legitimately round to opposite sides.)
for n in $SAMPLE_CTS; do
    # (fixtures are generated in parallel; "wait $pid" propagates failures)
    pids=""
    for fam in geno gt ds hds gp poly; do
        (python3 make_vcf.py $fam $n tmp_${fam}_$n.vcf tmp_${fam}_${n}_canon.vcf tmp_${fam}_${n}_bcfexp.vcf
         python3 ../utils/vcf_to_bcf.py tmp_${fam}_$n.vcf tmp_${fam}_$n.bcf
         python3 ../utils/vcf_to_bcf.py --legacy-phase tmp_${fam}_${n}_canon.vcf tmp_${fam}_${n}_legacy.bcf) &
        pids="$pids $!"
    done
    for pid in $pids; do
        wait $pid
    done
    for flags in "${HALF_CALL_FLAGS[@]}"; do
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
# thing pinning it down.  Records v8/v9/v10 are all 0/0 except a final
# "/0|1"/"|0|1"/"|1/0"; v7 (haploid only) starts with "|.", "1", "0", "|1",
# "/0".
$PLINK2 --vcf tmp_gt_33.vcf --vcf-half-call m --export vcf --out tmp_spot > /dev/null
python3 - << 'PYEOF'
rows = {}
for line in open("tmp_spot.vcf"):
    if not line.startswith("#"):
        f = line.rstrip("\n").split("\t")
        rows[f[2]] = f[9:]
# hets: only "|0|1" is phased; "/0|1", "|1/0", "/1|0", and "|0/1" aren't
assert rows["v1"][1:6] == ["0/1", "0|1", "0/1", "0/1", "0/1"], rows["v1"][1:6]
assert rows["v8"][-1] == "0/1", rows["v8"][-1]
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
