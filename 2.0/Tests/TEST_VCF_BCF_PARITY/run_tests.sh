#!/bin/bash

# --vcf/--bcf parity under VCF import flags, with position sweeps.
#
# The VCF and BCF readers are separate hand-optimized implementations of the
# same import logic, so running identical data through both and diffing the
# results checks each against the other without hand-written expectations.
# BCF fixtures are written by Tests/utils/vcf_to_bcf.py instead of plink2's
# own --export bcf, which never emits fields like GQ/DP; a BcfScanGt() crash
# on --vcf-min-gq/--vcf-min-dp with a phased call in the first sample went
# unnoticed for years because no BCF test could reach that code.
#
# make_vcf.py puts each special genotype, alone, at the first and last sample
# and on both sides of every SIMD-vector / genovec-word boundary, for sample
# counts straddling those boundaries; see that script for the cases.
#
# For each (family, sample count, flags) combination, the two runs must
# either both fail with the same exit code, or both succeed and export
# identical VCF bodies.

set -eo pipefail

PLINK2="$1/plink2 $2 $3"
SAMPLE_CTS="1 2 3 7 8 9 15 16 17 31 32 33 63 64 65 127 128 129"

case_ct=0

# $1: family, $2: sample count, $3: --vcf/--bcf modifiers, $4: other flags,
# $5: expected outcome ("ok" or "fail")
check() {
    local fam=$1 n=$2 mods=$3 flags=$4 want=$5
    local rc_v=0 rc_b=0
    $PLINK2 --vcf tmp_${fam}_$n.vcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_v > /dev/null 2>&1 || rc_v=$?
    $PLINK2 --bcf tmp_${fam}_$n.bcf $mods $flags --export vcf vcf-dosage=HDS-force --out tmp_b > /dev/null 2>&1 || rc_b=$?
    case_ct=$((case_ct + 1))
    local desc="family=$fam samples=$n mods=[$mods] flags=[$flags]"
    if [ $rc_v -ne $rc_b ]; then
        echo "exit code mismatch (vcf $rc_v, bcf $rc_b): $desc"
        exit 1
    fi
    if [ $want = ok ] && [ $rc_v -ne 0 ]; then
        echo "unexpected failure (exit $rc_v): $desc"
        exit 1
    fi
    if [ $want = fail ] && [ $rc_v -eq 0 ]; then
        echo "unexpected success: $desc"
        exit 1
    fi
    if [ $rc_v -eq 0 ] && ! diff -q <(grep -v '^##' tmp_v.vcf) <(grep -v '^##' tmp_b.vcf) > /dev/null; then
        echo "output mismatch: $desc"
        diff <(grep -v '^##' tmp_v.vcf) <(grep -v '^##' tmp_b.vcf) | cut -c 1-300 | head -n 20
        exit 1
    fi
}

QUAL_FILTERS=("" "--vcf-min-gq 20" "--vcf-min-dp 10" "--vcf-max-dp 40" "--vcf-min-gq 20 --vcf-min-dp 10 --vcf-max-dp 40")

for n in $SAMPLE_CTS; do
    for fam in geno ds dsbad gp; do
        python3 make_vcf.py $fam $n tmp_${fam}_$n.vcf
        python3 ../utils/vcf_to_bcf.py tmp_${fam}_$n.vcf tmp_${fam}_$n.bcf
    done

    for hc in h m r; do
        for q in "${QUAL_FILTERS[@]}"; do
            check geno $n "" "--vcf-half-call $hc --polyploid-mode m $q" ok
        done
    done
    for extra in "--vcf-require-gt" "--vcf-require-gt --vcf-min-gq 20" "--vcf-ref-n-missing" "--var-min-qual 20" "--import-max-alleles 3" "--import-max-alleles 3 --vcf-min-dp 10"; do
        check geno $n "" "--vcf-half-call r --polyploid-mode m $extra" ok
    done
    # Default --vcf-half-call / --polyploid-mode modes are 'error'.
    check geno $n "" "--polyploid-mode m" fail
    check geno $n "" "--polyploid-mode m --vcf-min-gq 20" fail
    check geno $n "" "--vcf-half-call h" fail
    check geno $n "" "--vcf-half-call h --vcf-min-dp 10" fail

    for mods in "" "dosage=DS" "dosage=HDS"; do
        check ds $n "$mods" "" ok
    done
    check dsbad $n "" "" ok
    check dsbad $n "dosage=DS" "" ok
    check dsbad $n "dosage=HDS" "" fail
    for mods in "" "dosage=GP" "dosage=GP-force"; do
        check gp $n "$mods" "" ok
    done
    check gp $n "dosage=GP" "--import-dosage-certainty 0.6" ok

    rm -f tmp_*_$n.vcf tmp_*_$n.bcf
done

echo "$case_ct parity cases passed."
