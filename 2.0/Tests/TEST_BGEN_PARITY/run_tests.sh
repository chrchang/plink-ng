#!/bin/bash

# --bgen/--vcf parity with position sweeps.
#
# make_data.py quantizes each probability to the BGEN scale and writes both a
# BGEN file (via Tests/utils/bgen_writer.py, not plink2's own --export bgen,
# which only emits a few fixed layouts) and a VCF containing the decimal
# value of exactly the encoded rationals.  The .bgen and VCF importers are
# separate implementations (integer magic-number division vs. floating
# point), so the resulting .pgen files must be byte-identical; .pvar
# positions/IDs/alleles and sample IDs must match too.
#
# Unphased data is compared against VCF GP (dosage=GP), phased against VCF
# HDS (dosage=HDS).  Coverage:
# * layout 2 at bit depths 1..28 (all of them at 33 samples; 1, 3, 8, 16, 28
#   elsewhere), plus layout 1;
# * at 8 and 16 bits (the SIMD fast paths), no compression / zlib / zstd,
#   no sample ID block (--sample), and 'ref-last' with the allele order
#   reversed and probabilities mirrored;
# * fixed- and variable-ploidy records (haploid, ploidy 0, missing samples);
# * --hard-call-threshold, --dosage-erase-threshold and
#   --import-dosage-certainty (including a certainty exactly equal to an
#   encoded probability);
# * each special genotype alone at the first/last sample and around
#   SIMD-vector / genovec-word boundaries, for sample counts straddling them.
#
# Known differences that the fixtures avoid:
# * Layout-1 import rounds exact dosage ties to even, while VCF GP import
#   rounds them up (deliberate).
# * --dosage-erase-threshold on phased dosages with equal haplotypes near 0.5
#   (e.g. 0.498|0.498): .bgen import applies the unphased rule (dosage 0.996
#   is within the threshold of 1, so it's erased to a het hardcall), while
#   VCF/BCF HDS conversion applies the phased rule (each haplotype must be
#   within threshold/2 of an integer, so the dosage is kept).  The phased
#   family is therefore not run with --dosage-erase-threshold here.

set -eo pipefail

PLINK2="$1/plink2 $2 $3"
SAMPLE_CTS="1 2 3 4 5 7 8 9 15 16 17 31 32 33 63 64 65"
FLAGS=("" "--hard-call-threshold 0.05" "--dosage-erase-threshold 0.05" "--import-dosage-certainty 0.6" "--import-dosage-certainty 0.9")
EXTRA_FLAGS=("" "--import-dosage-certainty 0.6")

case_ct=0
mismatch_ct=0

iids() {
    awk 'NR == 1 {for (i = 1; i <= NF; i++) if (($i == "IID") || ($i == "#IID")) c = i; next} {print $c}' $1
}

# $1: VCF/BGEN prefix, $2: dosage field, $3: encoding, $4: flags
check() {
    local pre=$1 dfield=$2 enc=$3 flags=$4
    local bgen_args="--bgen $pre.$enc.bgen ref-first"
    case $enc in
        reflast) bgen_args="--bgen $pre.$enc.bgen ref-last" ;;
        noids) bgen_args="$bgen_args --sample $pre.noids.sample" ;;
    esac
    if ! $PLINK2 $bgen_args $flags --make-pgen --out tmp_b > /dev/null 2>&1; then
        echo "--bgen failed: $pre $enc [$flags]"
        exit 1
    fi
    case_ct=$((case_ct + 1))
    # (.pvar and .psam contents don't depend on the flags, so they're only
    # compared once per encoding.)
    if ! cmp -s tmp_v.pgen tmp_b.pgen || \
       { [ -z "$flags" ] && \
         { ! diff -q <(grep -v '^#' tmp_v.pvar | cut -f 1-5) <(grep -v '^#' tmp_b.pvar | cut -f 1-5) > /dev/null || \
           ! diff -q <(iids tmp_v.psam) <(iids tmp_b.psam) > /dev/null; }; }; then
        echo "mismatch: $pre $enc [$flags]"
        $PLINK2 --pfile tmp_v --export vcf vcf-dosage=HDS-force --out tmp_ve > /dev/null 2>&1
        $PLINK2 --pfile tmp_b --export vcf vcf-dosage=HDS-force --out tmp_be > /dev/null 2>&1
        diff <(grep -v '^##' tmp_ve.vcf) <(grep -v '^##' tmp_be.vcf) | cut -c 1-300 | head -n 20 || true
        if [ -z "$PARITY_KEEP_GOING" ]; then
            exit 1
        fi
        mismatch_ct=$((mismatch_ct + 1))
    fi
}

# $1: family, $2: sample count, $3: bits or L1, $4: encodings
run_set() {
    local fam=$1 n=$2 bits=$3 encs=$4
    local pre=tmp_${fam}_${n}_$bits
    local dfield=GP
    [ $fam = hds ] && dfield=HDS
    [ -n "$PREGENERATED" ] || python3 make_data.py $fam $n $bits $pre $encs
    for flags in "${FLAGS[@]}"; do
        if [ $fam = hds ] && [[ $flags == --dosage-erase-threshold* ]]; then
            continue
        fi
        if ! $PLINK2 --vcf $pre.vcf dosage=$dfield $flags --make-pgen --out tmp_v > /dev/null 2>&1; then
            echo "--vcf failed: $pre [$flags]"
            exit 1
        fi
        for enc in ${encs//,/ }; do
            [ -f $pre.$enc.bgen ] || continue
            if [[ $enc =~ ^(zlib|l1zlib|l1none)$ ]]; then
                check $pre $dfield $enc "$flags"
            else
                for ef in "${EXTRA_FLAGS[@]}"; do
                    if [ "$ef" = "$flags" ]; then
                        check $pre $dfield $enc "$flags"
                    fi
                done
            fi
        done
    done
    rm -f $pre.*
}

# The weekly workflow (auto_checkcommits_tests_linuxmac.yml) sets
# PLINK2_FULL_TESTS=1 for the full sweep below.  Otherwise, only 33 samples
# (straddling the 32-sample genovec word) are run, with one bit depth per
# importer regime: the 1-bit case, the 8- and 16-bit fast paths, layout 1,
# and 25 bits for both the generic unpacker and the >24-bit rearrangement
# of the magic-number division.  Each family gets one of the two fast paths.
# (These fixtures are generated in parallel first; "wait $pid" propagates
# failures.)
if [ -z "$PLINK2_FULL_TESTS" ]; then
    QUICK_SETS=("gp 33 1 zlib" "gp 33 16 zlib,none,zstd,noids,reflast" "hds 33 8 zlib,none,zstd,noids,reflast" "hds 33 25 zlib" "gp 33 L1 l1zlib,l1none")
    pids=""
    for set_args in "${QUICK_SETS[@]}"; do
        read -r fam n bits encs <<< "$set_args"
        python3 make_data.py $fam $n $bits tmp_${fam}_${n}_$bits $encs &
        pids="$pids $!"
    done
    for pid in $pids; do
        wait $pid
    done
    PREGENERATED=1
    for set_args in "${QUICK_SETS[@]}"; do
        run_set $set_args
    done
    SAMPLE_CTS=""
fi

for n in $SAMPLE_CTS; do
    if [ $n = 33 ]; then
        bit_list=$(seq 1 28)
    else
        bit_list="1 3 8 16 28"
    fi
    for fam in gp hds; do
        for bits in $bit_list; do
            if [ $bits = 8 ] || [ $bits = 16 ]; then
                run_set $fam $n $bits zlib,none,zstd,noids,reflast
            else
                run_set $fam $n $bits zlib
            fi
        done
    done
    run_set gp $n L1 l1zlib,l1none
done

if [ $mismatch_ct -ne 0 ]; then
    echo "$mismatch_ct of $case_ct parity cases failed."
    exit 1
fi
echo "$case_ct parity cases passed."
