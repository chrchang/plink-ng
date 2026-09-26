#!/bin/bash

# A .pgen record whose type byte says "1-bit" but whose body is not a valid
# 1-bit record must be rejected, not counted.
#
# The first body byte of a 1-bit record packs two genotype codes, and
# CountparseOnebitSubset() uses both of them to index a four-element
# genocounts[] array.  Before the accompanying fix, a corrupted byte indexed
# that array with a value of up to 66, and --freq still exited 0 with a report
# full of wrong numbers.

set -exo pipefail

$1/plink2 $2 $3 --dummy 400 600 --seed 1 --make-pgen --out tmp_data

# Sanity check: the unmodified fileset is fine.
$1/plink2 $2 $3 --pfile tmp_data --freq --out plink2_ok
test "$(grep -cv '^#' plink2_ok.afreq)" -eq 600

# Byte 92 is inside the record-type array, and setting its low nibble to 1
# marks a record as 1-bit whose body then decodes to an out-of-range genotype
# code.  (If a .pgen layout change moves that array, this offset is what needs
# updating; the scan that found it just tried every header byte.)  The high
# nibble keeps its original value, 1: 0x91 would also flag the next record as
# multiallelic, which the .pvar allele-count check below rejects first.
python3 -c "
b = bytearray(open('tmp_data.pgen', 'rb').read())
b[92] = 0x11
open('tmp_bad.pgen', 'wb').write(bytes(b))
"
cp tmp_data.pvar tmp_bad.pvar
cp tmp_data.psam tmp_bad.psam

# --validate has always caught this one.
if $1/plink2 $2 $3 --pfile tmp_bad --validate --out plink2_val 2> tmp_val_err.txt; then
    echo "expected --validate to reject the corrupted .pgen"
    exit 1
fi
grep -q "Invalid 1-bit genotype record" tmp_val_err.txt

# --freq has to reject it too, rather than reporting counts.
if $1/plink2 $2 $3 --pfile tmp_bad --freq --out plink2_bad 2> tmp_err.txt; then
    echo "expected --freq to reject the corrupted .pgen"
    exit 1
fi
grep -q "Failed to unpack" tmp_err.txt
test ! -e plink2_bad.afreq

# A .pgen record with multiallelic hardcalls must not be paired with a .pvar
# line that lists only two alleles for that variant.
#
# plink2 never stores allele counts in the .pgen header, so the .pvar is the
# only source for them.  Before the accompanying fix, a record whose type byte
# says "multiallelic hardcalls present" was read as biallelic when the .pvar
# disagreed, and --make-pgen's MakePgenThread() aborted on an assertion.
# Variant 1 has alleles A,C,G and a 2 in its genotypes; variant 2 is A/C.
for sep in '/' '|'; do
    # '/' gives 4-bit record types, '|' (phased) 8-bit ones.
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' > tmp_ma.vcf
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n' >> tmp_ma.vcf
    printf '1\t10\tv1\tA\tC,G\t.\t.\t.\tGT\t0%s2\t1%s2\t0%s0\t2%s2\n' "$sep" "$sep" "$sep" "$sep" >> tmp_ma.vcf
    printf '1\t20\tv2\tA\tC\t.\t.\t.\tGT\t0%s1\t1%s1\t0%s0\t0%s1\n' "$sep" "$sep" "$sep" "$sep" >> tmp_ma.vcf
    $1/plink2 $2 $3 --vcf tmp_ma.vcf --make-pgen --out tmp_ma
    $1/plink2 $2 $3 --pfile tmp_ma --make-pgen --out plink2_ma_ok

    # v1 and v2 swap ALT columns: the .pvar still has one multiallelic
    # variant, just not the right one.
    awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C" } $3 == "v2" { $5 = "C,G" } { print }' tmp_ma.pvar > tmp_ma_swap.pvar
    # And no multiallelic variant at all.
    awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C" } { print }' tmp_ma.pvar > tmp_ma_bi.pvar
    for p in swap bi; do
        if $1/plink2 $2 $3 --pgen tmp_ma.pgen --pvar tmp_ma_$p.pvar --psam tmp_ma.psam --make-pgen --out plink2_ma_$p 2> tmp_err.txt; then
            echo "expected --make-pgen to reject tmp_ma_$p.pvar"
            exit 1
        fi
        if [ "$p" = bi ] && [ "$sep" = '|' ]; then
            # already caught by the global check
            grep -q "contains multiallelic variants, while .pvar does not" tmp_err.txt
        else
            grep -q "Variant #1 has multiallelic hardcalls in the .pgen file, but only two" tmp_err.txt
        fi
    done
done

# Storage modes 12 and 14 (twelfth header byte) are single-sample encodings.
# A multi-sample file claiming one of them must be rejected by the header
# check instead of tripping an assert (or, with NDEBUG, reading dosage tracks
# far past each 2-6 byte record).
for mode in 12 14; do
    python3 -c "
b = bytearray(open('tmp_data.pgen', 'rb').read())
b[11] = (b[11] & 0xf0) | $mode
open('tmp_mode.pgen', 'wb').write(bytes(b))
"
    cp tmp_data.pvar tmp_mode.pvar
    cp tmp_data.psam tmp_mode.psam
    if $1/plink2 $2 $3 --pfile tmp_mode --validate --out plink2_mode 2> tmp_mode_err.txt; then
        echo "expected --validate to reject storage mode $mode with 400 samples"
        exit 1
    fi
    grep -q "single-sample storage mode" tmp_mode_err.txt
done
