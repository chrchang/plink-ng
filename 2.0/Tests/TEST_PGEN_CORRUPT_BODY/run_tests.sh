#!/bin/bash

# Corrupted .pgen record bodies must be reported as malformed, not read out of
# bounds, by commands that don't run --validate.
#
# Each case below changes a few bytes of a valid fileset.  --validate has
# always rejected all of them; before the accompanying fixes, the listed
# command segfaulted, hit an assertion, or exited 0 with wrong results.
# The fixes are per-record or per-load checks; see the commit message.

set -exo pipefail

PLINK2=("$1/plink2" $2 $3)

python3 make_vcfs.py
"${PLINK2[@]}" --vcf tmp_hc.vcf --make-pgen --out tmp_hc
"${PLINK2[@]}" --vcf tmp_ds.vcf dosage=HDS --make-pgen --out tmp_ds

# The byte offsets below assume these record types and lengths; if the writer
# changes, they need updating.
python3 -c "
import struct
for name, want in (('hc', [(0x11, 8), (0x04, 4), (0x08, 11), (0x18, 16), (0x09, 8)]), ('ds', [(0x40, 54), (0xd1, 103)])):
    d = open('tmp_' + name + '.pgen', 'rb').read()
    vct = struct.unpack_from('<I', d, 3)[0]
    got = [(d[20 + i], d[20 + vct + i]) for i in range(vct)]
    assert got == want, (name, got)
"

# 12 of the 24 samples, for the sample-subsetting code paths.
awk 'NR > 1 && NR % 2 == 0 { print $1 }' tmp_hc.psam > tmp_keep.txt

run_case() {
    # usage: run_case <name> <base fileset> '<edits>' '<command>' <expected message>
    local name=$1 base=$2 edits=$3 cmd=$4 msg=$5
    python3 pgen_edit.py tmp_$base.pgen tmp_$name.pgen $edits
    cp tmp_$base.pvar tmp_$name.pvar
    cp tmp_$base.psam tmp_$name.psam
    if "${PLINK2[@]}" --pfile tmp_$name --validate --out plink2_${name}_val > /dev/null; then
        echo "$name: expected --validate to reject the edited .pgen"
        exit 1
    fi
    local rc=0
    "${PLINK2[@]}" --pfile tmp_$name $cmd --out plink2_$name 2> tmp_${name}_err.txt || rc=$?
    # 0 means the corruption went unnoticed; > 128 means a signal.
    if [ $rc -eq 0 ] || [ $rc -gt 128 ]; then
        echo "$name: expected a clean error, got exit code $rc"
        exit 1
    fi
    grep -q "$msg" tmp_${name}_err.txt
}

# Header: an LD-compressed first record makes GetLdbaseVidx() search for its
# reference variant before the start of the record-type array.
run_case ld_first hc 'vrtype:0=0x02' '--freq' 'is LD-compressed'
# Header: a phased-dosage flag without the dosage flags.
run_case dphase_no_dosage ds 'vrtype:1=0x91' '--freq' 'phased dosage bit set'
# Difflist: the second sample index is out of range.  64-bit builds checked
# only at 64-entry group boundaries, after IsSet() had already used it.
run_case difflist_index hc 'byte:1:3=0x7f' "--keep tmp_keep.txt --freq" 'Failed to unpack'
# Multiallelic hardcalls, bitarray form: empty subset, and a set bit past the
# end of the subset.  Both made the popcount disagree with the bits expanded.
run_case aux1a_empty hc 'byte:2:7=0 byte:2:8=0' "--keep tmp_keep.txt --make-pgen" 'Failed to unpack'
run_case aux1a_trailing hc 'byte:2:8=0x80' '--make-pgen' 'Failed to unpack'
# Allele codes out of range: ALT4 in a 5-allele het-ref track (codes 0..2),
# and ALT1/ALT1 in the ALTx/ALTy track.  These index count arrays and patch
# consumers.
run_case aux1a_code hc 'byte:3:9=0x1b' '--export vcf' 'Failed to unpack'
run_case aux1b_code hc 'byte:3:11=0x00' '--export vcf' 'Failed to unpack'
# Sample-list form of the ref/ALT2 track (v500): the one entry names sample
# 4, which is homozygous ref.  Callers pair list lengths with genotype counts,
# and a list longer than the matching genotype count made --hardy count
# about 4 billion genotypes.
run_case aux1a_list_sample hc 'byte:4:7=0x04' '--export vcf' 'Failed to unpack'
# Hardcall phase: a set bit past the first het_ct + 1 bits.
run_case aux2_trailing hc 'byte:3:13=0xff' '--make-pgen' 'Failed to unpack'
# Dosage out of range (0x9000 > 32768), and a phased-dosage delta that puts a
# haplotype dosage outside [0, 1].  The VCF printers index digit tables with
# these values.
run_case dosage_value ds 'u16:0:6=0x9000' '--export vcf vcf-dosage=HDS-force' 'Failed to unpack'
run_case dphase_value ds 'u16:1:61=0x7000' '--export vcf vcf-dosage=HDS-force' 'Failed to unpack'
# The same bad dosage through --r2: the dosage readers now reject it, and
# PgrGetInv1D() used to invert dosage_main before looking at the error code,
# with dosage_ct never set.
run_case dosage_value_r2 ds 'u16:0:6=0x9000' '--r2-unphased --ld-window-r2 0' 'Failed to unpack'
