#!/bin/bash

# A zstd-compressed BGEN v1.3 genotype block that decompresses cleanly, but to
# fewer bytes than the block's declared uncompressed length, must be rejected
# with an error.
#
# ZSTD_decompress() returns the number of bytes it wrote.  Both .bgen import
# threads treated any return value other than the declared length as a zstd
# error, and asserted as much; a well-formed frame that is simply too short
# tripped the assertion.

set -exo pipefail

# Dosages are present, so the dosage/phase scan stops at the first variant
# and the conversion pass decodes the remaining ones itself.
$1/plink2 $2 $3 --dummy 8 5 dosage-freq=0.5 --seed 1 --export bgen-1.3 --out tmp_data

# Sanity check: the unmodified file imports.
$1/plink2 $2 $3 --bgen tmp_data.bgen ref-first --sample tmp_data.sample --make-pgen --out plink2_ok

# tmp_bad<N>.bgen: variant N's declared uncompressed length is one byte longer
# than what its zstd frame actually holds.
python3 -c "
import struct
src = open('tmp_data.bgen', 'rb').read()
offset = struct.unpack_from('<I', src, 0)[0]
flags = struct.unpack_from('<I', src, struct.unpack_from('<I', src, 4)[0])[0]
assert flags & 3 == 2, 'expected zstd compression'
variant_ct = struct.unpack_from('<I', src, 8)[0]
for target in (1, variant_ct):
    b = bytearray(src)
    pos = offset + 4
    for vidx in range(1, variant_ct + 1):
        for _ in range(3):
            pos += 2 + struct.unpack_from('<H', b, pos)[0]
        pos += 4
        allele_ct = struct.unpack_from('<H', b, pos)[0]
        pos += 2
        for _ in range(allele_ct):
            pos += 4 + struct.unpack_from('<I', b, pos)[0]
        block_len, uncompressed_len = struct.unpack_from('<II', b, pos)
        if vidx == target:
            struct.pack_into('<I', b, pos + 4, uncompressed_len + 1)
        pos += 4 + block_len
    open('tmp_bad%d.bgen' % target, 'wb').write(bytes(b))
"

# Variant 1 is caught by the dosage/phase scan, variant 5 by the conversion
# pass.  (The error message's vidx is 0-based.)
for v in 1 5; do
    if $1/plink2 $2 $3 --bgen tmp_bad$v.bgen ref-first --sample tmp_data.sample --make-pgen --out plink2_bad$v 2> tmp_err$v.txt; then
        echo "expected --bgen to reject tmp_bad$v.bgen"
        exit 1
    fi
    grep -q "uncompressed_byte_ct mismatch, vidx=$((v - 1))\." tmp_err$v.txt
done
