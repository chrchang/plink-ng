#!/bin/bash

# --validate must reject a sparse (difflist) genotype record whose list
# contains a sample in the record's common genotype category.
#
# The .pgen spec requires every difflist entry to change its sample's
# category.  --validate used to accept this entry, and sparse readers such as
# SampleCountsThread(), which index arrays with (entry ^ common) - 1, then
# segfaulted on a file that had passed it.

set -exo pipefail

# 400 samples, so that the rare variant is stored as a difflist.
python3 -c "
n = 400
print('##fileformat=VCFv4.3')
print('##contig=<ID=1,length=1000>')
print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join('s%d' % i for i in range(n)))
print('1\t10\tv10\tA\tC\t.\t.\t.\tGT\t' + '\t'.join('0/1' if i in (5, 300) else '0/0' for i in range(n)))
print('1\t20\tv20\tA\tC\t.\t.\t.\tGT\t' + '\t'.join(['0/1', '1/1', '0/0', '0/0'][i % 4] for i in range(n)))
" > tmp_data.vcf
$1/plink2 $2 $3 --vcf tmp_data.vcf --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --validate --out plink2_val_ok

# v10 is stored as a difflist (record type 4, common genotype hom-ref): a
# 1-byte length (2), s5's 2-byte sample index, then one byte of 2-bit
# genotypes (0x05: het, het).  Change s5's entry to hom-ref, the common value.
python3 -c "
import struct
d = bytearray(open('tmp_data.pgen', 'rb').read())
assert d[2] == 0x10 and struct.unpack_from('<II', d, 3) == (2, 400)
assert (d[11] & 15) < 4 and (d[20] & 15) == 4, 'expected a difflist record'
rec = struct.unpack_from('<Q', d, 12)[0]
assert d[rec:rec + 4] == bytes([2, 5, 0, 5])
d[rec + 3] = 4
open('tmp_common.pgen', 'wb').write(bytes(d))
"
cp tmp_data.pvar tmp_common.pvar
cp tmp_data.psam tmp_common.psam
if $1/plink2 $2 $3 --pfile tmp_common --validate --out plink2_val_common; then
    echo "--validate accepted a difflist entry equal to the common genotype"
    exit 1
fi
grep -q "Invalid genotype difflist for (0-based) variant #0" plink2_val_common.log
