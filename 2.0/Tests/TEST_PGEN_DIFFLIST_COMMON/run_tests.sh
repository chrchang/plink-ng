#!/bin/bash

# A sparse (difflist) genotype record whose list contains a sample in the
# record's common genotype category must be read like the same data stored
# normally.
#
# The .pgen spec describes these lists as covering samples outside the common
# category, but --validate accepts such an entry, and a full decode just
# writes the common value back.  PgrGetDifflistOrGenovec() handed it to its
# sparse callers as-is, and SampleCountsThread() indexes its count arrays with
# (entry ^ common) - 1: --sample-counts segfaulted.

set -exo pipefail

# 400 samples, so that the reader returns the rare variant in sparse form.
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
$1/plink2 $2 $3 --pfile tmp_common --validate --out plink2_val

# Reference: the same genotypes written normally.
$1/plink2 $2 $3 --pfile tmp_common --make-pgen --out tmp_redo
printf 's5\ns6\ns300\ns301\n' > tmp_keep.txt
for x in common redo; do
    $1/plink2 $2 $3 --pfile tmp_$x --sample-counts --out plink2_$x
    $1/plink2 $2 $3 --pfile tmp_$x --keep tmp_keep.txt --sample-counts --out plink2_keep_$x
done
diff plink2_common.scount plink2_redo.scount
diff plink2_keep_common.scount plink2_keep_redo.scount
