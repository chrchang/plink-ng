#!/bin/bash

# A .pgen record may carry an empty dosage list (record type 0x20 with a
# zero-length sample list): the spec allows it and --validate accepts it.
# ParseAndSaveDeltalistAsBitarr() returned before clearing the dosage-present
# bitarray in that case, and the sample-subsetting path of ParseDosage16()
# popcounts that bitarray without checking the list length first.  After a
# variant with dosages, the stale bits made it copy dosages that don't exist,
# walking off the end of the array: --export vcf with --keep segfaulted.
#
# Such a record must read exactly like the same hardcalls without a dosage
# track.

set -exo pipefail

# v1: every sample has a dosage (fixed-width encoding), which fills the
# reader's dosage-present workspace.  v2: hardcalls only.
cat > tmp_data.vcf <<VCF
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6	s7	s8
1	10	v1	A	C	.	.	.	GT:DS	0/0:0.1	0/1:0.9	0/1:1.1	1/1:1.8	0/0:0.2	0/1:0.8	0/0:0.3	1/1:1.7
1	20	v2	A	G	.	.	.	GT	0/0	0/1	1/1	./.	0/1	0/0	0/0	0/1
VCF
$1/plink2 $2 $3 --vcf tmp_data.vcf dosage=DS --make-pgen --out tmp_data

# Rewrite v2 as an empty dosage list: record type 0x20, then 2-bit hardcalls
# for the 8 samples, then a zero-length list.
python3 -c "
import struct
d = open('tmp_data.pgen', 'rb').read()
assert d[2] == 0x10 and struct.unpack_from('<II', d, 3) == (2, 8)
ctrl = d[11]
assert (ctrl & 15) == 4, 'expected 8-bit record types, 1-byte record lengths'
fpos = struct.unpack_from('<Q', d, 12)[0]
vrtypes = bytearray(d[20:22])
lens = list(d[22:24])
assert fpos == 24 and vrtypes[0] == 0x40
rec0 = d[fpos:fpos + lens[0]]
geno = [0, 1, 2, 3, 1, 0, 0, 1]
genovec = bytes(sum(g << (2 * k) for k, g in enumerate(geno[i:i + 4])) for i in (0, 4))
rec1 = genovec + b'\x00'
vrtypes[1] = 0x20
out = d[:20] + bytes(vrtypes) + bytes([len(rec0), len(rec1)]) + rec0 + rec1
open('tmp_empty.pgen', 'wb').write(out)
"
cp tmp_data.pvar tmp_empty.pvar
cp tmp_data.psam tmp_empty.psam
$1/plink2 $2 $3 --pfile tmp_empty --validate --out plink2_val

printf 's2\ns3\ns5\ns8\n' > tmp_keep.txt
for x in data empty; do
    $1/plink2 $2 $3 --pfile tmp_$x --export vcf vcf-dosage=DS-force --out plink2_full_$x
    $1/plink2 $2 $3 --pfile tmp_$x --keep tmp_keep.txt --export vcf vcf-dosage=DS-force --out plink2_keep_$x
    $1/plink2 $2 $3 --pfile tmp_$x --keep tmp_keep.txt --freq --out plink2_keep_$x
done
for x in full keep; do
    diff <(grep -v '^##' plink2_${x}_data.vcf) <(grep -v '^##' plink2_${x}_empty.vcf)
done
diff plink2_keep_data.afreq plink2_keep_empty.afreq
