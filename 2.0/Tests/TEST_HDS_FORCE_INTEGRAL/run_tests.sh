#!/bin/bash

# --export vcf vcf-dosage=HDS-force must handle an unphased dosage of exactly 0
# or 2.
#
# plink2's importer stores such a dosage as a plain hardcall, but the .pgen
# spec permits it in the dosage track (it is required in the fixed-width
# encoding, and appears there whenever the hardcall is missing), and
# --validate accepts it.  PrintHaploidNonintDosage() assumed its argument was
# strictly between 0 and 32768: 32768 indexed past the end of
# u32toa_trunc4()'s digit table, and 0 printed "0." instead of "0".

set -exo pipefail

# Every sample has a dosage, so the importer uses the fixed-width dosage
# encoding (record type 0x40).  s2's hardcall is missing.
cat > tmp_data.vcf <<VCF
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4
1	10	v1	A	C	.	.	.	GT:DS	0/0:0.1	./.:1.3	0/1:0.9	1/1:1.8
VCF
$1/plink2 $2 $3 --vcf tmp_data.vcf dosage=DS --make-pgen --out tmp_data

# Set s2's stored dosage to exactly 2 (32768) in tmp_two, and to exactly 0 in
# tmp_zero.
python3 -c "
import struct
d = bytearray(open('tmp_data.pgen', 'rb').read())
assert d[2] == 0x10 and struct.unpack_from('<II', d, 3) == (1, 4)
assert d[20] == 0x40, 'expected a fixed-width dosage record'
rec = struct.unpack_from('<Q', d, 12)[0]
# 1 byte of 2-bit hardcalls for 4 samples, then one uint16 dosage per sample
for name, val in (('two', 32768), ('zero', 0)):
    struct.pack_into('<H', d, rec + 1 + 2 * 1, val)
    open('tmp_' + name + '.pgen', 'wb').write(bytes(d))
"
for x in two zero; do
    cp tmp_data.pvar tmp_$x.pvar
    cp tmp_data.psam tmp_$x.psam
    $1/plink2 $2 $3 --pfile tmp_$x --validate --out plink2_val_$x
    $1/plink2 $2 $3 --pfile tmp_$x --export vcf vcf-dosage=HDS-force --out plink2_$x
done
test "$(grep -v '^#' plink2_two.vcf | cut -f 11)" = "./.:2:1,1"
test "$(grep -v '^#' plink2_zero.vcf | cut -f 11)" = "./.:0:0,0"
