#!/bin/bash

# --info-to-cols, checked against the source VCF rather than against another
# plink2 run.

set -exo pipefail

cat > tmp_in.vcf <<'VCFEOF'
##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">
##INFO=<ID=UNUSED,Number=1,Type=String,Description="Never present">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	100	rs1	A	G	.	.	DP=30;AF=0.25;DB;VT=SNP	GT	0/0	0/1	1/1
1	200	rs2	C	T	.	.	DP=12;VT=SNP	GT	0/1	0/1	0/0
1	300	rs3	G	A,C	.	.	AF=0.1,0.2;DB	GT	0/1	1/2	0/0
1	400	rs4	T	A	.	.	.	GT	0/0	0/0	0/1
1	500	rs5	C	G	.	.	VT	GT	0/1	0/0	0/0
VCFEOF

$1/plink2 $2 $3 --vcf tmp_in.vcf --make-pgen --out tmp_data

# 1. Explicit keys, in the order asked for.
$1/plink2 $2 $3 --pfile tmp_data --info-to-cols DP,AF,DB,VT --out plink2_keys
head -n 1 plink2_keys.vinfo | grep -qx '#CHROM	POS	ID	REF	ALT	DP	AF	DB	VT'
python3 check_vinfo.py tmp_in.vcf plink2_keys.vinfo

# 2. 'all' takes every declared key, including one no variant carries.
$1/plink2 $2 $3 --pfile tmp_data --info-to-cols all --out plink2_all
head -n 1 plink2_all.vinfo | grep -qx '#CHROM	POS	ID	REF	ALT	DP	AF	DB	VT	UNUSED'
python3 check_vinfo.py tmp_in.vcf plink2_all.vinfo
# The key no variant carries must be all-NA rather than absent.
test "$(awk 'NR > 1 {print $10}' plink2_all.vinfo | sort -u)" = "NA"

# 3. A key that is not declared at all still works, reported by presence.
$1/plink2 $2 $3 --pfile tmp_data --info-to-cols DP,NOSUCHKEY --out plink2_undeclared
test "$(awk 'NR > 1 {print $7}' plink2_undeclared.vinfo | sort -u)" = "NA"

# 4. Zstd output round-trips.
$1/plink2 $2 $3 --pfile tmp_data --info-to-cols all zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.vinfo.zst > plink2_zs.vinfo
diff -q plink2_zs.vinfo plink2_all.vinfo

# 5. Reading the VCF directly gives the same table as reading the .pgen.
$1/plink2 $2 $3 --vcf tmp_in.vcf --info-to-cols all --out plink2_fromvcf
diff -q plink2_fromvcf.vinfo plink2_all.vinfo

# 6. Variant filters apply.
$1/plink2 $2 $3 --pfile tmp_data --extract-if-info "DP > 20" --info-to-cols DP --out plink2_filtered
test "$(awk 'NR > 1' plink2_filtered.vinfo | wc -l)" -eq 1

# 7. A fileset with no INFO column is an error rather than an empty table.
#    PLINK 1 binary filesets have no INFO at all, so --make-bed gives one.
$1/plink2 $2 $3 --pfile tmp_data --max-alleles 2 --make-bed --out tmp_noinfo
if $1/plink2 $2 $3 --bfile tmp_noinfo --info-to-cols DP --out plink2_noinfo 2> tmp_err.txt; then
    echo "expected --info-to-cols to fail without an INFO column"
    exit 1
fi
grep -q "requires an INFO column" tmp_err.txt
