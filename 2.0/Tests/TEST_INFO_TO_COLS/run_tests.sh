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
1	600	rs6	A	G	.	.	DP=1;DP=2;;AF=0.5;	GT	0/1	0/0	0/0
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

# 6. Without filters, the variant file is scanned directly.  That covers .pvar,
#    zstd-compressed .pvar, gzipped VCF, and CRLF input, and must agree with
#    the regular path, which a no-op --chr filter forces.
$1/plink2 $2 $3 --pvar tmp_data.pvar --info-to-cols all --out plink2_frompvar
diff -q plink2_frompvar.vinfo plink2_all.vinfo
$1/plink2 $2 $3 --pfile tmp_data --make-just-pvar zs --out tmp_zs
$1/plink2 $2 $3 --pvar tmp_zs.pvar.zst --info-to-cols all --out plink2_frompvarzst
diff -q plink2_frompvarzst.vinfo plink2_all.vinfo
gzip -c tmp_in.vcf > tmp_in.vcf.gz
$1/plink2 $2 $3 --vcf tmp_in.vcf.gz --info-to-cols all --out plink2_fromvcfgz
diff -q plink2_fromvcfgz.vinfo plink2_all.vinfo
awk '{printf "%s\r\n", $0}' tmp_in.vcf > tmp_crlf.vcf
$1/plink2 $2 $3 --vcf tmp_crlf.vcf --info-to-cols all --out plink2_fromcrlf
diff -q plink2_fromcrlf.vinfo plink2_all.vinfo
$1/plink2 $2 $3 --vcf tmp_in.vcf --chr 1 --info-to-cols all --out plink2_regular
diff -q plink2_regular.vinfo plink2_all.vinfo
# VCF prohibits duplicate INFO keys; if one appears anyway, its first value is
# reported.
test "$(awk '$3 == "rs6" {print $6}' plink2_all.vinfo)" = "1"

# 7. Values and alleles too long for the direct-write fast case.
awk 'BEGIN {
  s = "x";
  while (length(s) < 200000) s = s s;
  print "##fileformat=VCFv4.2";
  print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">";
  print "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type\">";
  print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1";
  print "1\t100\trs1\tA\tG\t.\t.\tDP=1;VT=" s "\tGT\t0/1";
  print "1\t200\trs2\tA" s "\tG,C" s "\t.\t.\tVT=" s s ";DP=2\tGT\t0/1";
  print "1\t300\trs3\tA\tG\t.\t.\tDP=3\tGT\t0/1";
}' > tmp_long.vcf
$1/plink2 $2 $3 --vcf tmp_long.vcf --info-to-cols DP,VT --out plink2_long
python3 check_vinfo.py tmp_long.vcf plink2_long.vinfo
$1/plink2 $2 $3 --vcf tmp_long.vcf --chr 1 --info-to-cols DP,VT --out plink2_long_regular
diff -q plink2_long_regular.vinfo plink2_long.vinfo

# 8. Variant filters apply.
$1/plink2 $2 $3 --pfile tmp_data --extract-if-info "DP > 20" --info-to-cols DP --out plink2_filtered
test "$(awk 'NR > 1' plink2_filtered.vinfo | wc -l)" -eq 1

# 9. A fileset with no INFO column is an error rather than an empty table.
#    PLINK 1 binary filesets have no INFO at all, so --make-bed gives one.
$1/plink2 $2 $3 --pfile tmp_data --max-alleles 2 --make-bed --out tmp_noinfo
if $1/plink2 $2 $3 --bfile tmp_noinfo --info-to-cols DP --out plink2_noinfo 2> tmp_err.txt; then
    echo "expected --info-to-cols to fail without an INFO column"
    exit 1
fi
grep -q "requires an INFO column" tmp_err.txt
