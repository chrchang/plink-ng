#!/bin/bash

# Small reporting commands with no other test coverage: --sample-counts,
# --pgen-info (including standalone --pgen/--pgi use), and
# --genotyping-rate (including --y-nosex-missing-stats).  Expected values
# are hand-counted.
#
# Regression check: standalone --pgen-info reported a maximum allele count of
# 2 for a multiallelic .pgen with neither phase nor dosage data, since
# variant record types were only scanned for the multiallelic bit when they
# were stored with 8 bits.

set -exo pipefail

plink2="$1/plink2 $2 $3"
pgen_compress="$1/pgen_compress"

# ts1: transition; tv1: transversion; ind1: indel; sym1: symbolic; ma1:
# triallelic (C>T transition, C>G transversion).  U1 has unknown sex, so it
# is treated as female.  On chrY, M1 and F1 are missing and U1 has a
# haploid ALT call.
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##contig=<ID=X,length=200000000>
##contig=<ID=Y,length=60000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M1	F1	U1
1	100	ts1	A	G	.	.	.	GT	0/1	1/1	0/0
1	200	tv1	A	C	.	.	.	GT	1/1	0/1	./.
1	300	ind1	AT	A	.	.	.	GT	0/1	0/0	0/0
1	400	sym1	A	<DEL>	.	.	.	GT	0/0	0/1	0/0
1	500	ma1	C	T,G	.	.	.	GT	1/2	2/2	0/1
X	3000000	x1	G	A	.	.	.	GT	1	0/1	0/0
Y	3000000	y1	C	T	.	.	.	GT	.	.	1
VCF
printf '#IID\tSEX\nM1\t1\nF1\t2\nU1\tNA\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

# --sample-counts, default columns.
# - A heterozygous ALT1/ALT2 call counts each ALT separately in the SNP
#   subtype columns, so M1's ma1 call adds 2 to HET_SNP_CT and 1 each to the
#   transition and transversion counts; a homozygous-ALT call adds 1.
# - ind1 is a singleton for M1, sym1 for F1.  x1's ALT appears in M1
#   (haploid) and F1 (diploid): that is a diploid singleton for F1 only.
# - U1's chrY call counts in the INCL_FEMALE_Y columns.
$plink2 --pfile tmp_data --sample-counts --out tmp_sc > /dev/null
diff tmp_sc.scount - <<'EOF'
#IID	HOM_REF_CT	HOM_ALT_SNP_CT	HET_SNP_CT	DIPLOID_TRANSITION_CT	DIPLOID_TRANSVERSION_CT	DIPLOID_NONSNP_NONSYMBOLIC_CT	DIPLOID_SINGLETON_CT	HAP_REF_INCL_FEMALE_Y_CT	HAP_ALT_INCL_FEMALE_Y_CT	MISSING_INCL_FEMALE_Y_CT
M1	1	1	3	2	2	1	1	0	1	1
F1	1	2	2	2	2	0	2	0	0	1
U1	4	0	1	1	0	0	0	0	1	1
EOF

# Every other column.  Haploid calls now count toward the transition
# column, x1 is no longer a singleton (2 carriers), and chrY is excluded
# for F1 and U1.
$plink2 --pfile tmp_data --sample-counts cols=fid,sid,sex,hom,homref,homalt,het,refalt,het2alt,ts,tv,nonsnpsymb,symbolic,nonsnp,single,hapref,hapalt,missing --out tmp_sc_all > /dev/null
diff tmp_sc_all.scount - <<'EOF'
#FID	IID	SID	SEX	HOM_CT	HOM_REF_CT	HOM_ALT_CT	HET_CT	HET_REF_ALT_CT	HET_2ALT_CT	TRANSITION_CT	TRANSVERSION_CT	NONSNP_NONSYMBOLIC_CT	SYMBOLIC_CT	NONSNP_CT	SINGLETON_CT	HAP_REF_CT	HAP_ALT_CT	MISSING_CT
0	M1	0	1	2	1	1	3	2	1	3	2	1	0	1	1	0	1	1
0	F1	0	2	3	1	2	3	3	0	2	2	0	1	1	1	0	0	0
0	U1	0	NA	4	4	0	1	1	0	1	0	0	0	0	0	0	0	1
EOF

# 'zs' and a +/- column descriptor.
$plink2 --pfile tmp_data --sample-counts zs cols=-missingwfemaley,+missing,-homref --out tmp_sc_zs > /dev/null
zstd -dc tmp_sc_zs.scount.zst | diff - <(paste <(cut -f 1,3-10 tmp_sc.scount) <(cut -f 19 tmp_sc_all.scount))
if $plink2 --pfile tmp_data --sample-counts cols=hom,+het --out tmp_sc_bad > /dev/null 2>&1; then
    exit 1
fi

# --genotyping-rate weights every variant equally; on chrY, only males count
# by default, and unknown-sex samples also count with
# --y-nosex-missing-stats.
# - Default: (17/3 [non-Y] + 0/1 [chrY]) / 7 = 0.809524.
# - --y-nosex-missing-stats: (17/3 + 1/2) / 7 = 0.880952.
# - Males only: (6/1 + 0/1) / 7 = 0.857143.
# - Females only: chrY is skipped, and F1 has no other missing calls.
# - chr1 only: 14/15.
# $1: expected rate line suffix; remaining args: plink2 flags.
check_rate() {
    local expected=$1
    shift
    $plink2 --pfile tmp_data "$@" --genotyping-rate --out tmp_gr > /dev/null
    grep -qxF "Total (hardcall) genotyping rate $expected" tmp_gr.log
}
check_rate "is 0.809524."
check_rate "is 0.880952." --y-nosex-missing-stats
check_rate "in remaining samples is 0.857143." --keep-males
check_rate "in remaining samples is exactly 1." --keep-females
check_rate "is 0.933333." --chr 1

# Prints the --pgen-info report from log file $1.
pgen_info_report() {
    awk '/^--pgen-info on/ {f = 1} /^$/ || /^End time/ {f = 0} f' "$1"
}

# --pgen-info on the full fileset.
# $1: output prefix; remaining args: plink2 flags.  Prints the report.
pgen_info() {
    local out=$1
    shift
    $plink2 "$@" --pgen-info --out $out > /dev/null
    pgen_info_report $out.log
}
diff <(pgen_info tmp_pi --pfile tmp_data) - <<'EOF'
--pgen-info on tmp_data.pgen:
  Variants: 7
  Samples: 3
  REF alleles are all known
  Maximum allele count for a single variant: 3
  No hardcalls are explicitly phased
  No dosages present
EOF

# Standalone (.pgen only): the maximum allele count isn't stored, but the
# multiallelic variant is still detectable.
diff <(pgen_info tmp_pi_sa --pgen tmp_data.pgen) - <<'EOF'
--pgen-info on tmp_data.pgen:
  Variants: 7
  Samples: 3
  REF alleles are all known
  Maximum allele count for a single variant: >2, not explicitly stored
  No hardcalls are explicitly phased
  No dosages present
EOF
$plink2 --pfile tmp_data --chr 1 --snps ts1-ind1 --make-pgen --out tmp_bi > /dev/null
pgen_info tmp_pi_bi --pgen tmp_bi.pgen | grep -qxF '  Maximum allele count for a single variant: 2'

# Phased hardcalls, unphased dosages, and a mix of known and provisional REF
# alleles.  S1's v2 dosage (0.25) has no hardcall.
cat > tmp_dosage.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
1	100	v1	A	G	.	.	PR	GT:DS	0|1:1	1|1:2	0/0:0	0/1:1
1	200	v2	A	C	.	.	.	GT:DS	0/0:0.25	0/1:1	./.:.	1/1:2
1	300	v3	A	T	.	.	.	GT:DS	0/0:0	0/1:1	1/1:2	./.:.
VCF
$plink2 --vcf tmp_dosage.vcf dosage=DS --make-pgen vzs --out tmp_dosage > /dev/null
diff <(pgen_info tmp_pi_d --pfile tmp_dosage vzs) - <<'EOF'
--pgen-info on tmp_dosage.pgen:
  Variants: 3
  Samples: 4
  REF alleles are a mix of known and provisional
  Maximum allele count for a single variant: 2
  Explicitly phased hardcalls present
  Dosage present, none explicitly phased
EOF
# Dosage may be multiallelic even when hardcalls aren't.
pgen_info tmp_pi_dsa --pgen tmp_dosage.pgen | grep -qxF '  Maximum allele count for a single variant: not explicitly stored'
# 3 of 12 hardcalls are missing, but only 2 of 12 dosages.
$plink2 --pfile tmp_dosage vzs --genotyping-rate --out tmp_grd > /dev/null
grep -qxF 'Total (hardcall) genotyping rate is 0.75.' tmp_grd.log
$plink2 --pfile tmp_dosage vzs --genotyping-rate dosage --out tmp_grd2 > /dev/null
grep -qxF 'Total (dosage) genotyping rate is 0.833333.' tmp_grd2.log

# Phased dosages.
cat > tmp_hds.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haploid dosages">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	v1	A	G	.	.	.	GT:HDS	0|1:0.1,0.9	1|1:1,1
VCF
$plink2 --vcf tmp_hds.vcf dosage=HDS --make-pgen --out tmp_hds > /dev/null
pgen_info tmp_pi_h --pfile tmp_hds | grep -qxF '  Explicitly phased dosages present'

# A PLINK 1 fileset stores all REF alleles as provisional; --pgen-info
# skips a .bed (with a warning), but not a .pgen written in a separate-index
# layout.  --pgi points to a renamed index file.
$plink2 --pfile tmp_bi --make-bed --out tmp_bi > /dev/null
pgen_info tmp_pi_bed --bfile tmp_bi > /dev/null
grep -q '^Warning: Skipping --pgen-info since a .bed file was provided' tmp_pi_bed.log
$pgen_compress -i tmp_bi.bed tmp_si.pgen 3 > /dev/null
test -s tmp_si.pgen.pgi
diff <(pgen_info tmp_pi_si --pgen tmp_si.pgen) - <<'EOF'
--pgen-info on tmp_si.pgen:
  Variants: 3
  Samples: 3
  REF alleles are all provisional
  Maximum allele count for a single variant: 2
  No hardcalls are explicitly phased
  No dosages present
EOF
mv tmp_si.pgen.pgi tmp_si_renamed.pgi
if $plink2 --pgen tmp_si.pgen --pgen-info --out tmp_pi_nopgi > /dev/null 2>&1; then
    exit 1
fi
# Same report as with the default index name.
diff <(pgen_info tmp_pi_pgi --pgen tmp_si.pgen --pgi tmp_si_renamed.pgi) <(pgen_info_report tmp_pi_si.log)
# --pgi also works for regular loading: frequencies match the .bed's.
cp tmp_bi.bim tmp_si.pvar
cp tmp_bi.fam tmp_si.psam
$plink2 --pgen tmp_si.pgen --pgi tmp_si_renamed.pgi --pvar tmp_si.pvar --psam tmp_si.psam --freq --out tmp_si > /dev/null
$plink2 --bfile tmp_bi --freq --out tmp_bi > /dev/null
diff tmp_si.afreq tmp_bi.afreq
