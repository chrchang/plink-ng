#!/bin/bash

# Sample-ID and phenotype-loading flags with no other test coverage:
# --select-sid-representatives, --indv, --iid-sid, --strict-sid0,
# --not-pheno, --not-covar, --mpheno, --family, --family-missing-catname,
# --within/--mwithin, --missing-catname, --no-categorical,
# --neg9-pheno-really-missing, and --no-input-missing-phenotype.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Prints "FID:IID:SID" for each sample in .psam file $1 (missing FID/SID
# columns are printed as empty).
psam_ids() {
    awk -F '\t' 'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
        {printf "%s:%s:%s ", ("#FID" in col)? $col["#FID"] : "", $col[("#FID" in col)? "IID" : "#IID"], ("SID" in col)? $col["SID"] : ""}' "$1"
    echo
}

# $1: expected IDs, as printed by psam_ids; remaining args: plink2 flags.
check_samples() {
    local expected=$1
    shift
    rm -f tmp_s.psam
    $plink2 "$@" --make-just-psam --out tmp_s > /dev/null
    test "$(psam_ids tmp_s.psam)" = "$expected "
}

# --- SID-based flags ---
# FID+IID groups with multiple SIDs: (f1, i1) with s1/s2/s10 (all male),
# (f3, p3) with x/y, and (f3, c3) with u/v (p3's children).  Missingness:
# - s1: v1's dosage (0.5) has no hardcall.
# - s2: heterozygous haploid call on chrX.
# - s10: one missing call.
# - x and v: one missing call each.
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##contig=<ID=X,length=200000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	a	b	c	d	e	f	g	h
1	100	v1	A	G	.	.	.	GT:DS	0/1:0.5	0/1:1	0/0:0	0/0:0	0/1:1	0/0:0	1/1:2	0/1:1
1	200	v2	A	C	.	.	.	GT:DS	0/0:0	0/0:0	./.:.	0/1:1	0/0:0	0/0:0	0/0:0	0/0:0
1	300	v3	A	T	.	.	.	GT:DS	0/0:0	0/1:1	0/0:0	./.:.	0/0:0	0/1:1	0/0:0	1/1:2
1	400	v4	C	T	.	.	.	GT:DS	1/1:2	0/0:0	0/0:0	0/0:0	0/0:0	0/0:0	./.:.	0/0:0
X	3000000	x1	G	A	.	.	.	GT:DS	0:0	0/1:1	0:0	0:0	1:1	0/1:1	0/0:0	0/1:1
VCF
printf '#IID\tSEX\na\t1\nb\t1\nc\t1\nd\t1\ne\t1\nf\t2\ng\t2\nh\t2\n' > tmp_sex.psam
$plink2 --vcf tmp_data.vcf dosage=DS --psam tmp_sex.psam --make-pgen --out tmp_sid > /dev/null
printf '#FID\tIID\tSID\tPAT\tMAT\tSEX\nf1\ti1\ts1\t0\t0\t1\nf1\ti1\ts2\t0\t0\t1\nf1\ti1\ts10\t0\t0\t1\nf3\tp3\tx\t0\t0\t1\nf3\tp3\ty\t0\t0\t1\nf3\tc3\tu\tp3\t0\t2\nf3\tc3\tv\tp3\t0\t2\nf4\ti4\t0\t0\t0\t2\n' > tmp_sid.psam

others="f3:p3:y f3:c3:u f4:i4:0"
# Default 'hc' mode: s1's dosage-only call counts as missing.
check_samples "f1:i1:s2 $others" --pfile tmp_sid --select-sid-representatives first
check_samples "f1:i1:s2 $others" --pfile tmp_sid --select-sid-representatives hc last
# 'hh-missing': all three (f1, i1) samples tie, so the SID tiebreak decides.
# Natural sort: s1 < s2 < s10; ASCII sort: s1 < s10 < s2.
check_samples "f1:i1:s1 $others" --pfile tmp_sid --select-sid-representatives hh-missing first
check_samples "f1:i1:s10 $others" --pfile tmp_sid --select-sid-representatives hh-missing last
check_samples "f1:i1:s1 $others" --pfile tmp_sid --select-sid-representatives hh-missing first-ascii
check_samples "f1:i1:s2 $others" --pfile tmp_sid --select-sid-representatives hardcall last-ascii
check_samples "f1:i1:s2 $others" --pfile tmp_sid --select-sid-representatives hh-missing last-ascii
# 'dosage': s1 and s2 tie.
check_samples "f1:i1:s1 $others" --pfile tmp_sid --select-sid-representatives dosage first
check_samples "f1:i1:s2 $others" --pfile tmp_sid --select-sid-representatives dosage last
# 'sid-only' ignores missingness.
check_samples "f1:i1:s1 f3:p3:x f3:c3:u f4:i4:0" --pfile tmp_sid --select-sid-representatives sid-only first
check_samples "f1:i1:s10 f3:p3:y f3:c3:v f4:i4:0" --pfile tmp_sid --select-sid-representatives sid-only last
check_samples "f1:i1:s2 f3:p3:y f3:c3:v f4:i4:0" --pfile tmp_sid --select-sid-representatives sid-only last-ascii
# 'parents-only': only p3 is a parent.
check_samples "f1:i1:s1 f1:i1:s2 f1:i1:s10 f3:p3:y f3:c3:u f3:c3:v f4:i4:0" --pfile tmp_sid --select-sid-representatives first parents-only
check_samples "f1:i1:s1 f1:i1:s2 f1:i1:s10 f3:p3:x f3:c3:u f3:c3:v f4:i4:0" --pfile tmp_sid --select-sid-representatives sid-only first parents-only
if $plink2 --pfile tmp_sid --select-sid-representatives hc --make-just-psam --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# --indv: 3 parts are FID-IID-SID; 1 part is an IID with FID 0.  Without
# a SID, all SIDs match (as with a --keep file without a SID column).
check_samples "f1:i1:s1 f1:i1:s2 f1:i1:s10" --pfile tmp_sid --indv f1 i1
check_samples "f3:c3:u f3:c3:v" --pfile tmp_sid --indv f3 c3
printf 'f1\ti1\n' > tmp_keep_nosid.txt
check_samples "f1:i1:s1 f1:i1:s2 f1:i1:s10" --pfile tmp_sid --keep tmp_keep_nosid.txt
check_samples "f1:i1:s10" --pfile tmp_sid --indv f1 i1 s10
check_samples "f3:c3:v" --pfile tmp_sid --indv f3 c3 v
# (The SID column is dropped when every remaining SID is 0.)
check_samples "f4:i4:" --pfile tmp_sid --indv f4 i4 0
if $plink2 --pfile tmp_sid --indv f1 i1 s3 --make-just-psam --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# --iid-sid: two-token IDs are IID-SID instead of FID-IID, for --id-delim
# and --indv.
cat > tmp_delim.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	i1_s1	i1_s2	i2_s1
1	100	v1	A	G	.	.	.	GT	0/1	0/0	1/1
VCF
check_samples "i1:s1: i1:s2: i2:s1:" --vcf tmp_delim.vcf --id-delim _
check_samples ":i1:s1 :i1:s2 :i2:s1" --vcf tmp_delim.vcf --id-delim _ --iid-sid
check_samples ":i1:s2" --vcf tmp_delim.vcf --id-delim _ --iid-sid --indv i1 s2
check_samples "i1:s2:" --vcf tmp_delim.vcf --id-delim _ --indv i1 s2
check_samples ":i1:s1 :i1:s2" --vcf tmp_delim.vcf --id-delim _ --iid-sid --indv i1

# --strict-sid0: the dataset has no SID column, but the --keep/--remove/
# --update-sex files and --indv do.  By default, their SIDs are ignored.
printf '#FID\tIID\tPAT\tMAT\tSEX\nf1\ti1\t0\t0\t1\nf1\ti2\t0\t0\t1\nf3\tp3\t0\t0\t1\nf4\ti4\t0\t0\t2\n' > tmp_nosid.psam
printf 'a\nb\nd\nh\n' > tmp_nosid_keep.txt
$plink2 --vcf tmp_data.vcf --chr 1 --keep tmp_nosid_keep.txt --make-pgen --out tmp_nosid_g > /dev/null
cp tmp_nosid_g.pgen tmp_nosid.pgen
cp tmp_nosid_g.pvar tmp_nosid.pvar
printf '#FID\tIID\tSID\nf1\ti1\t0\nf1\ti2\tz\nf3\tp3\tz\nf4\ti4\t0\n' > tmp_keep_sid.txt
check_samples "f1:i1: f1:i2: f3:p3: f4:i4:" --pfile tmp_nosid --keep tmp_keep_sid.txt
check_samples "f1:i1: f4:i4:" --pfile tmp_nosid --keep tmp_keep_sid.txt --strict-sid0
check_samples "f1:i2: f3:p3:" --pfile tmp_nosid --remove tmp_keep_sid.txt --strict-sid0
check_samples "f1:i2:" --pfile tmp_nosid --indv f1 i2 z
check_samples "f1:i2:" --pfile tmp_nosid --indv f1 i2 0 --strict-sid0
if $plink2 --pfile tmp_nosid --indv f1 i2 z --strict-sid0 --make-just-psam --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
sex_col() {
    awk -F '\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "SEX") c = i; next} {printf "%s ", $c}' "$1"
}
printf '#FID\tIID\tSID\tSEX\nf1\ti1\t0\t2\nf1\ti2\tz\t2\n' > tmp_update_sex.txt
$plink2 --pfile tmp_nosid --update-sex tmp_update_sex.txt --make-just-psam --out tmp_us > /dev/null
test "$(sex_col tmp_us.psam)" = "2 2 1 2 "
$plink2 --pfile tmp_nosid --update-sex tmp_update_sex.txt --strict-sid0 --make-just-psam --out tmp_us0 > /dev/null
test "$(sex_col tmp_us0.psam)" = "2 1 1 2 "

# --- Phenotype loading ---
cat > tmp_ph.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	a	b	c	d	e
1	100	v1	A	G	.	.	.	GT	0/1	0/0	1/1	0/1	0/0
VCF
$plink2 --vcf tmp_ph.vcf --make-pgen --out tmp_ph_g > /dev/null
for prefix in tmp_ph tmp_fam_text tmp_fam_num tmp_plain; do
    cp tmp_ph_g.pgen $prefix.pgen
    cp tmp_ph_g.pvar $prefix.pvar
done
# P1 has -9 and -8.5; P2 has -9 (binary if that's missing); P3 is
# categorical with a 'NONE' entry.
printf '#FID\tIID\tSEX\tP1\tP2\tP3\nfA\ta\t1\t1.5\t-9\tx\nfA\tb\t2\t-9\t2\ty\nfB\tc\t1\t2.5\t1\tNONE\nfB\td\t2\t-8.5\t2\tx\nfC\te\t1\t3\t1\tz\n' > tmp_ph.psam
cp tmp_ph.psam tmp_fam_text.psam
printf '#FID\tIID\tSEX\nfA\ta\t1\nfA\tb\t2\nfB\tc\t1\nfB\td\t2\nfC\te\t1\n' > tmp_plain.psam
printf '#FID\tIID\tSEX\n7\ta\t1\n7\tb\t2\n0\tc\t1\n8\td\t2\n9\te\t1\n' > tmp_fam_num.psam

# $1: expected values (columns 2 and 4+ of the output .psam, comma-separated
# within a line, space-separated between lines); remaining args: flags.
check_psam() {
    local expected=$1
    shift
    rm -f tmp_p.psam
    $plink2 "$@" --make-just-psam --out tmp_p > /dev/null
    test "$(cut -f 2,4- tmp_p.psam | tr '\t\n' ', ')" = "$expected "
}

# -9 near another value in [-10, -8] is an error unless
# --neg9-pheno-really-missing or --no-input-missing-phenotype is given.
if $plink2 --pfile tmp_ph --make-just-psam --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
check_psam "IID,P1,P2,P3 a,1.5,NA,x b,NA,2,y c,2.5,1,NONE d,-8.5,2,x e,3,1,z" --pfile tmp_ph --neg9-pheno-really-missing
check_psam "IID,P1,P2,P3 a,1.5,-9,x b,-9,2,y c,2.5,1,NONE d,-8.5,2,x e,3,1,z" --pfile tmp_ph --no-input-missing-phenotype
check_psam "IID,P1,P2,P3 a,1.5,NA,NA b,NA,2,NA c,2.5,1,NA d,-8.5,2,NA e,3,1,NA" --pfile tmp_ph --neg9-pheno-really-missing --no-categorical

# --missing-catname: 'x' instead of 'NONE' marks a missing category.  (The
# output uses the same string, so check via --require-pheno.)
check_psam "IID,P1,P2,P3 a,1.5,NA,x b,NA,2,y d,-8.5,2,x e,3,1,z" --pfile tmp_ph --neg9-pheno-really-missing --require-pheno P3
check_psam "IID,P1,P2,P3 b,NA,2,y c,2.5,1,NONE e,3,1,z" --pfile tmp_ph --neg9-pheno-really-missing --require-pheno P3 --missing-catname x

# --not-pheno and --not-covar take names (comma- or space-separated).
check_psam "IID,P2 a,NA b,2 c,1 d,2 e,1" --pfile tmp_ph --neg9-pheno-really-missing --not-pheno P1,P3
check_psam "IID,P1 a,1.5 b,NA c,2.5 d,-8.5 e,3" --pfile tmp_ph --neg9-pheno-really-missing --not-pheno P2 P3
printf 'FID IID Q1 Q2 Q3\nfA a 1 2 3\nfA b 4 5 6\nfB c 7 8 9\nfB d 1 1 1\nfC e 2 2 2\n' > tmp_pheno.txt
check_psam "IID,P1,P2,P3,Q2 a,1.5,NA,x,2 b,NA,2,y,5 c,2.5,1,NONE,8 d,-8.5,2,x,1 e,3,1,z,2" --pfile tmp_ph --neg9-pheno-really-missing --pheno tmp_pheno.txt --not-pheno Q1 Q3
$plink2 --pfile tmp_ph --neg9-pheno-really-missing --covar tmp_pheno.txt --not-covar Q2 --make-just-psam --write-covar --out tmp_nc > /dev/null
test "$(head -n 1 tmp_nc.cov)" = "$(printf '#FID\tIID\tQ1\tQ3')"
$plink2 --pfile tmp_ph --neg9-pheno-really-missing --covar tmp_pheno.txt --covar-name Q1-Q2 --not-covar Q2 --make-just-psam --write-covar --out tmp_nc2 > /dev/null
test "$(head -n 1 tmp_nc2.cov)" = "$(printf '#FID\tIID\tQ1')"

# --mpheno n is --pheno-col-nums n+2.
check_psam "IID,Q2 a,2 b,5 c,8 d,1 e,2" --pfile tmp_ph --pheno tmp_pheno.txt --mpheno 2
check_psam "IID,Q2 a,2 b,5 c,8 d,1 e,2" --pfile tmp_ph --pheno tmp_pheno.txt --pheno-col-nums 4

# --family: categorical phenotype from FID.  Numeric FIDs get a 'C' prefix.
check_psam "IID,P1,P2,P3,CATPHENO a,1.5,NA,x,fA b,NA,2,y,fA c,2.5,1,NONE,fB d,-8.5,2,x,fB e,3,1,z,fC" --pfile tmp_fam_text --neg9-pheno-really-missing --family
check_psam "IID,P1,P2,P3,FAMCAT a,1.5,NA,x,fA b,NA,2,y,fA c,2.5,1,NONE,NONE d,-8.5,2,x,NONE e,3,1,z,fC" --pfile tmp_fam_text --neg9-pheno-really-missing --family FAMCAT --family-missing-catname fB
check_psam "IID,CATPHENO a,C7 b,C7 c,C0 d,C8 e,C9" --pfile tmp_fam_num --family
check_psam "IID,CATPHENO a,C7 b,C7 c,NONE d,C8 e,C9" --pfile tmp_fam_num --family --family-missing-catname 0
$plink2 --pfile tmp_fam_text --neg9-pheno-really-missing --family --missing-catname fC --require-pheno CATPHENO --make-just-psam --out tmp_fmc > /dev/null
test "$(tail -n +2 tmp_fmc.psam | cut -f 2 | tr '\n' ' ')" = "a b c d "
# Mixed numeric and non-numeric FIDs.
sed 's/^fC\t/7\t/' tmp_plain.psam > tmp_mixed.psam
if $plink2 --pgen tmp_plain.pgen --pvar tmp_plain.pvar --psam tmp_mixed.psam --family --make-just-psam --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
if $plink2 --pfile tmp_fam_num --family-missing-catname 9 --make-just-psam --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# --within/--mwithin: categories from column n+2 (default n=1).
printf 'fA a c1 5\nfA b c2 6\nfB c c1 NA\nfB d c3 7\nfC e NA 8\n' > tmp_within.txt
check_psam "IID,CATPHENO a,c1 b,c2 c,c1 d,c3 e,NONE" --pfile tmp_plain --within tmp_within.txt
check_psam "IID,CATPHENO a,c1 b,c2 c,c1 d,c3 e,NONE" --pfile tmp_plain --within tmp_within.txt --mwithin 1
check_psam "IID,WCAT a,C5 b,C6 c,NONE d,C7 e,C8" --pfile tmp_plain --within tmp_within.txt WCAT --mwithin 2
