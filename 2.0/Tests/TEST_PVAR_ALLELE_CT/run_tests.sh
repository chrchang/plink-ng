#!/bin/bash

# A .pvar that understates a multiallelic variant's allele count must be
# rejected, not decoded with the wrong allele-code width.
#
# plink2 never stores allele counts in the .pgen header, so the .pvar is their
# only source.  A record written with 5 alleles but read with 3 has its
# allele codes read at the wrong width; before this change --freq,
# --geno-counts, --make-pgen and --export vcf all exited 0 with wrong
# genotypes.  The record parsers now check that a multiallelic-hardcall
# record ends where its last track does.
#
# The load-time check for a multiallelic record that the .pvar calls
# biallelic also has to cover the other filesets plink2 reads: the
# --pgen-diff and --flip-scan-ref-pfile filesets, and --pmerge inputs.

set -exo pipefail

expect_fail() {
    # usage: expect_fail <expected stderr substring> <plink2 args...>
    local msg=$1
    shift
    if $1/plink2 "${@:2}" 2> tmp_err.txt; then
        echo "expected failure: ${*:2}"
        exit 1
    fi
    grep -q "$msg" tmp_err.txt
}

for sep in '/' '|'; do
    if [ "$sep" = "/" ]; then p=u; else p=ph; fi
    cat > tmp_$p.vcf <<VCF
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6
1	10	v1	A	C,G,T,AA	.	.	.	GT	0${sep}4	1${sep}2	0${sep}0	3${sep}4	0${sep}1	2${sep}2
1	20	v2	A	C	.	.	.	GT	0${sep}1	1${sep}1	0${sep}0	0${sep}1	0${sep}0	1${sep}1
VCF
    $1/plink2 $2 $3 --vcf tmp_$p.vcf --make-pgen --out tmp_$p
    # The same records, with v1 declared as having 3 alleles.
    awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C,G" } { print }' tmp_$p.pvar > tmp_${p}3.pvar
    cmds=("--make-pgen" "--export vcf")
    if [ "$p" = "u" ]; then
        # These don't parse the phase track, so for phased records the
        # misread can't be seen from the end of the record.
        cmds+=("--freq" "--geno-counts")
    fi
    for c in "${cmds[@]}"; do
        expect_fail "Failed to unpack" $1 $2 $3 --pgen tmp_$p.pgen --pvar tmp_${p}3.pvar --psam tmp_$p.psam $c --out plink2_${p}3
    done
done

# Second filesets: v1 declared biallelic.
awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C" } { print }' tmp_u.pvar > tmp_ubi.pvar
cp tmp_u.pgen tmp_ubi.pgen
cp tmp_u.psam tmp_ubi.psam
expect_fail "in the --pgen-diff fileset has multiallelic hardcalls" $1 $2 $3 --pfile tmp_u --pgen-diff tmp_ubi --out plink2_diff
# (--flip-scan-ref-pfile gets the same check, but is still gated as under
# development.)

# --pmerge only reads genotype records one at a time when it merges
# same-position records, which master reaches through duplicate records inside
# one fileset.
cat > tmp_dup.vcf <<VCF
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6
1	10	v1	A	C,G	.	.	.	GT	0/2	0/1	./.	./.	0/0	1/2
1	10	v1	A	C	.	.	.	GT	./.	./.	0/1	1/1	0/0	./.
VCF
$1/plink2 $2 $3 --vcf tmp_dup.vcf --make-pgen --out tmp_dup
awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C" } { print }' tmp_dup.pvar > tmp_dupbi.pvar
cat > tmp_other.vcf <<VCF
##fileformat=VCFv4.3
##contig=<ID=2,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6
2	10	w1	A	C	.	.	.	GT	0/0	0/1	1/1	0/0	0/1	0/0
VCF
$1/plink2 $2 $3 --vcf tmp_other.vcf --make-pgen --out tmp_other
# Sanity check: with the correct .pvar this merge runs.
$1/plink2 $2 $3 --pfile tmp_other --pmerge tmp_dup --merge-mode nm-match --make-pgen --out plink2_merge_ok
expect_fail "in a --pmerge\[-list\] fileset has multiallelic hardcalls" $1 $2 $3 --pfile tmp_other --pmerge tmp_dup.pgen tmp_dupbi.pvar tmp_dup.psam --merge-mode nm-match --make-pgen --out plink2_merge_bad
