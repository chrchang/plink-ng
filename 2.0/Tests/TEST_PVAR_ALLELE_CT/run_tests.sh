#!/bin/bash

# A .pvar that understates a multiallelic variant's allele count must be
# rejected, not decoded with the wrong allele-code width.
#
# (Exception: we can tolerate this when the .pvar entry is biallelic, and phase
# information is never requested from the .pgen .  See later comments on PR
# #539.)
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
    # (sample subset for the subsetting code paths)
    printf 's1\ns2\ns4\ns6\n' > tmp_keep.txt
    # The same records, with v1 declared as having 3 alleles.
    awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C,G" } { print }' tmp_$p.pvar > tmp_${p}3.pvar
    if [ "$sep" = "/"]; then
        for c in "--make-pgen" "--export vcf" "--freq" "--geno-counts" "--keep tmp_keep.txt --freq" "--keep tmp_keep.txt --geno-counts"; do
            expect_fail ".pvar entry for (0-based) variant" $1 $2 $3 --pgen tmp_$p.pgen --pvar tmp_${p}3.pvar --psam tmp_$p.psam $c --out plink2_${p}3
        done
    else
        # Check is opportunistic, not exhaustive.  --geno-counts does not try
        # to scan to the end of the multiallelic .pgen record, so it doesn't
        # trigger the error.
        for c in "--make-pgen" "--export vcf" "--freq" "--keep tmp_keep.txt --freq"; do
            expect_fail ".pvar entry for (0-based) variant" $1 $2 $3 --pgen tmp_$p.pgen --pvar tmp_${p}3.pvar --psam tmp_$p.psam $c --out plink2_${p}3
        done
    fi
done

# --flip-scan-ref-pfile gets the same check.  It is still gated as under
# development, so this only runs once the gate is lifted.  (--flip-scan wants
# at least 50 samples.)
python3 -c "
import random
r = random.Random(3)
n = 60
print('##fileformat=VCFv4.3')
print('##contig=<ID=1,length=100000>')
print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join('s%d' % i for i in range(n)))
print('1\t10\tv1\tA\tC,G,T,AA\t.\t.\t.\tGT\t' + '\t'.join('%d/%d' % (r.randrange(5), r.randrange(5)) for i in range(n)))
for v in range(2, 8):
    print('1\t%d\tv%d\tA\tC\t.\t.\t.\tGT\t' % (v * 10, v) + '\t'.join('%d/%d' % (r.randrange(2), r.randrange(2)) for i in range(n)))
" > tmp_flip.vcf
$1/plink2 $2 $3 --vcf tmp_flip.vcf --make-pgen --out tmp_flipref
$1/plink2 $2 $3 --pfile tmp_flipref --max-alleles 2 --make-pgen --out tmp_flipmain
awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C" } { print }' tmp_flipref.pvar > tmp_flipbi.pvar
if $1/plink2 $2 $3 --pfile tmp_flipmain --flip-scan --flip-scan-ref-pfile tmp_flipref --out plink2_flip_ok 2> tmp_err.txt; then
    expect_fail "Variant #1 in the --flip-scan reference fileset" $1 $2 $3 --pfile tmp_flipmain --flip-scan --flip-scan-ref-pfile tmp_flipref.pgen tmp_flipbi.pvar tmp_flipref.psam --out plink2_flip_bad
else
    grep -q "under development" tmp_err.txt
fi

# --pmerge only reads genotype records one at a time when it merges
# same-position records, which master reaches through duplicate records inside
# one fileset.
cat > tmp_dup.vcf <<VCF
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6
1	10	v1	A	C,G,T,AC	.	.	.	GT	0/2	0/4	./.	./.	0/0	1/2
1	10	v1	A	C	.	.	.	GT	./.	./.	0/1	1/1	0/0	./.
VCF
$1/plink2 $2 $3 --vcf tmp_dup.vcf --make-pgen --out tmp_dup
awk 'BEGIN { FS = OFS = "\t" } $3 == "v1" { $5 = "C,G" } { print }' tmp_dup.pvar > tmp_duptri.pvar
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
expect_fail ".pvar entry for (0-based) variant" $1 $2 $3 --pfile tmp_other --pmerge tmp_dup.pgen tmp_duptri.pvar tmp_dup.psam --merge-mode nm-match --make-pgen --out plink2_merge_bad
