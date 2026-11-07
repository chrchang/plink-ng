#!/bin/bash

# --set-all-var-ids followed by --recover-var-ids must restore the original
# IDs, including when variants have been filtered out.
#
# Bug fixed in #666: when every remaining variant on a chromosome at or
# after a recovery-file position had been filtered out, that position was
# compared against the first variant of the *next* chromosome (or past the
# end of the variant array, for the last chromosome).  A next-chromosome
# variant with the same position and alleles was then also matched to the
# filtered variant's record, and the run failed with a "conflicting
# matching-position-and-alleles records" error.
#
# In a position-sorted recovery file, the scan has already moved past a
# chromosome's trailing filtered variants by the time their records come
# up, so this needs the filtered variant's record to come first: here, the
# original .pvar is unsorted, and the dataset was sorted (with --sort-vars)
# before its IDs were overwritten.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# chr1's last variant after sorting (rs3, bp 300) shares its position and
# alleles with chr2's first (rs4); likewise chr2's last (rs6, bp 500) and
# chr3's first (rs7).  rs3 and rs6 come first in their chromosomes here.
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##contig=<ID=2,length=1000>
##contig=<ID=3,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	300	rs3	G	A	.	.	.	GT	0/0	0/1
1	100	rs1	A	G	.	.	.	GT	0/1	0/0
1	200	rs2	C	T	.	.	.	GT	1/1	0/1
2	500	rs6	A	C	.	.	.	GT	0/0	1/1
2	300	rs4	G	A	.	.	.	GT	0/1	1/1
2	400	rs5	T	C	.	.	.	GT	0/1	0/0
3	500	rs7	A	C	.	.	.	GT	1/1	0/0
3	600	rs8	G	T	.	.	.	GT	0/1	0/1
VCF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_orig > /dev/null
$plink2 --pfile tmp_orig --sort-vars --set-all-var-ids '@:#:$r:$a' --make-pgen --out tmp_renamed > /dev/null
[[ "$(awk '!/^#/ {printf "%s ", $3}' tmp_renamed.pvar)" == "1:100:A:G 1:200:C:T 1:300:G:A 2:300:G:A 2:400:T:C 2:500:A:C 3:500:A:C 3:600:G:T " ]] || exit 1

ids() {
    awk '!/^#/ {printf "%s ", $3}' "$1"
}

# 1. No filtering: every ID comes back.
$plink2 --pfile tmp_renamed --recover-var-ids tmp_orig.pvar --make-just-pvar --out tmp_all > /dev/null
[[ "$(ids tmp_all.pvar)" == "rs1 rs2 rs3 rs4 rs5 rs6 rs7 rs8 " ]] || exit 1

# 2. The last variant of each chromosome is excluded, by its overwritten ID.
#    (--exclude-snps is applied before --recover-var-ids.  Load-time filters
#    like --snps-only remove variants entirely, and --exclude is applied
#    afterward, so neither reaches the affected code.)
$plink2 --pfile tmp_renamed --exclude-snps 1:300:G:A,2:500:A:C,3:600:G:T --recover-var-ids tmp_orig.pvar --make-just-pvar --out tmp_filtered > /dev/null
[[ "$(ids tmp_filtered.pvar)" == "rs1 rs2 rs4 rs5 rs7 " ]] || exit 1
