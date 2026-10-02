#!/bin/bash

# --make-pgen multiallelics=+ (variant-join).
# 1. Round trip: multiallelics=- followed by multiallelics=+ must reproduce
#    the original variants, up to ALT allele order, with the same genotypes
#    and phase.
# 2. Small fixtures for the .pvar fields, conflicting calls, and the
#    join-mode rules.

set -exo pipefail

# Canonical form of each data line: CHROM POS REF, sorted ALT alleles, then
# one genotype per sample spelled with allele strings.  Homozygous and
# unphased calls are written as 'x/y' with x <= y; phased hets as 'x|y'.
canon_vcf() {
    grep -v '^#' $1 | awk -F'\t' '{
        n = split($5, alts, ",")
        alleles[0] = $4
        for (i = 1; i <= n; i++) {
            alleles[i] = alts[i]
        }
        # insertion sort of the ALT strings
        for (i = 2; i <= n; i++) {
            x = alts[i]
            for (j = i - 1; (j >= 1) && (alts[j] > x); j--) {
                alts[j + 1] = alts[j]
            }
            alts[j + 1] = x
        }
        line = $1 " " $2 " " $4 " " alts[1]
        for (i = 2; i <= n; i++) {
            line = line "," alts[i]
        }
        for (i = 10; i <= NF; i++) {
            gt = $i
            if (substr(gt, 1, 1) == ".") {
                line = line " ."
                continue
            }
            sep = (index(gt, "|"))? "|" : "/"
            split(gt, ab, /[|\/]/)
            a = alleles[ab[1]]
            b = alleles[ab[2]]
            if ((a == b) || (sep == "/")) {
                if (a > b) {
                    x = a
                    a = b
                    b = x
                }
                sep = "/"
            }
            line = line " " a sep b
        }
        print line
    }'
}

# Random phased multiallelic data.  Each site is all-SNP, all-deletion
# (multi-base REF), or all-insertion, so multiallelics=+ (i.e. +both) inverts
# multiallelics=-.  Some deletions share a position with a SNP; their REF
# alleles differ, so they must not be joined.
awk -v seed=1 'BEGIN {
    srand(seed)
    sample_ct = 50
    print "##fileformat=VCFv4.3"
    print "##contig=<ID=1,length=100000000>"
    print "##contig=<ID=2,length=100000000>"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (s = 1; s <= sample_ct; s++) {
        header = header "\ts" s
    }
    print header
    split("A C G T", bases, " ")
    pos = 0
    vidx = 0
    for (site = 1; site <= 300; site++) {
        chr = (site <= 150)? 1 : 2
        if (site == 151) {
            pos = 0
        }
        pos += 1 + int(rand() * 40)
        kind = int(rand() * 3)
        rec_ct = 1
        if ((kind == 0) && (rand() < 0.2)) {
            rec_ct = 2
        }
        for (rec = 0; rec < rec_ct; rec++) {
            if (rec) {
                kind = 1
            }
            delete alleles
            delete seen
            if (kind == 0) {
                ref = bases[1 + int(rand() * 4)]
            } else if (kind == 1) {
                ref = ""
                for (k = 0; k < 3; k++) {
                    ref = ref bases[1 + int(rand() * 4)]
                }
            } else {
                ref = bases[1 + int(rand() * 4)]
            }
            alleles[0] = ref
            seen[ref] = 1
            alt_ct = 1 + int(rand() * ((kind == 0)? 3 : 4))
            for (a = 1; a <= alt_ct; a++) {
                do {
                    if (kind == 0) {
                        alt = bases[1 + int(rand() * 4)]
                    } else if (kind == 1) {
                        alt = substr(ref, 1, 1 + int(rand() * 2))
                        extra = int(rand() * 3)
                        for (k = 0; k < extra; k++) {
                            alt = alt bases[1 + int(rand() * 4)]
                        }
                    } else {
                        alt = ref
                        extra = 1 + int(rand() * 3)
                        for (k = 0; k < extra; k++) {
                            alt = alt bases[1 + int(rand() * 4)]
                        }
                    }
                } while (alt in seen)
                seen[alt] = 1
                alleles[a] = alt
            }
            alt_str = alleles[1]
            for (a = 2; a <= alt_ct; a++) {
                alt_str = alt_str "," alleles[a]
            }
            ++vidx
            line = chr "\t" pos "\tv" vidx "\t" ref "\t" alt_str "\t.\t.\t.\tGT"
            ref_prob = rand()
            phase_prob = rand()
            for (s = 1; s <= sample_ct; s++) {
                if (rand() < 0.05) {
                    line = line "\t./."
                    continue
                }
                h1 = (rand() < ref_prob)? 0 : (1 + int(rand() * alt_ct))
                h2 = (rand() < ref_prob)? 0 : (1 + int(rand() * alt_ct))
                line = line "\t" h1 ((rand() < phase_prob)? "|" : "/") h2
            }
            print line
        }
    }
}' > rand.vcf

$1/plink2 $2 $3 --vcf rand.vcf --make-pgen --out tmp_orig
$1/plink2 $2 $3 --pfile tmp_orig --make-pgen multiallelics=- --out tmp_split
for mode in + +any; do
    $1/plink2 $2 $3 --pfile tmp_split --make-pgen multiallelics=$mode --out tmp_join
    $1/plink2 $2 $3 --pfile tmp_orig --export vcf --out tmp_orig
    $1/plink2 $2 $3 --pfile tmp_join --export vcf --out tmp_join
    canon_vcf tmp_orig.vcf > tmp_orig.txt
    canon_vcf tmp_join.vcf > tmp_join.txt
    diff tmp_orig.txt tmp_join.txt
done


# Same round trip with a sample subset and a new sample order.
awk 'NR > 1 && NR <= 31 {print $1}' tmp_orig.psam > tmp_keep.txt
awk 'NR > 1 {print $1}' tmp_orig.psam | sort -r > tmp_order.txt
$1/plink2 $2 $3 --pfile tmp_orig --keep tmp_keep.txt --indiv-sort f tmp_order.txt --make-pgen --out tmp_orig2
$1/plink2 $2 $3 --pfile tmp_split --keep tmp_keep.txt --indiv-sort f tmp_order.txt --make-pgen multiallelics=+ --out tmp_join2
$1/plink2 $2 $3 --pfile tmp_orig2 --export vcf --out tmp_orig2
$1/plink2 $2 $3 --pfile tmp_join2 --export vcf --out tmp_join2
canon_vcf tmp_orig2.vcf > tmp_orig2.txt
canon_vcf tmp_join2.vcf > tmp_join2.txt
diff tmp_orig2.txt tmp_join2.txt

# erase-phase
$1/plink2 $2 $3 --pfile tmp_orig --make-pgen erase-phase --out tmp_orig3
$1/plink2 $2 $3 --pfile tmp_split --make-pgen multiallelics=+ erase-phase --out tmp_join3
$1/plink2 $2 $3 --pfile tmp_orig3 --export vcf --out tmp_orig3
$1/plink2 $2 $3 --pfile tmp_join3 --export vcf --out tmp_join3
canon_vcf tmp_orig3.vcf > tmp_orig3.txt
canon_vcf tmp_join3.vcf > tmp_join3.txt
diff tmp_orig3.txt tmp_join3.txt
if grep -q '|' tmp_join3.txt; then
    exit 1
fi

# .pvar fields, conflicting calls, and join modes.
# At position 100, C (frequency 5/10) comes before G (5/12).  s1 has two
# copies of C and one of G, which is a conflict; s3's phased hets put G on
# the first haplotype; s4's put both ALT alleles on the second haplotype, so
# its call loses phase.  At position 300, all frequencies are zero, so ALT
# alleles are natural-sorted.
cat > fields.vcf << 'EOF2'
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##FILTER=<ID=q10,Description="x">
##FILTER=<ID=q20,Description="x">
##INFO=<ID=AC,Number=A,Type=Integer,Description="x">
##INFO=<ID=AF,Number=A,Type=Float,Description="x">
##INFO=<ID=DP,Number=1,Type=Integer,Description="x">
##INFO=<ID=RV,Number=R,Type=Integer,Description="x">
##INFO=<ID=FL,Number=0,Type=Flag,Description="x">
##FORMAT=<ID=GT,Number=1,Type=String,Description="x">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6
1	100	a	A	C	30	q10	AC=1;DP=10;AF=0.1;FL;RV=5,1	GT	1/1	0/1	0|1	0|1	./.	0/0
1	100	b	A	G	20	q20;q10	AC=2;DP=10;AF=0.2;RV=5,2	GT	0/1	0/1	1|0	0|1	0/1	0/0
1	100	c	A	AT	.	PASS	DP=7	GT	0/0	0/0	0/0	0/0	0/0	0/1
1	200	d	A	C	.	.	.	GT	0/1	0/0	0/0	0/0	0/0	0/0
1	200	d	A	G	.	.	.	GT	0/0	0/1	0/0	0/0	0/0	0/0
1	300	e1	A	G	.	.	.	GT	0/0	0/0	0/0	0/0	0/0	0/0
1	300	e2	A	AT	.	.	.	GT	0/0	0/0	0/0	0/0	0/0	0/0
1	300	e3	A	C	.	.	.	GT	0/0	0/0	0/0	0/0	0/0	0/0
1	300	e4	A	AC	.	.	.	GT	0/0	0/0	0/0	0/0	0/0	0/0
EOF2

cat > expected_both.txt << 'EOF2'
1	100	a;b	A	C,G	20	q10;q20	AC=1,2;DP=10;AF=0.1,0.2;FL;RV=5,1,2	. 1/2 2|1 1/2 . 0/0
1	100	c	A	AT	.	PASS	DP=7	0/0 0/0 0/0 0/0 0/0 0/1
1	200	d	A	C,G	.	.	.	0/1 0/2 0/0 0/0 0/0 0/0
1	300	e3;e1	A	C,G	.	.	.	0/0 0/0 0/0 0/0 0/0 0/0
1	300	e4;e2	A	AC,AT	.	.	.	0/0 0/0 0/0 0/0 0/0 0/0
EOF2

cat > expected_snps.txt << 'EOF2'
1	100	a;b	A	C,G	20	q10;q20	AC=1,2;DP=10;AF=0.1,0.2;FL;RV=5,1,2	. 1/2 2|1 1/2 . 0/0
1	100	c	A	AT	.	PASS	DP=7	0/0 0/0 0/0 0/0 0/0 0/1
1	200	d	A	C,G	.	.	.	0/1 0/2 0/0 0/0 0/0 0/0
1	300	e3;e1	A	C,G	.	.	.	0/0 0/0 0/0 0/0 0/0 0/0
1	300	e2	A	AT	.	.	.	0/0 0/0 0/0 0/0 0/0 0/0
1	300	e4	A	AC	.	.	.	0/0 0/0 0/0 0/0 0/0 0/0
EOF2

cat > expected_any.txt << 'EOF2'
1	100	a;b;c	A	C,G,AT	20	q10;q20	AC=1,2,.;DP=.;AF=0.1,0.2,.;FL;RV=5,1,2,.	. 1/2 2|1 1/2 . 0/3
1	200	d	A	C,G	.	.	.	0/1 0/2 0/0 0/0 0/0 0/0
1	300	e4;e2;e3;e1	A	AC,AT,C,G	.	.	.	0/0 0/0 0/0 0/0 0/0 0/0
EOF2

$1/plink2 $2 $3 --vcf fields.vcf --make-pgen --out tmp_fields
for mode in both snps any; do
    $1/plink2 $2 $3 --pfile tmp_fields --make-pgen multiallelics=+$mode varid-join --out tmp_fields_$mode
    grep -q 'Warning: 1 genotype call was set to missing' tmp_fields_$mode.log
    if [ $mode = any ]; then
        grep -q "Warning: 1 INFO value was set to '.'" tmp_fields_$mode.log
    fi
    $1/plink2 $2 $3 --pfile tmp_fields_$mode --export vcf --out tmp_fields_$mode
    # .pvar columns, then genotypes with '.|.'/'./.' collapsed to '.'
    grep -v '^#' tmp_fields_$mode.pvar > tmp_pvar.txt
    grep -v '^#' tmp_fields_$mode.vcf | cut -f 10- | sed 's/\.[|/]\./\./g; s/0|0/0\/0/g' | tr '\t' ' ' > tmp_gt.txt
    paste tmp_pvar.txt tmp_gt.txt > tmp_got.txt
    diff expected_$mode.txt tmp_got.txt
done

# Error cases.
# Duplicate ALT allele at the same position and REF.
sed 's/^1	100	b	A	G	/1	100	b	A	C	/' fields.vcf > dup.vcf
$1/plink2 $2 $3 --vcf dup.vcf --make-pgen --out tmp_dup
if $1/plink2 $2 $3 --pfile tmp_dup --make-pgen multiallelics=+ --out tmp_dup_join; then
    exit 1
fi
grep -q "Error: Variants 'a' and 'b' have the same position, REF allele, and ALT" tmp_dup_join.log

# Unsorted input.
(grep '^#' tmp_fields.pvar; grep -v '^#' tmp_fields.pvar | sort -k2,2nr) > tmp_unsorted.pvar
if $1/plink2 $2 $3 --pgen tmp_fields.pgen --pvar tmp_unsorted.pvar --psam tmp_fields.psam --make-pgen multiallelics=+ --out tmp_unsorted_join; then
    exit 1
fi
grep -q 'Error: Variant-join requires a sorted .pvar.' tmp_unsorted_join.log

# Dosages need 'erase-dosage'.
cat > dosage.vcf << 'EOF2'
##fileformat=VCFv4.3
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="x">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="x">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
1	100	a	A	C	.	.	.	GT:DS	0/1:0.9	0/0:0.1
1	100	b	A	G	.	.	.	GT:DS	0/0:0	0/1:1.2
EOF2
$1/plink2 $2 $3 --vcf dosage.vcf dosage=DS --make-pgen --out tmp_dosage
if $1/plink2 $2 $3 --pfile tmp_dosage --make-pgen multiallelics=+ --out tmp_dosage_join; then
    exit 1
fi
grep -q 'Error: Variant-join does not support dosages yet.' tmp_dosage_join.log
$1/plink2 $2 $3 --pfile tmp_dosage --make-pgen multiallelics=+ erase-dosage --out tmp_dosage_join
