#!/bin/bash

# "--export bcf vcf-dosage=HDS-force" on unphased dosages used to write each
# dosage's HDS pair 64 samples too far to the right, at twice the right value,
# and (on chrX) with a second value for males.  Check that a BCF round trip
# through HDS matches a VCF round trip.

set -exo pipefail

# 100 samples (more than 64), unphased GT:DS on chr1 and chrX; a Park-Miller
# generator keeps the file independent of the awk implementation.
awk 'BEGIN {
    OFS = "\t"
    print "##fileformat=VCFv4.3"
    print "##contig=<ID=1,length=1000000>"
    print "##contig=<ID=X,length=100000000>"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage\">"
    line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (s = 0; s < 100; s++) {
        line = line "\ts" s
    }
    print line
    x = 12345
    for (v = 0; v < 8; v++) {
        chr = (v < 4)? "1" : "X"
        line = chr "\t" (5000000 + v * 100) "\tv" v "\tA\tC\t.\t.\t.\tGT:DS"
        for (s = 0; s < 100; s++) {
            x = (x * 16807) % 2147483647
            ds = (x % 2001) / 1000
            gt = (ds < 0.5)? "0/0" : ((ds < 1.5)? "0/1" : "1/1")
            # haploid males on chrX: dosage on 0..1 scale
            if ((chr == "X") && (s % 2 == 0)) {
                ds = (x % 1001) / 1000
                gt = (ds < 0.5)? "0" : "1"
            }
            line = line "\t" gt ":" ds
        }
        print line
    }
}' > dos.vcf
awk 'BEGIN {
    print "#IID\tSEX"
    for (s = 0; s < 100; s++) {
        print "s" s "\t" ((s % 2 == 0)? 1 : 2)
    }
}' > dos.psam

$1/plink2 $2 $3 --vcf dos.vcf dosage=DS --psam dos.psam --make-pgen --out tmp_data
$1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-dosage=HDS-force --out tmp_v
$1/plink2 $2 $3 --vcf tmp_v.vcf dosage=HDS --psam dos.psam --export vcf vcf-dosage=DS-force --out tmp_want
$1/plink2 $2 $3 --pfile tmp_data --export bcf vcf-dosage=HDS-force --out tmp_b
$1/plink2 $2 $3 --bcf tmp_b.bcf dosage=HDS --psam dos.psam --export vcf vcf-dosage=DS-force --out tmp_got
diff <(grep -v '^##' tmp_want.vcf) <(grep -v '^##' tmp_got.vcf)
