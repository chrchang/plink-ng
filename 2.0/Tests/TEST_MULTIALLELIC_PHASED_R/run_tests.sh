#!/bin/bash

# Phased r between a biallelic variant and the major allele of a multiallelic
# variant, when both are in perfect LD.
#
# The phaseinfo bit of an x/y het (0 < x < y) refers to y.  Get1MP() used to
# return it unchanged when counting x, so every 1|2 or 2|1 het had its phase
# flipped and --r-phased reported |r| well below 1 whenever ALT1 (or ALT2 of a
# variant with more alleles) was the major allele.
#
# Every haplotype has one of three types, and each variant maps a type to one
# allele, so all pairs below are in perfect LD.  chr1 is dense in ALTx/ALTy
# hets; chr2 has only a few, which its multiallelic records store as a sparse
# list (aux1b mode 1).  m3 covers 1-bit allele codes, m5 2-bit, m7 4-bit and
# m20 8-bit ones; m5, m7 and m20 have ALT2+ as the major allele and also have
# hets whose other allele is lower.

set -exo pipefail

python3 -c "
n = 400
# (type 0 share, type 1 share) per chromosome; type 0 is the major one
shares = {1: (0.5, 0.25), 2: (0.7, 0.01)}
# allele for types 0, 1, 2
variants = {1: [('b', 2, (1, 0, 0)), ('m3', 3, (1, 2, 0)), ('m5', 5, (2, 3, 1)),
                ('m7', 7, (2, 6, 4)), ('m20', 20, (5, 17, 3))],
            2: [('b2', 2, (1, 0, 0)), ('s3', 3, (1, 2, 0)), ('s5', 5, (1, 3, 0))]}
seed = 12345
def rand():
    global seed
    seed = (seed * 1103515245 + 12345) % 2147483648
    return seed / 2147483648.0
print('##fileformat=VCFv4.3')
print('##contig=<ID=1,length=1000>')
print('##contig=<ID=2,length=1000>')
print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join('s%d' % i for i in range(n)))
for chrom in (1, 2):
    s0, s1 = shares[chrom]
    types = []
    for i in range(2 * n):
        x = rand()
        types.append(0 if x < s0 else (1 if x < s0 + s1 else 2))
    for pos, (vid, allele_ct, amap) in enumerate(variants[chrom]):
        alts = ','.join('A' * (k + 1) + 'C' for k in range(allele_ct - 1))
        gts = ['%d|%d' % (amap[types[2 * i]], amap[types[2 * i + 1]]) for i in range(n)]
        print('%d\t%d\t%s\tG\t%s\t.\t.\t.\tGT\t' % (chrom, 10 * (pos + 1), vid, alts) + '\t'.join(gts))
" > tmp_data.vcf
$1/plink2 $2 $3 --vcf tmp_data.vcf --make-pgen --out tmp_data
python3 -c "print('\n'.join('s%d' % i for i in range(0, 400, 3)))" > tmp_keep.txt

for keep in "" "--keep tmp_keep.txt"; do
    $1/plink2 $2 $3 --pfile tmp_data $keep --r-phased --ld-window-r2 0 --out plink2_r
    # 10 pairs on chr1 and 3 on chr2, all with |r| = 1
    awk 'NR > 1 { ++ct; r = $9 + 0; if (r < 0) { r = -r; } if (r < 0.99999) { print "bad r: " $0; exit 1; } } END { if (ct != 13) { print "pair count " ct; exit 1; } }' plink2_r.vcor
done
