#!/bin/bash

# --split-par when every chrX variant is in PAR1.
#
# The PAR2 boundary used to be located with a binary search over the empty
# range of remaining chrX variants, which read the next chromosome's first
# position (or, when chrX was last, one past the end) and could report a
# PAR2 start beyond chrX.  The first variant of the following chromosome was
# then labeled chrX.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n' > tmp_par1.vcf
printf 'X\t100000\tx1\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\n' >> tmp_par1.vcf
printf 'X\t100100\tx2\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\n' >> tmp_par1.vcf
printf 'Y\t1480\ty1\tA\tG\t.\t.\t.\tGT\t0\t1\n' >> tmp_par1.vcf
printf 'Y\t1510\ty2\tA\tG\t.\t.\t.\tGT\t0\t1\n' >> tmp_par1.vcf

$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_par1.vcf --split-par b37 --make-just-pvar --out tmp_par1
test "$(grep -v '^#' tmp_par1.pvar | cut -f1 | tr '\n' ' ')" = "PAR1 PAR1 Y Y "

# chrX last in the file.
head -4 tmp_par1.vcf > tmp_par1x.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_par1x.vcf --split-par b37 --make-just-pvar --out tmp_par1x
test "$(grep -v '^#' tmp_par1x.pvar | cut -f1 | tr '\n' ' ')" = "PAR1 PAR1 "
