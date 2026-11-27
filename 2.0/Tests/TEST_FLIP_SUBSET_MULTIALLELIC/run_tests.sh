#!/bin/bash

# --flip-subset on multiallelic variants, without sample subsetting.
#
# Flipping can turn 0/1 hets into rare-alt hets (here, A/C -> T/G = 2/3).
# When no sample subsetting was needed, the rare-alt patch entries were
# rewritten in place in the load buffer, so the larger result overwrote the
# following variants' records.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# 200 samples.  v1 (flipped): every sample is an A/C het, except for two
# unflipped samples with rare-alt genotypes (so that the variant has small
# rare-alt patch sections in the load buffer).  v2 and v3 (not
# flipped; non-SNP alleles): mixed rare-alt genotypes.
awk 'BEGIN {
  OFS = "\t";
  print "##fileformat=VCFv4.2";
  print "##contig=<ID=1,length=1000>";
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
  hdr = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (i = 1; i <= 200; ++i) hdr = hdr "\tS" i;
  print hdr;
  split("0/2 1/3 2/3 0/0 3/3 0/1 1/2", g, " ");
  l1 = "1\t100\tv1\tA\tC,G,T\t.\t.\t.\tGT";
  l2 = "1\t200\tv2\tA\tCC,CG,CT\t.\t.\t.\tGT";
  l3 = "1\t300\tv3\tAA\tCC,CG,CT\t.\t.\t.\tGT";
  for (i = 1; i <= 200; ++i) {
    l1 = l1 "\t" ((i == 2)? "0/2" : ((i == 200)? "1/2" : "0/1"));
    l2 = l2 "\t" g[(i % 7) + 1];
    l3 = l3 "\t" g[((i * 3) % 7) + 1];
  }
  print l1; print l2; print l3;
}' > tmp_data.vcf
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
printf 'v1\nv2\nv3\n' > tmp_flip.txt
awk 'BEGIN {for (i = 1; i <= 200; i += 2) print "S" i}' > tmp_subset.txt
$plink2 --pfile tmp_data --flip tmp_flip.txt permissive --flip-subset tmp_subset.txt --make-pgen --out tmp_flipped > /dev/null
$plink2 --pfile tmp_flipped --validate > /dev/null
$plink2 --pfile tmp_flipped --export vcf --out tmp_flipped > /dev/null
# v1: flipped samples (odd-numbered) become 2/3; v2 and v3 are unchanged.
# (The VCF exporter writes alleles in ascending order.)
awk '!/^#/ && $3 == "v1" {for (i = 10; i <= NF; ++i) printf("%s%s", $i, (i == NF)? "\n" : " ")}' tmp_flipped.vcf > tmp_v1.txt
awk 'BEGIN {for (i = 1; i <= 200; ++i) printf("%s%s", (i % 2)? "2/3" : ((i == 2)? "0/2" : ((i == 200)? "1/2" : "0/1")), (i == 200)? "\n" : " ")}' > tmp_v1_expected.txt
diff -q tmp_v1.txt tmp_v1_expected.txt
diff -q <(awk '!/^#/ && $3 != "v1" {$8 = "."; print}' tmp_flipped.vcf) <(awk '!/^#/ && $3 != "v1" {$8 = "."; print}' tmp_data.vcf)
