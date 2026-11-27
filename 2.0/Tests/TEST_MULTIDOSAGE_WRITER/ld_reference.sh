#!/bin/bash
# Reference checks for --r2-[un]phased and --ld on multiallelic dosages.  At
# a multiallelic variant, these use the major allele's dosage (and
# dosage-phase) vs. the rest, which is exactly what the split variant for
# that allele (--make-pgen multiallelics=-) stores.  So, for fileset $2
# (variants with at most 5 alleles) and output prefix $3, this writes:
# * ${3}_perm: a copy with REF set to each multiallelic variant's least
#   frequent allele, so that the major allele is never REF (there's no
#   "REF vs. rest" split variant), and ${3}_maj.txt: <ID> <major allele>.
# * ${3}_split: its split, with IDs '@:#:$r:$a', and ${3}_split_ext.txt:
#   each variant's major-allele split variant (any one, for biallelic
#   variants), in the same order as ${3}_perm.
set -eo pipefail
P=$1; F=$2; O=$3
$P --pfile $F --max-alleles 5 --freq --out ${O}_f0 > /dev/null
awk 'BEGIN {OFS = "\t"}
  FNR == NR {
    if (FNR > 1) {
      n = split($5, f, ",");
      rf = 1;
      for (j = 1; j <= n; ++j) rf -= f[j];
      best = 0;
      bf = rf;
      for (j = 1; j <= n; ++j) if (f[j] < bf) { bf = f[j]; best = j }
      least[$2] = best;
    }
    next;
  }
  !/^#/ {
    n = split($5, a, ",");
    if ((n >= 2) && least[$3]) print $3, a[least[$3]];
  }' ${O}_f0.afreq $F.pvar > ${O}_newref.txt
$P --pfile $F --max-alleles 5 --ref-allele force ${O}_newref.txt 2 1 --make-pgen --out ${O}_perm > /dev/null
$P --pfile ${O}_perm --freq --out ${O}_pf > /dev/null
awk 'BEGIN {OFS = "\t"} NR > 1 {
  n = split($5, f, ",");
  split($4, a, ",");
  rf = 1;
  for (j = 1; j <= n; ++j) rf -= f[j];
  best = $3;
  bf = rf;
  for (j = 1; j <= n; ++j) if (f[j] > bf) { bf = f[j]; best = a[j] }
  if ((n >= 2) && (best == $3)) { print "major allele of " $2 " is still REF" > "/dev/stderr"; exit 1 }
  print $2, best;
}' ${O}_pf.afreq > ${O}_maj.txt
$P --pfile ${O}_perm --make-pgen multiallelics=- varid-dup --out ${O}_split0 > /dev/null
$P --pfile ${O}_split0 --set-all-var-ids '@:#:$r:$a' --make-pgen --out ${O}_split > /dev/null
awk 'BEGIN {OFS = "\t"} FNR == NR {maj[$1] = $2; next} !/^#/ {n = split($5, a, ","); print $1 ":" $2 ":" $4 ":" ((n == 1)? $5 : maj[$3])}' ${O}_maj.txt ${O}_perm.pvar > ${O}_split_ext.txt
