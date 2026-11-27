#!/bin/bash
# Reference checks for --glm on multiallelic dosages, built from plink2's
# biallelic code paths.  Writes, for fileset $2 and output prefix $3:
# * ${3}_multi.*, ${3}_multi1.*, ${3}_multis.*: --glm 'omit-ref' results on
#   the multiallelic variants (linear via GlmLinearSubbatchThread() for Q/Q2
#   and GlmLinearThread() for Q alone, logistic via GlmLogisticThreadD() and
#   (with 'single-prec-cc') GlmLogisticThreadF(), multinomial).
# * ${3}_ref.txt: reference rows (see glm_compare.awk) for each ALT allele of
#   up to $4 3-allele and $4 4-allele variants whose ALT alleles all have
#   nonzero dosages in at least 10 samples: --glm on the split biallelic
#   variant, with the other ALT alleles' dosages as covariates (the same
#   column space).
# * ${3}_dflt.*, ${3}_permref.txt: default allele omission (often a non-REF
#   allele, with a different column layout), vs. 'omit-ref' on a copy with
#   REF set to each variant's omitted allele (the same model).
# $5: extra flags (e.g. "--xchr-model 1").
set -eo pipefail
P=$1; F=$2; O=$3; MAXV=$4; X=$5
# deterministic phenotypes (not awk rand(), which differs across awks); Q and
# Q2 have the same missingness, so they're analyzed as a batch
awk 'BEGIN {OFS = "\t"; print "#IID", "Q", "Q2", "B", "M"} NR > 1 {h = (NR * 7919 + 13) % 1009; print $1, h / 1009, ((NR * 337) % 101) / 101, ((NR * 104729) % 7 < 3)? 2 : 1, "L" ((NR * 31) % 3)}' $F.psam > ${O}_ph.txt
G="$X --glm allow-no-covars omit-ref multinomial cols=+a1freq,+gcountcc"
$P --pfile $F --pheno ${O}_ph.txt $G --out ${O}_multi > /dev/null
$P --pfile $F --pheno ${O}_ph.txt --pheno-name Q $G --out ${O}_multi1 > /dev/null
$P --pfile $F --pheno ${O}_ph.txt --pheno-name B $G single-prec-cc --out ${O}_multis > /dev/null
# Covariate dosages must be on the scale glm uses internally: 0..2 for
# everyone (so they're exported with CHROM relabeled to 1; only the 3- and
# 4-allele variants are needed), except
# --xchr-model 1 chrX males (0..1).  (Other haploid chromosomes are rescaled
# uniformly, which doesn't matter.)
cp $F.pgen ${O}_auto.pgen
cp $F.psam ${O}_auto.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "1"; print}' $F.pvar > ${O}_auto.pvar
$P --pfile ${O}_auto --max-alleles 4 --export vcf vcf-dosage=DS-force --out ${O}_ds > /dev/null
halve_males=$(awk '!/^#/ {print ($1 == "X")? 1 : 0; exit}' $F.pvar)
case "$X" in *"--xchr-model 1"*) ;; *) halve_males=0 ;; esac
awk 'NR == 1 {for (i = 1; i <= NF; ++i) if ($i == "SEX") c = i; next} {print $1, (c && ($c == 1))? 1 : 0}' $F.psam > ${O}_males.txt
$P --pfile $F --make-pgen multiallelics=- varid-dup --out ${O}_split0 > /dev/null
$P --pfile ${O}_split0 --set-all-var-ids '@:#:$r:$a' --make-pgen --out ${O}_split > /dev/null
awk -v m=$MAXV '!/^#/ {
  n = split($5, a, ",");
  if ((n != 2) && (n != 3)) next;
  for (j = 1; j <= n; ++j) nz[j] = 0;
  for (i = 10; i <= NF; ++i) {
    split($i, f, ":");
    if (f[2] == ".") continue;
    split(f[2], d, ",");
    for (j = 1; j <= n; ++j) nz[j] += (d[j] != 0);
  }
  ok = 1;
  for (j = 1; j <= n; ++j) ok = ok && (nz[j] >= 10);
  if (ok && (++c[n] <= m)) print $3;
}' ${O}_ds.vcf > ${O}_vars.txt
rm -f ${O}_ref.txt
for v in $(cat ${O}_vars.txt); do
  chrom=$(awk -v v=$v '!/^#/ && ($3 == v) {print $1}' $F.pvar)
  line=$(grep -v '^#' ${O}_ds.vcf | awk -v v=$v '$3 == v')
  pos=$(echo "$line" | cut -f2)
  ref=$(echo "$line" | cut -f4)
  alts=$(echo "$line" | cut -f5)
  nalt=$(echo "$alts" | awk -F, '{print NF}')
  for a in $(seq 1 $nalt); do
    alt=$(echo "$alts" | cut -d, -f$a)
    grep -v '^##' ${O}_ds.vcf | awk -v v=$v -v a=$a -v nalt=$nalt -v halve_males=$halve_males 'BEGIN {OFS = "\t"}
      FNR == NR {male[$1] = $2; next}
      /^#CHROM/ {for (i = 10; i <= NF; ++i) id[i] = $i; ns = NF; next}
      $3 == v {
        hdr = "#IID";
        for (j = 1; j <= nalt; ++j) if (j != a) hdr = hdr "\tD" j;
        print hdr;
        for (i = 10; i <= ns; ++i) {
          split($i, f, ":");
          row = id[i];
          if (f[2] == ".") {
            for (j = 1; j <= nalt; ++j) if (j != a) row = row "\tNA";
          } else {
            split(f[2], d, ",");
            for (j = 1; j <= nalt; ++j) if (j != a) row = row "\t" ((halve_males && male[id[i]])? (d[j] / 2) : d[j]);
          }
          print row;
        }
      }' ${O}_males.txt - > ${O}_cov.txt
    echo "$chrom:$pos:$ref:$alt" > ${O}_ext.txt
    # (constant covariates are excluded with a warning, which is fine)
    $P --pfile ${O}_split --extract ${O}_ext.txt --pheno ${O}_ph.txt --covar ${O}_cov.txt $G hide-covar --out ${O}_bi > /dev/null 2>&1
    for ph in Q.glm.linear Q2.glm.linear B.glm.logistic.hybrid M.glm.multinomial.hybrid; do
      awk -v v=$v -v ph=$ph 'BEGIN {OFS = "\t"} !/^#/ {print ph, v, $0}' ${O}_bi.$ph >> ${O}_ref.txt
    done
  done
done
GD="$X --glm allow-no-covars multinomial cols=+a1freq,+gcountcc,+omitted"
$P --pfile $F --pheno ${O}_ph.txt $GD --out ${O}_dflt > /dev/null
awk 'BEGIN {OFS = "\t"} NR == 1 {for (i = 1; i <= NF; ++i) if ($i == "OMITTED") c = i; next} !seen[$3]++ {print $3, $c}' ${O}_dflt.Q.glm.linear > ${O}_omitted.txt
$P --pfile $F --ref-allele force ${O}_omitted.txt 2 1 --make-pgen --out ${O}_perm > /dev/null
$P --pfile ${O}_perm --pheno ${O}_ph.txt $GD omit-ref --out ${O}_permglm > /dev/null
rm -f ${O}_permref.txt
for ph in Q.glm.linear Q2.glm.linear B.glm.logistic.hybrid M.glm.multinomial.hybrid; do
  awk -v ph=$ph 'BEGIN {OFS = "\t"} !/^#/ {print ph, $3, $0}' ${O}_permglm.$ph >> ${O}_permref.txt
done
