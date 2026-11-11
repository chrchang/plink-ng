#!/bin/bash

# Regression tests for PLINK 1.9 fixes, checked against PLINK 1.07 where it
# has an equivalent command.  Each check failed before its fix.
#
# Usage: ./run_tests.sh
# Requires plink (1.9) to be in the system PATH.  The PLINK 1.07 binary is
# taken from $PLINK107, defaulting to the build in the repository's 1.07/
# directory; the checks that need it are skipped if it isn't there.

set -uo pipefail

plink107="${PLINK107:-../../../1.07/plink}"
have107=0
if [[ -x "$plink107" ]]; then
    have107=1
else
    echo "note: PLINK 1.07 not found at $plink107; skipping checks that need it"
fi

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
skip() {
    echo "skipped: $1"
}
# 32-bit LCG, small enough to be exact in awk's doubles.
lcg='function rnd() { s = (s * 69069 + 1) % 4294967296; return s / 4294967296 }'

# --genome excludes non-autosomal variants from marker_ct, but only skipped
# haploid/MT chromosomes after crossing a chromosome boundary, so a leading
# chrY (or MT) chromosome was still read.  The fileset starts with 50 chrY
# variants, followed by chromosomes 1 and 2; PLINK 1.07 skips chrY.
# 10 families of 3: two full siblings and one unrelated sample.
awk "$lcg"'
BEGIN {
  s = 77
  for (j = 0; j < 400; ++j) p[j] = 0.2 + 0.6 * rnd()
  for (f = 0; f < 10; ++f) {
    for (k = 0; k < 4; ++k) for (j = 0; j < 400; ++j) par[k, j] = (rnd() < p[j])? "A" : "G"
    for (c = 0; c < 3; ++c) {
      sex = 1 + (c % 2)
      line = "fam" f " c" c " 0 0 " sex " -9"
      for (j = 0; j < 400; ++j) {
        if (c < 2) {
          a = par[int(rnd() * 2), j]
          b = par[2 + int(rnd() * 2), j]
        } else {
          a = (rnd() < p[j])? "A" : "G"
          b = (rnd() < p[j])? "A" : "G"
        }
        if (j < 50) {
          if (sex == 1) b = a
          else a = b = "0"
        }
        line = line " " a " " b
      }
      print line
    }
  }
}' > tmp_gn.ped
awk 'BEGIN { for (j = 0; j < 400; ++j) print ((j < 50)? "Y" : ((j < 230)? 1 : 2)) " v" j " 0 " (100000 * (j + 1)) }' > tmp_gn.map
# --file would sort chrY after the autosomes, so assemble the .bed directly.
plink --file tmp_gn --chr y --make-bed --out tmp_gn_y > /dev/null 2>&1
plink --file tmp_gn --autosome --make-bed --out tmp_gn_a > /dev/null 2>&1
cat tmp_gn_y.bim tmp_gn_a.bim > tmp_gn.bim
cp tmp_gn_a.fam tmp_gn.fam
{ cat tmp_gn_y.bed; tail -c +4 tmp_gn_a.bed; } > tmp_gn.bed
if [[ $have107 -eq 1 ]]; then
    plink --bfile tmp_gn --genome --out tmp_gn19 > /dev/null
    "$plink107" --noweb --bfile tmp_gn --genome --out tmp_gn107 > /dev/null
    test "$(head -n 1 tmp_gn.bim | cut -f1)" = "24" && \
        test -s tmp_gn19.genome && cmp -s tmp_gn19.genome tmp_gn107.genome
    check "--genome skips a leading chrY, as PLINK 1.07 does" $?
else
    skip "--genome skips a leading chrY, as PLINK 1.07 does"
fi

# --cluster --match: after a covariate value equal to the --match missing
# code, every later covariate was checked with the wrong --match-type entry.
# C1 is a negative match ('-': members must differ) with missing values, and
# C2 is a positive match ('+': members must agree).  PLINK 1.07 has no missing
# code, so it gets a unique C1 value per missing sample instead, which under a
# negative match imposes no constraint, just as a missing value does.
plink --dummy 24 400 0.01 --seed 9 --make-bed --out tmp_mt > /dev/null
awk '{ i = NR - 1; print $1, $2, (i % 5 == 2)? "NA" : "g" (i % 3), "h" (i % 2) }' tmp_mt.fam > tmp_mt19.match
awk '{ if ($3 == "NA") $3 = "u_" $2; print }' tmp_mt19.match > tmp_mt107.match
echo "- +" > tmp_mt.type
if [[ $have107 -eq 1 ]]; then
    plink --bfile tmp_mt --cluster --match tmp_mt19.match NA --match-type tmp_mt.type --out tmp_mt19 > /dev/null
    "$plink107" --noweb --bfile tmp_mt --cluster --match tmp_mt107.match --match-type tmp_mt.type --out tmp_mt107 > /dev/null
    mt_ok=0
    for f in 1 2 3; do
        if ! (test -s tmp_mt19.cluster$f && cmp -s tmp_mt19.cluster$f tmp_mt107.cluster$f); then
            mt_ok=1
        fi
    done
    check "--cluster --match with missing values matches PLINK 1.07" $mt_ok
else
    skip "--cluster --match with missing values matches PLINK 1.07"
fi

# --homozyg group-verbose printed the consensus alleles of reversed variants
# (those whose .bim A1 is the major allele) swapped.  12 of 40 samples share
# a long homozygous haplotype, and every third variant has A1 set to its
# major allele.  The two versions form pools differently, so
# compare_hom_verbose.awk compares the individual genotype cells and the
# consensus alleles rather than whole files.
awk "$lcg"'
BEGIN {
  s = 31
  for (j = 0; j < 300; ++j) { p[j] = 0.25 + 0.5 * rnd(); hap[j] = (rnd() < 0.5)? "A" : "G" }
  for (i = 0; i < 40; ++i) {
    line = "f" i " i" i " 0 0 " (1 + i % 2) " " (1 + (i < 20))
    for (j = 0; j < 300; ++j) {
      if (i < 12 && j >= 80 + (i % 4) * 5 && j < 220 - (i % 3) * 7) {
        g = hap[j] " " hap[j]
      } else {
        a = (rnd() < p[j])? "A" : "G"
        b = (rnd() < p[j])? "A" : "G"
        g = a " " b
      }
      line = line " " g
    }
    print line
  }
}' > tmp_hz.ped
awk 'BEGIN { for (j = 0; j < 300; ++j) print "1 v" j " 0 " (10000 * (j + 1)) }' > tmp_hz.map
plink --file tmp_hz --make-bed --out tmp_hz0 > /dev/null
awk 'NR % 3 == 0 { print $2, $6 }' tmp_hz0.bim > tmp_hz_a1.txt
cut -d' ' -f1 tmp_hz_a1.txt > tmp_hz_reversed.txt
plink --bfile tmp_hz0 --a1-allele tmp_hz_a1.txt 2 1 --keep-allele-order --make-bed --out tmp_hz > /dev/null
hz="--homozyg-snp 50 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-het 1 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05"
if [[ $have107 -eq 1 ]]; then
    rm -f tmp_hz19.hom.overlap.S*.verbose tmp_hz107.hom.overlap.S*.verbose
    plink --bfile tmp_hz --homozyg group-verbose $hz --out tmp_hz19 > /dev/null
    "$plink107" --noweb --bfile tmp_hz --homozyg-group --homozyg-verbose $hz --out tmp_hz107 > /dev/null
    n107=$(ls tmp_hz107.hom.overlap.S*.verbose 2> /dev/null | wc -l | tr -d ' ')
    test "$n107" -gt 0 && \
        awk -v n107="$n107" -v reversed=tmp_hz_reversed.txt -f compare_hom_verbose.awk \
            tmp_hz107.hom.overlap.S*.verbose tmp_hz19.hom.overlap.S*.verbose
    check "--homozyg group-verbose alleles match PLINK 1.07" $?
else
    skip "--homozyg group-verbose alleles match PLINK 1.07"
fi

# --dosage --score with --q-score-range: the averaged score of every range
# used the first range's missing-dosage counts as its denominator.  PLINK
# 1.07 can't check this: it never averages, and for a missing dosage it keeps
# the previous variant's value.  Instead, each range's .profile must match a
# plain --score run whose score file holds only that range's variants.
awk "$lcg"'
BEGIN {
  s = 5
  printf "SNP A1 A2" > "tmp_ds.dose"
  for (i = 0; i < 30; ++i) {
    printf " f%d i%d", i, i > "tmp_ds.dose"
    print "f" i, "i" i, 0, 0, 1 + i % 2, (i % 3 == 0)? 2 : 1 > "tmp_ds.fam"
  }
  print "" > "tmp_ds.dose"
  for (j = 0; j < 40; ++j) {
    line = "v" j " A G"
    for (i = 0; i < 30; ++i) {
      u = rnd()
      d = 2 * rnd()
      line = line " " ((u < 0.15)? "NA" : sprintf("%.3f", d))
    }
    print line > "tmp_ds.dose"
    print "1 v" j " 0 " (1000 * (j + 1)) > "tmp_ds.map"
    print "v" j, "A", sprintf("%.3f", 4 * rnd() - 2) > "tmp_ds.score"
    print "v" j, (j % 4) + 0.5 > "tmp_ds.q"
  }
}'
printf 'r1 0 1\nr2 1 2\nr3 2 4\n' > tmp_ds.ranges
ds="--dosage tmp_ds.dose format=1 --fam tmp_ds.fam --map tmp_ds.map"
plink $ds --score tmp_ds.score no-sum no-mean-imputation --q-score-range tmp_ds.ranges tmp_ds.q --out tmp_ds > /dev/null
ds_ok=0
while read -r r lo hi; do
    awk -v lo="$lo" -v hi="$hi" 'NR == FNR { if ($2 >= lo && $2 <= hi) keep[$1] = 1; next } $1 in keep' tmp_ds.q tmp_ds.score > tmp_ds_$r.score
    plink $ds --score tmp_ds_$r.score no-sum no-mean-imputation --out tmp_ds_$r > /dev/null
    if ! (test -s tmp_ds.$r.profile && cmp -s tmp_ds.$r.profile tmp_ds_$r.profile); then
        ds_ok=1
    fi
done < tmp_ds.ranges
check "--dosage --q-score-range averages use each range's own counts" $ds_ok

exit $failed
