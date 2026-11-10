#!/bin/bash

# Regression tests for PLINK 1.9 association, regression, LD and epistasis
# fixes that change numeric output.  Each check failed before its fix.
#
# Usage: ./run_tests.sh {plink2 build dir} {up to 2 plink2 args}
# Requires plink (1.9) to be in the system PATH.

set -uo pipefail

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}

# --- Association and regression ---

# Genotype x covariate interaction columns on non-X haploid chromosomes used
# diploid 0/1/2 coding, while the main genotype term used haploid 0/1 coding.
# Compare with an explicit model (G coded 0/1, covariate C, and G*C).  Males
# are hemizygous on chrY; females are missing there.
awk 'BEGIN {
  s = 3
  for (i = 0; i < 120; ++i) {
    sex = (i < 90)? 1 : 2
    s = (s * 69069 + 1) % 4294967296
    q = (int(s / 65536) % 1000) / 100.0
    s = (s * 69069 + 1) % 4294967296
    b = 1 + int(s / 65536) % 2
    gl = ""
    for (j = 0; j < 4; ++j) {
      s = (s * 69069 + 1) % 4294967296
      r = int(s / 65536) % 3
      g = (sex == 2)? "0 0" : ((r == 0)? "A A" : "C C")
      if (j == 0 && r == 0) q += 2
      gl = gl " " g
    }
    print "f" i " i" i " 0 0 " sex " " q gl > "tmp_hap.ped"
    print "f" i, "i" i, b > "tmp_hap_cc.pheno"
    s = (s * 69069 + 1) % 4294967296
    print "f" i, "i" i, (int(s / 65536) % 100) / 10.0 > "tmp_hap_c.txt"
  }
}'
printf 'Y y1 0 100\nY y2 0 200\nY y3 0 300\nY y4 0 400\n' > tmp_hap.map
{ echo "FID IID C1"; cat tmp_hap_c.txt; } > tmp_hap.cov
plink --file tmp_hap --make-bed --out tmp_hap > /dev/null
plink --bfile tmp_hap --recode A --out tmp_hap > /dev/null
hap_ok=0
for model in linear logistic; do
    if [[ $model == logistic ]]; then
        pheno="--pheno tmp_hap_cc.pheno"
        out=assoc.logistic
        extra=beta
    else
        pheno=""
        out=assoc.linear
        extra=""
    fi
    plink --bfile tmp_hap $pheno --$model interaction $extra --covar tmp_hap.cov --out tmp_hap_$model > /dev/null
    for j in 0 1 2 3; do
        awk -v c=$((j + 7)) 'NR == FNR { if (FNR > 1) cov[$2] = $3; next }
          FNR == 1 { print "FID IID G C GC"; next }
          $c == "NA" { print $1, $2, -9, -9, -9; next }
          { g = $c / 2; print $1, $2, g, cov[$2], g * cov[$2] }' tmp_hap.cov tmp_hap.raw > tmp_hap_ref.cov
        plink --bfile tmp_hap $pheno --$model no-snp $extra --covar tmp_hap_ref.cov --out tmp_hap_ref > /dev/null
        awk -v snp=y$((j + 1)) '$2 == snp { print $7, $8 }' tmp_hap_$model.$out > tmp_hap_a.txt
        awk 'NR > 1 { print $3, $4 }' tmp_hap_ref.$out > tmp_hap_b.txt
        if ! (test -s tmp_hap_a.txt && cmp -s tmp_hap_a.txt tmp_hap_b.txt); then
            hap_ok=1
        fi
    done
done
check "chrY genotype x covariate interaction uses haploid coding" $hap_ok

# A joint test (--tests) of parameters that only exist on chrX (SEX and
# ADDxSEX, with 'interaction') has no constraints on autosomes.  The
# autosomal results then came out as nan / p = 0 instead of matching a run
# without --tests.
awk 'BEGIN {
  s = 9
  for (i = 0; i < 200; ++i) {
    sex = 1 + (i % 2)
    s = (s * 69069 + 1) % 4294967296
    q = (int(s / 65536) % 1000) / 100.0
    gl = ""
    for (j = 0; j < 6; ++j) {
      s = (s * 69069 + 1) % 4294967296
      r = int(s / 65536) % 3
      if (j >= 3 && sex == 1) {
        g = (r == 0)? "A A" : "C C"
      } else {
        g = (r == 0)? "A A" : ((r == 1)? "A C" : "C C")
      }
      if (r == 0) q += 0.7 * (j + 1) / 3
      gl = gl " " g
    }
    print "f" i " i" i " 0 0 " sex " " q gl > "tmp_jt.ped"
    print "f" i, "i" i, (q > 6)? 2 : 1 > "tmp_jt_cc.pheno"
  }
}'
printf '1 a1 0 100\n1 a2 0 200\n1 a3 0 300\nX x1 0 100\nX x2 0 200\nX x3 0 300\n' > tmp_jt.map
jt_ok=0
for model in linear logistic; do
    pheno=""
    if [[ $model == logistic ]]; then
        pheno="--pheno tmp_jt_cc.pheno"
    fi
    plink --file tmp_jt $pheno --$model interaction --out tmp_jt1 > /dev/null
    plink --file tmp_jt $pheno --$model interaction --tests 2,3 --out tmp_jt2 > /dev/null
    a=$(awk '$1 == 1' tmp_jt1.assoc.$model)
    b=$(awk '$1 == 1' tmp_jt2.assoc.$model)
    if ! (test -n "$a" && test "$a" = "$b"); then
        jt_ok=1
    fi
    # The chrX joint test itself must produce a statistic.
    if ! awk '$1 == 23 && $5 == "USER_2DF" && $NF != "NA" { ++n } END { exit (n != 3) }' tmp_jt2.assoc.$model; then
        jt_ok=1
    fi
done
check "--tests of chrX-only parameters leaves autosomal results unchanged" $jt_ok

# --logistic no-snp intercept: the INTERCEPT row lacked its two trailing NA
# columns and newline, so the next record was appended to it.
plink --dummy 200 2 --seed 4 --make-bed --out tmp_ns > /dev/null
awk 'BEGIN { s = 8; print "FID IID C1 C2" }
  { s = (s * 69069 + 1) % 4294967296; a = (int(s / 65536) % 100) / 10
    s = (s * 69069 + 1) % 4294967296; print $1, $2, a, int(s / 65536) % 2 }' tmp_ns.fam > tmp_ns.cov
ns_ok=0
for mods in "intercept" "intercept beta --ci 0.95"; do
    plink --bfile tmp_ns --logistic no-snp $mods --covar tmp_ns.cov --out tmp_ns > /dev/null
    if ! awk 'NR == 1 { nf = NF } NF != nf { bad = 1 } $1 == "INTERCEPT" { ++n }
              END { exit (bad || n != 1 || NR != 4) }' tmp_ns.assoc.logistic; then
        ns_ok=1
    fi
    if [[ "$(tail -c 1 tmp_ns.assoc.logistic | od -An -c | tr -d ' ')" != '\n' ]]; then
        ns_ok=1
    fi
done
check "--logistic no-snp intercept rows are complete" $ns_ok

# Fisher-Freeman-Halton mid-p (--model fisher-midp, GENO row) missed tables
# tied with the observed one.  Reference values are from Fisher23LnP() in
# github.com/chrchang/stats: 1/220 = 0.004545 and 0.7381.
fisher_midp() {
    # args: case AA/AC/CC counts, then control AA/AC/CC counts
    awk -v c="$*" 'BEGIN {
      split(c, n, " "); g[1] = "A A"; g[2] = "A C"; g[3] = "C C"; i = 0
      for (k = 1; k <= 6; ++k) {
        for (j = 0; j < n[k]; ++j) {
          ++i; print "f" i, "i" i, 0, 0, 1, (k <= 3)? 2 : 1, g[(k - 1) % 3 + 1]
        }
      }
    }' > tmp_fisher.ped
    echo "1 v1 0 100" > tmp_fisher.map
    plink --file tmp_fisher --model fisher-midp --cell 0 --out tmp_fisher > /dev/null
    awk '$5 == "GENO" { print $NF }' tmp_fisher.model
}
p1=$(fisher_midp 3 0 0 0 3 6)
p2=$(fisher_midp 0 1 10 1 0 10)
awk -v p1="$p1" -v p2="$p2" 'BEGIN {
  d1 = p1 / 0.00454545 - 1; d2 = p2 / 0.738095 - 1
  exit !(p1 != "" && p2 != "" && d1 < 1e-3 && d1 > -1e-3 && d2 < 1e-3 && d2 > -1e-3)
}'
check "--model fisher-midp 2x3 mid-p counts tied tables" $?

# --- LD and epistasis ---

# chrX LD (EM phasing with a haploid second locus, and the --r2 dprime chrX
# offsets) disagreed with plink2 for nearly every pair involving chrX.  Check
# windowed, inter-chr (X-autosome pairs in both orders) and --ld-snp-list
# reports against plink2 --r2-phased.
awk -f make_ld.awk
plink --file tmp_ld --make-bed --out tmp_ld > /dev/null
awk 'NR % 7 == 3 { print $2 }' tmp_ld.bim > tmp_ld_snps.txt
win="--ld-window 99999 --ld-window-kb 100000 --ld-window-r2 0"
ld_ok=0
plink --bfile tmp_ld --r2 dprime $win --out tmp_ld19_w > /dev/null
$1/plink2 ${2:-} ${3:-} --bfile tmp_ld --r2-phased cols=+dprime $win --out tmp_ld2_w > /dev/null
awk -f compare_ld.awk tmp_ld19_w.ld tmp_ld2_w.vcor || ld_ok=1
plink --bfile tmp_ld --r2 inter-chr dprime --ld-window-r2 0 --out tmp_ld19_i > /dev/null
$1/plink2 ${2:-} ${3:-} --bfile tmp_ld --r2-phased inter-chr cols=+dprime --ld-window-r2 0 --out tmp_ld2_i > /dev/null
awk -f compare_ld.awk tmp_ld19_i.ld tmp_ld2_i.vcor || ld_ok=1
plink --bfile tmp_ld --r2 dprime $win --ld-snp-list tmp_ld_snps.txt --out tmp_ld19_s > /dev/null
$1/plink2 ${2:-} ${3:-} --bfile tmp_ld --r2-phased cols=+dprime $win --ld-snp-list tmp_ld_snps.txt --out tmp_ld2_s > /dev/null
awk -f compare_ld.awk tmp_ld19_s.ld tmp_ld2_s.vcor || ld_ok=1
check "--r2 dprime on chrX matches plink2" $ld_ok

# Heterozygous calls in males on chrX must be treated as missing.  --show-tags
# applied haploid_fix() to the wrong buffer, and --blocks built its haploid
# masks from the wrong sample set (founders, rather than founders with a
# phenotype), so both gave different results from a fileset where those calls
# were already set to missing.  The chrX data has 4 blocks of 5 variants, and
# about 15% of male calls are heterozygous.
# (plink2 can't serve as the reference here: its --show-tags refuses haploid
# chromosomes, including chrX, and its chrX --r2-unphased codes males 0/2 with
# half weight where PLINK 1.9 codes them 0/1, so the r^2 values differ.)
awk 'BEGIN {
  s = 21; split("00000 11111 01011", h, " ")
  for (i = 0; i < 150; ++i) {
    sex = 1 + (i % 2)
    line = "f" i " i" i " 0 0 " sex " " ((i % 7 == 2)? -9 : 1 + int(i / 2) % 2)
    for (b = 0; b < 4; ++b) {
      s = (s * 69069 + 1) % 4294967296; h1 = h[1 + int(s / 65536) % 3]
      s = (s * 69069 + 1) % 4294967296; h2 = h[1 + int(s / 65536) % 3]
      for (k = 1; k <= 5; ++k) {
        a = substr(h1, k, 1); c = substr(h2, k, 1)
        if (sex == 1) {
          s = (s * 69069 + 1) % 4294967296
          g = (int(s / 65536) % 100 < 15)? "A C" : ((a == "0")? "A A" : "C C")
        } else {
          g = ((a == "0")? "A" : "C") " " ((c == "0")? "A" : "C")
        }
        line = line " " g
      }
    }
    print line
  }
}' > tmp_hh.ped
awk 'BEGIN { for (j = 0; j < 20; ++j) print "X x" j " 0 " (1000 * (j + 1)) }' > tmp_hh.map
plink --file tmp_hh --make-bed --out tmp_hh > /dev/null 2>&1
plink --bfile tmp_hh --set-hh-missing --make-bed --out tmp_hh_fixed > /dev/null 2>&1
printf 'x2\nx7\nx12\nx17\n' > tmp_hh_snps.txt
tags_ok=0
for r2 in 0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9; do
    plink --bfile tmp_hh --show-tags tmp_hh_snps.txt --tag-r2 $r2 --list-all --out tmp_tags1 > /dev/null 2>&1
    plink --bfile tmp_hh_fixed --show-tags tmp_hh_snps.txt --tag-r2 $r2 --list-all --out tmp_tags2 > /dev/null 2>&1
    if ! (test -s tmp_tags1.tags.list && cmp -s tmp_tags1.tags.list tmp_tags2.tags.list); then
        tags_ok=1
    fi
done
check "--show-tags ignores heterozygous haploid calls" $tags_ok
plink --bfile tmp_hh --blocks --out tmp_blocks1 > /dev/null 2>&1
plink --bfile tmp_hh_fixed --blocks --out tmp_blocks2 > /dev/null 2>&1
test -s tmp_blocks1.blocks.det && cmp -s tmp_blocks1.blocks.det tmp_blocks2.blocks.det
check "--blocks ignores heterozygous haploid calls" $?
# With a different set of founders missing phenotypes, the same mask mismatch
# made --blocks loop forever, even with no heterozygous haploid calls.  perl's
# alarm() limits it to 60 seconds; the check is skipped without perl.
if command -v perl > /dev/null; then
    awk '{ i = substr($1, 2) + 0; $6 = (i % 7 == 3)? -9 : 1 + int(i / 2) % 2; print }' tmp_hh_fixed.fam > tmp_hh_m3.fam
    cp tmp_hh_fixed.bed tmp_hh_m3.bed
    cp tmp_hh_fixed.bim tmp_hh_m3.bim
    perl -e 'alarm 60; exec @ARGV' plink --bfile tmp_hh_m3 --blocks --out tmp_blocks3 > /dev/null 2>&1
    check "--blocks finishes when some founders lack phenotypes" $?
fi

# --r2 square bin wrote each row with the size of the last variant block
# rather than the full row, so when --memory forced several blocks the matrix
# was truncated.  The binary matrix must match the text one.
plink --dummy 100000 700 0.01 --seed 2 --make-bed --out tmp_sq > /dev/null
plink --bfile tmp_sq --r2 square bin --memory 64 --out tmp_sq_bin > /dev/null
plink --bfile tmp_sq --r2 square --memory 64 --out tmp_sq_txt > /dev/null
od -An -v -tf8 tmp_sq_bin.ld.bin | tr -s ' ' '\n' | grep . > tmp_sq_bin.txt
tr '\t' '\n' < tmp_sq_txt.ld | grep . > tmp_sq_txt.txt
test "$(wc -l < tmp_sq_bin.txt)" -eq 490000 && test "$(wc -l < tmp_sq_txt.txt)" -eq 490000 && \
paste tmp_sq_bin.txt tmp_sq_txt.txt | awk '
  { d = $1 - $2; if (d < 0) d = -d; a = ($2 < 0)? -$2 : $2; if (d > 1e-5 + 1e-5 * a) bad = 1 }
  END { exit bad }'
check "--r2 square bin matches the text matrix across blocks" $?
rm -f tmp_sq.bed tmp_sq_bin.txt tmp_sq_txt.txt tmp_sq_txt.ld tmp_sq_bin.ld.bin

exit $failed
