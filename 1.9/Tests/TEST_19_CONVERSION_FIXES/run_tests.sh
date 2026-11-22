#!/bin/bash

# Regression tests for PLINK 1.9 output-formatting, data conversion and merge
# fixes.  Each check failed before its fix (several crashed).
#
# Usage: ./run_tests.sh {plink 1.9 dir} {plink2 build dir} {up to 2 plink2 args}
# Requires python3 in the system PATH.  plink (1.9) is taken from {plink 1.9
# dir}, or from the system PATH if that is empty or omitted.  {plink2 build
# dir} defaults to ../../../2.0/build_dynamic.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}
plink2="${2:-../../../2.0/build_dynamic}/plink2 ${3:-} ${4:-}"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
# Runs plink and prints its exit code.  A crash gives 128 or more.
rc_of() {
    plink "$@" > /dev/null 2>&1
    echo $?
}
# 32-bit LCG, small enough to be exact in awk's doubles.
lcg='function rnd() { s = (s * 69069 + 1) % 4294967296; return s / 4294967296 }'
geno='function geno(a, b) { r = int(rnd() * 3); return (r == 0)? a " " a : ((r == 1)? a " " b : b " " b) }'

python3 make_bcf.py
python3 make_bgen.py

# --- dtoa_e ---

# dtoa_e() (used for --make-grm-gz and --ibc values) had a threshold with a
# dropped '9', so values in [9.999995e-4, 9.9999995e-4) were printed as
# 1.000000e-03.  With one variant and --read-freq, two samples with no copies
# of an allele of frequency p have relationship 2p / (1 - p); p is chosen to
# make that 9.999999e-4.
printf 'f1 i1 0 0 1 1 C C\nf2 i2 0 0 1 1 C C\nf3 i3 0 0 1 1 A C\n' > tmp_de.ped
echo "1 v1 0 100" > tmp_de.map
p=$(awk 'BEGIN { v = 9.999999e-4; printf "%.17g", v / (2 + v) }')
printf ' CHR  SNP   A1   A2          MAF  NCHROBS\n   1   v1    A    C  %s        6\n' "$p" > tmp_de.frq
plink --file tmp_de --read-freq tmp_de.frq --make-grm-gz --out tmp_de > /dev/null
expected=$(awk -v p="$p" 'BEGIN { printf "%.6e", 2 * p / (1 - p) }')
test "$expected" = "9.999999e-04" && \
    test "$(gzip -dc tmp_de.grm.gz | awk '$1 == 2 && $2 == 1 { print $4 }')" = "$expected"
check "--make-grm-gz prints values just below 1e-3 correctly" $?

# --- Merging ---

# --merge with a text fileset sized its buffer from the variant count without
# the negative-bp lines that merge_main() still indexes, and overran it.  The
# result must match merging the same data as a binary fileset.
awk "$lcg $geno"'
BEGIN {
  s = 3
  for (i = 0; i < 20; ++i) {
    line = "a" i " a" i " 0 0 1 1"
    for (j = 0; j < 50; ++j) line = line " " geno("A", "C")
    print line > "tmp_mga.ped"
  }
  for (j = 0; j < 50; ++j) print "1 v" j " 0 " (1000 * (j + 1)) > "tmp_mga.map"
  for (i = 0; i < 20; ++i) {
    line = "b" i " b" i " 0 0 2 1"
    for (j = 0; j < 400; ++j) line = line " " geno("A", "C")
    print line > "tmp_mgb.ped"
  }
  # 350 of the 400 variants have negative bp, so they are skipped.
  for (j = 0; j < 400; ++j) {
    if (j % 8 == 0) print "1 v" (j / 8) " 0 " (1000 * (j / 8 + 1)) > "tmp_mgb.map"
    else print "1 n" j " 0 " (-(j + 1)) > "tmp_mgb.map"
  }
}'
plink --file tmp_mga --make-bed --out tmp_mga > /dev/null
plink --file tmp_mgb --make-bed --out tmp_mgb > /dev/null 2>&1
mg_ok=0
test "$(rc_of --bfile tmp_mga --merge tmp_mgb.ped tmp_mgb.map --make-bed --out tmp_mg_text)" -eq 0 || mg_ok=1
plink --bfile tmp_mga --bmerge tmp_mgb --make-bed --out tmp_mg_bin > /dev/null
cmp -s tmp_mg_text.bed tmp_mg_bin.bed && cmp -s tmp_mg_text.bim tmp_mg_bin.bim || mg_ok=1
check "--merge of a text fileset with negative-bp variants" $mg_ok

# --merge-mode 7 with --merge-equal-pos looked up alleles by the merged-in
# variant's own index, which --merge-equal-pos leaves empty when its ID sorts
# after the matching ID in the first fileset.  The .diff report (nonmissing
# discordances) must match plink2 --pgen-diff on the same data with the IDs
# already matching.  compare_diff.awk converts both to allele letters.
awk "$lcg"'
BEGIN {
  s = 11; split("A C G T", base, " ")
  for (j = 0; j < 30; ++j) {
    k = 1 + int(rnd() * 4); a1[j] = base[k]; a2[j] = base[1 + k % 4]
    line = (1 + (j >= 15)) " rs" ((j * 7) % 30) " 0 " (1000 * (j + 1))
    print line > "tmp_mdA.map"; print line > "tmp_mdB.map"
  }
  for (i = 0; i < 15; ++i) {
    la = "s" i " s" i " 0 0 1 1"; lb = la
    for (j = 0; j < 30; ++j) {
      r = int(rnd() * 3)
      ga = (r == 0)? a1[j] " " a1[j] : ((r == 1)? a1[j] " " a2[j] : a2[j] " " a2[j])
      if (rnd() < 0.1) r = (r + 1) % 3
      gb = (r == 0)? a1[j] " " a1[j] : ((r == 1)? a1[j] " " a2[j] : a2[j] " " a2[j])
      la = la " " ga; lb = lb " " gb
    }
    print la > "tmp_mdA.ped"; print lb > "tmp_mdB.ped"
  }
}'
plink --file tmp_mdA --make-bed --out tmp_mdA > /dev/null
plink --file tmp_mdB --make-bed --out tmp_mdB > /dev/null
awk 'BEGIN { OFS = "\t" } { $2 = "z_" $2; print }' tmp_mdB.bim > tmp_mdBz.bim
cp tmp_mdB.bed tmp_mdBz.bed
cp tmp_mdB.fam tmp_mdBz.fam
md_ok=0
test "$(rc_of --bfile tmp_mdA --bmerge tmp_mdBz --merge-equal-pos --merge-mode 7 --out tmp_md19)" -eq 0 || md_ok=1
$plink2 --bfile tmp_mdA --pgen-diff tmp_mdB.bed tmp_mdB.bim tmp_mdB.fam --out tmp_md2 > /dev/null
awk -v strip=z_ -f compare_diff.awk tmp_md19.diff 2> /dev/null | sort > tmp_md19.txt
awk -f compare_diff.awk tmp_mdB.bim tmp_md2.pdiff | sort > tmp_md2.txt
test -s tmp_md2.txt && cmp -s tmp_md19.txt tmp_md2.txt || md_ok=1
check "--merge-mode 7 --merge-equal-pos .diff matches plink2 --pgen-diff" $md_ok

# --- Text and Oxford imports ---

# --lfile with a .map that is unsorted and has variants excluded at load time
# (negative bp) filtered the variant indices twice.  It must match --file on
# the same genotypes, and plink2 --lfile (which keeps the .map's variant
# order, so that comparison is per call).
awk "$lcg"'
BEGIN {
  s = 13
  for (j = 0; j < 20; ++j) {
    pos = 1000 * (1 + (j * 7) % 20)
    if (j == 3 || j == 11) pos = -pos
    print "1 v" j " 0 " pos > "tmp_lg.map"
  }
  for (i = 0; i < 12; ++i) {
    line = "f" i " i" i " 0 0 1 1"
    print "f" i, "i" i, 0, 0, 1, 1 > "tmp_lg.fam"
    for (j = 0; j < 20; ++j) {
      r = int(rnd() * 4)
      g = (r == 0)? "A A" : ((r == 1)? "A C" : ((r == 2)? "C C" : "0 0"))
      line = line " " g
      split(g, a, " "); print "f" i, "i" i, "v" j, a[1], a[2] > "tmp_lg.lgen"
    }
    print line > "tmp_lg.ped"
  }
}'
lg_ok=0
test "$(rc_of --lfile tmp_lg --make-bed --out tmp_lg_l)" -eq 0 || lg_ok=1
plink --file tmp_lg --make-bed --out tmp_lg_p > /dev/null
cmp -s tmp_lg_l.bed tmp_lg_p.bed && cmp -s tmp_lg_l.bim tmp_lg_p.bim || lg_ok=1
plink --bfile tmp_lg_l --recode --out tmp_lg_l > /dev/null 2>&1
$plink2 --lfile tmp_lg --export ped --out tmp_lg_2 > /dev/null
awk -f ped_long.awk tmp_lg_l.map tmp_lg_l.ped 2> /dev/null | sort > tmp_lg_l.txt
awk -f ped_long.awk tmp_lg_2.map tmp_lg_2.ped | sort > tmp_lg_2.txt
test -s tmp_lg_2.txt && cmp -s tmp_lg_l.txt tmp_lg_2.txt || lg_ok=1
check "--lfile with an unsorted .map matches --file and plink2" $lg_ok

# A .bgen v1.1 SNP block that decompressed to fewer than 6 * sample count
# bytes was accepted, leaving probabilities uninitialized.  A valid file must
# still match --gen on the same probabilities.
bg_ok=0
plink --bgen ok.bgen --sample ok.sample --make-bed --out tmp_bg_ok > /dev/null
plink --gen ok.gen --sample ok.sample --make-bed --out tmp_gen_ok > /dev/null
cmp -s tmp_bg_ok.bed tmp_gen_ok.bed || bg_ok=1
rc=$(rc_of --bgen short.bgen --sample ok.sample --make-bed --out tmp_bg_short)
test "$rc" -gt 0 && test "$rc" -lt 128 || bg_ok=1
grep -q 'Invalid compressed SNP block' tmp_bg_short.log || bg_ok=1
check "--bgen v1.1 rejects a short SNP block" $bg_ok

# --- VCF ---

# Sites with 10 or more ALT alleles use a separate parser.  It ran the GP
# check even when there was no GP field, so nearly every call became missing;
# it dereferenced NULL when a sample's field was shorter than the GQ
# position; and with --vcf-half-call haploid it set ALT for a half-call of a
# different ALT allele.  The site has 11 ALT alleles, and PLINK 1.9 keeps the
# most common one (GA, allele 11); calls involving other ALT alleles become
# missing.
vcf_header() {
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
    printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">\n'
    printf '##FORMAT=<ID=GP,Number=G,Type=Float,Description="p">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    for i in 1 2 3 4 5 6 7 8; do printf '\tS%s' $i; done
    printf '\n'
}
alts="C,G,T,AC,AG,AT,CA,CC,CG,CT,GA"
{ vcf_header; printf '1\t200\tmu\tA\t%s\t.\tPASS\t.\tGT\t11/11\t11/1\t1/1\t0/11\t11/11\t1/1\t11/0\t0/0\n' $alts; } > tmp_many.vcf
{ vcf_header; printf '1\t200\tmu\tA\tGA\t.\tPASS\t.\tGT\t1/1\t./.\t./.\t0/1\t1/1\t./.\t1/0\t0/0\n'; } > tmp_many_ref.vcf
plink --vcf tmp_many.vcf --make-bed --out tmp_many > /dev/null
plink --vcf tmp_many_ref.vcf --make-bed --out tmp_many_ref > /dev/null
cmp -s tmp_many.bed tmp_many_ref.bed && cmp -s tmp_many.bim tmp_many_ref.bim
check "--vcf site with 11 ALT alleles and no GP field" $?

# S3's field has no GQ value; it must be treated like GQ '.'.
awk 'BEGIN { FS = OFS = "\t" } /^1\t200/ { $9 = "GT:GQ"; for (k = 10; k <= NF; ++k) if (k != 12) $k = $k ":50" } { print }' tmp_many.vcf > tmp_gq.vcf
awk 'BEGIN { FS = OFS = "\t" } /^1\t200/ { $12 = $12 ":." } { print }' tmp_gq.vcf > tmp_gq_ref.vcf
gq_ok=0
test "$(rc_of --vcf tmp_gq.vcf --vcf-min-gq 20 --make-bed --out tmp_gq)" -eq 0 || gq_ok=1
plink --vcf tmp_gq_ref.vcf --vcf-min-gq 20 --make-bed --out tmp_gq_ref > /dev/null
cmp -s tmp_gq.bed tmp_gq_ref.bed || gq_ok=1
check "--vcf-min-gq with a sample field shorter than the GQ position" $gq_ok

# A './1' half-call must be loaded like the haploid call '1' (allele 1 is
# not the one kept, so it is missing).  The GP field and --vcf-min-gp make
# the old code reach the half-call path.
gp=$(awk 'BEGIN { for (k = 0; k < 78; ++k) printf "%s1", (k? "," : "") }')
{
    vcf_header
    printf '1\t200\tmu\tA\t%s\t.\tPASS\t.\tGT:GP' $alts
    for g in 11/11 11/0 ./1 0/11 11/11 0/0 11/0 0/0; do printf '\t%s:%s' $g "$gp"; done
    printf '\n'
} > tmp_hc.vcf
sed "s|	\./1:|	1:|" tmp_hc.vcf > tmp_hc_ref.vcf
plink --vcf tmp_hc.vcf --vcf-half-call haploid --vcf-min-gp 0.5 --make-bed --out tmp_hc > /dev/null
plink --vcf tmp_hc_ref.vcf --vcf-min-gp 0.5 --make-bed --out tmp_hc_ref > /dev/null
cmp -s tmp_hc.bed tmp_hc_ref.bed && cmp -s tmp_hc.bim tmp_hc_ref.bim
check "--vcf-half-call haploid on a site with 11 ALT alleles" $?

# --biallelic-only list only opened the .skip.3allele file if no variant had
# been skipped yet, so a 3-allele variant after a QUAL-filtered one crashed.
# (PLINK 1.9 names the file <out>-temporary.skip.3allele.)
{
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n'
    printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
    printf '1\t100\tlowq\tA\tC\t5\tPASS\t.\tGT\t0/1\t1/1\n'
    printf '1\t200\ttri\tA\tC,G\t50\tPASS\t.\tGT\t0/2\t1/1\n'
    printf '1\t300\tok\tA\tC\t50\tPASS\t.\tGT\t0/1\t0/0\n'
} > tmp_sk.vcf
rm -f tmp_sk-temporary.skip.3allele
sk_ok=0
test "$(rc_of --vcf tmp_sk.vcf --vcf-min-qual 10 --biallelic-only list --make-bed --out tmp_sk)" -eq 0 || sk_ok=1
test "$(cat tmp_sk-temporary.skip.3allele 2> /dev/null)" = "tri" || sk_ok=1
test "$(cut -f2 tmp_sk.bim 2> /dev/null)" = "ok" || sk_ok=1
check "--biallelic-only list after a skipped variant" $sk_ok

# --- BCF ---

# Each .bcf from make_bcf.py must load like the matching .vcf.
bcf_matches_vcf() {
    plink --bcf "$1.bcf" $2 --make-bed --out "tmp_b_$1" > /dev/null 2>&1 || return 1
    plink --vcf "$1.vcf" $2 --make-bed --out "tmp_v_$1" > /dev/null 2>&1 || return 1
    cmp -s "tmp_b_$1.bim" "tmp_v_$1.bim" && cmp -s "tmp_b_$1.bed" "tmp_v_$1.bed"
}
# An ID used by both INFO and FORMAT (DP) has one string dictionary entry,
# but was counted twice, shifting the GT index.
bcf_matches_vcf dup_id ""
check "--bcf with an ID shared by INFO and FORMAT" $?
# '##contig=<ID=1>' (no comma) is legal.
bcf_matches_vcf contig_id_only ""
check "--bcf with an ID-only contig line" $?
# With no samples, the allele buffer size divided by zero.
bcf_matches_vcf no_samples "--allow-no-samples"
check "--bcf with no samples" $?
# '##FILTER=<ID=PASS>' (no comma) was mis-parsed.  Either a clean error or
# the same result as the .vcf is acceptable.
rc=$(rc_of --bcf filter_no_comma.bcf --make-bed --out tmp_b_filter_no_comma)
if [[ "$rc" -eq 0 ]]; then
    bcf_matches_vcf filter_no_comma ""
else
    test "$rc" -lt 128
fi
check "--bcf with a FILTER line without a comma" $?

# --- --tucc ---

# make_pseudocontrols() sought to the wrong variants and indexed genotypes by
# the unfiltered variant index, so --tucc with excluded variants wrote the
# wrong genotypes.  The pseudo-case/control genotypes must match plink2
# --tucc.  10 trios with affected children; 4 of 30 variants excluded.
awk "$lcg"'
BEGIN {
  s = 19
  for (j = 0; j < 30; ++j) print "1 v" j " 0 " (1000 * (j + 1)) > "tmp_tu.map"
  for (f = 0; f < 10; ++f) {
    for (j = 0; j < 30; ++j) for (k = 0; k < 4; ++k) h[k, j] = (rnd() < 1 / 3)? "C" : "A"
    pa = "t" f " pa 0 0 1 1"; ma = "t" f " ma 0 0 2 1"; ch = "t" f " ch pa ma " (1 + f % 2) " 2"
    for (j = 0; j < 30; ++j) {
      u = int(rnd() * 2)
      m = 2 + int(rnd() * 2)
      pa = pa " " h[0, j] " " h[1, j]; ma = ma " " h[2, j] " " h[3, j]; ch = ch " " h[u, j] " " h[m, j]
    }
    print pa > "tmp_tu.ped"; print ma > "tmp_tu.ped"; print ch > "tmp_tu.ped"
  }
}'
printf 'v0\nv7\nv8\nv20\n' > tmp_tu_ex.txt
plink --file tmp_tu --make-bed --out tmp_tu > /dev/null
plink --bfile tmp_tu --exclude tmp_tu_ex.txt --tucc write-bed --out tmp_tu19 > /dev/null
plink --bfile tmp_tu19.tucc --recode --out tmp_tu19 > /dev/null
$plink2 --bfile tmp_tu --exclude tmp_tu_ex.txt --tucc --out tmp_tu2 > /dev/null
$plink2 --pfile tmp_tu2.tucc --export ped --out tmp_tu2 > /dev/null
awk -f ped_long.awk tmp_tu19.map tmp_tu19.ped | sort > tmp_tu19.txt
awk -f ped_long.awk tmp_tu2.map tmp_tu2.ped | sort > tmp_tu2.txt
test "$(wc -l < tmp_tu2.txt)" -eq 520 && cmp -s tmp_tu19.txt tmp_tu2.txt
check "--tucc with excluded variants matches plink2" $?

exit $failed
