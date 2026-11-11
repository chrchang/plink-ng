#!/bin/bash

# Regression tests for PLINK 1.9 command-line and output bugs.  Each check
# failed before its fix.
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
rc_of() {
    plink "$@" --out tmp_out > tmp_out.txt 2>&1 < /dev/null
    echo $?
}

plink --dummy 8 4 --seed 1 --make-bed --out tmp_data > /dev/null

# Invalid --recode parameter must be rejected, not ignored.
test "$(rc_of --bfile tmp_data --recode hvjunk)" -eq 5
check "--recode hvjunk is rejected" $?

# --set-r2-phase was swallowed by the --set-r2 prefix match.
rc_of --bfile tmp_data --set-r2-phase > /dev/null
grep -q 'set-r2-phase is provisionally retired' tmp_out.txt
check "--set-r2-phase is recognized" $?

# 'mat' was never reachable, since the 'parentdt2' test was duplicated.
rc_of --bfile tmp_data --tdt poo mat > /dev/null
! grep -q "Invalid --tdt parameter" tmp_out.txt
check "--tdt poo mat is recognized" $?

cnv_rc=$(rc_of --cnv-list tmp_nonexistent.cnv --cnv-make-map)
# Stable builds (-DSTABLE_BUILD, the 1.9 default) reject --cnv-list up front.
if grep -q 'is either unfinished' tmp_out.txt; then
    echo "skipped: --cnv-list checks (stable build)"
else
    # --cnv-make-map discarded its error code.
    test "$cnv_rc" -ne 0
    check "--cnv-make-map failure gives a nonzero exit code" $?
    # --cnv-list copied its file name without a terminating null.
    grep -qx 'Error: Failed to open tmp_nonexistent.cnv.' tmp_out.txt
    check "--cnv-list file name is printed intact" $?
fi

# These checks were skipped or had the wrong exit code.
test "$(rc_of --bfile tmp_data --read-dists tmp_nonexistent.dist)" -eq 5
check "--read-dists requires a command using it" $?
test "$(rc_of --bfile tmp_data --meta-analysis-report-dups)" -eq 5
check "--meta-analysis-report-dups without --meta-analysis is a command-line error" $?

# --flip-subset with --indiv-sort flipped the wrong sample.
printf '1 b 0 0 1 1\n1 c 0 0 1 1\n1 a 0 0 1 1\n' > tmp_sort.fam
printf '1\trs1\t0\t100\tT\tA\n' > tmp_sort.bim
printf '\154\033\001\062' > tmp_sort.bed
printf 'rs1\n' > tmp_flip.txt
printf '1 a\n' > tmp_sub.txt
plink --bfile tmp_sort --indiv-sort ascii --flip tmp_flip.txt --flip-subset tmp_sub.txt --make-bed --out tmp_sorted > /dev/null
test "$(od -An -tx1 -j3 tmp_sorted.bed | tr -d ' \n')" = "08"
check "--flip-subset with --indiv-sort flips the right sample" $?

# --recode HV with multi-character alleles corrupted genotypes.
awk 'BEGIN {
  s = 12345
  a = "AAAAAAAAAAAA"; g = "GGGGGGGGGGGGGGGG"
  for (i = 0; i < 8; ++i) {
    line = "f" i " i" i " 0 0 1 1"
    for (j = 0; j < 24; ++j) {
      s = (s * 1103515245 + 12345) % 2147483648
      line = line " " ((int(s / 1024) % 2)? a : g)
    }
    print line
  }
}' > tmp_hv.ped
awk 'BEGIN { for (j = 1; j <= 12; ++j) print "1\trs" j "\t0\t" (j * 1000) }' > tmp_hv.map
plink --file tmp_hv --recode HV --out tmp_hv1 > /dev/null
plink --file tmp_hv --recode --out tmp_hv2 > /dev/null
test "$(cut -d' ' -f7- tmp_hv1.chr-1.ped)" = "$(cut -d' ' -f7- tmp_hv2.ped)"
check "--recode HV keeps multi-character alleles" $?


# The checks below cover fixes that change numeric output.  Random data uses a
# 32-bit LCG small enough to be exact in awk's doubles.

# --model --adjust labeled results with the wrong variant IDs after a skipped
# chrY/MT variant.  The .map is sorted by chromosome code, so the XY variants
# end up after the Y variants.
awk 'BEGIN {
  s = 7
  for (i = 0; i < 60; ++i) {
    sex = 1 + (i % 2)
    line = "f" i " i" i " 0 0 " sex " " (1 + (int(i / 3) % 2))
    for (j = 0; j < 10; ++j) {
      s = (s * 69069 + 1) % 4294967296
      r = int(s / 65536) % 4
      if (j == 3 || j == 4) {
        g = (sex == 2)? "0 0" : ((r < 2)? "A A" : "C C")
      } else if (j == 7) {
        g = (r < 2)? "A A" : "C C"
      } else {
        g = (r == 0)? "A A" : ((r == 3)? "C C" : "A C")
      }
      line = line " " g
    }
    print line
  }
}' > tmp_adj.ped
printf '1 v1 0 100\n1 v2 0 200\n1 v3 0 300\nY y1 0 100\nY y2 0 200\nXY v4 0 400\nXY v5 0 500\nMT m1 0 100\nXY v6 0 600\nXY v7 0 700\n' > tmp_adj.map
plink --file tmp_adj --model trend-only --adjust --out tmp_adj > /dev/null
awk 'NR == FNR { if (FNR > 1) p[$2] = $NF; next }
  FNR > 1 { ++n; d = $3 - p[$2]; if (!($2 in p) || d > 1e-6 || d < -1e-6) bad = 1 }
  END { exit (bad || n != 7) }' tmp_adj.model tmp_adj.model.trend.adjusted
check "--model --adjust reports each variant's own p-value" $?

# --ld's HWE printout on chrX subtracted males incorrectly.  These genotypes
# have multiple phasing solutions, so --ld prints HWE p-values; they must
# match --hardy, which uses only females on chrX.
{
  printf 'm1 m1 0 0 1 1 A A A A\nm2 m2 0 0 1 1 A A C C\n'
  printf 'm3 m3 0 0 1 1 C C A A\nm4 m4 0 0 1 1 C C A A\n'
  printf 'm5 m5 0 0 1 1 C C C C\nm6 m6 0 0 1 1 C C C C\n'
  printf 'w1 w1 0 0 2 1 A A A A\nw2 w2 0 0 2 1 A C A A\n'
  for i in $(seq 3 18); do printf 'w%s w%s 0 0 2 1 A C A C\n' $i $i; done
} > tmp_ldx.ped
printf 'X x1 0 100\nX x2 0 200\n' > tmp_ldx.map
plink --file tmp_ldx --ld x1 x2 --out tmp_ldx > /dev/null
plink --file tmp_ldx --hardy --out tmp_ldx_hardy > /dev/null
grep -A3 '^HWE exact test' tmp_ldx.log | awk '/: / { sub(":", "", $1); print $1, $2 }' > tmp_ldx_hwe.txt
awk 'NR == FNR { p[$1] = $2; next }
  FNR > 1 && $3 == "ALL" { ++n; d = ($NF - p[$2]) / $NF; if (!($2 in p) || d > 1e-3 || d < -1e-3) bad = 1 }
  END { exit (bad || n != 2) }' tmp_ldx_hwe.txt tmp_ldx_hardy.hwe
check "--ld HWE p-values on chrX match --hardy" $?

# --epistasis seeded each thread's best-chisq buffer from the wrong cells,
# which corrupted the .summary BEST_CHISQ column after the first pass.
# --memory 64 forces several passes over these variants.  Each BEST_CHISQ
# must equal the largest reported STAT for that variant, or be below the
# reporting threshold (chisq 6.635 for --epi1 0.01) if none was reported.
epi_summary_ok() {
    awk 'FNR == 1 { next }
      NR == FNR {
        s = $(NF - 1) + 0
        if (s > best[$2]) best[$2] = s
        if (s > best[$4]) best[$4] = s
        next
      }
      {
        if ($2 in best) {
          d = $6 - best[$2]
          if (d > 0.001 * best[$2] || d < -0.001 * best[$2]) bad = 1
        } else if ($6 > 6.635) {
          bad = 1
        }
      }
      END { exit bad }' "$1" "$1.summary"
}
plink --dummy 60 2500 0.02 scalar-pheno --seed 5 --make-bed --out tmp_epiq > /dev/null
plink --bfile tmp_epiq --epistasis --epi1 0.01 --memory 64 --threads 4 --out tmp_epiq > /dev/null
epi_summary_ok tmp_epiq.epi.qt
check "--epistasis (linear) .summary matches the pair report" $?
plink --dummy 60 3000 0.02 --seed 5 --make-bed --out tmp_epic > /dev/null
plink --bfile tmp_epic --epistasis --epi1 0.01 --memory 64 --threads 4 --out tmp_epic > /dev/null
epi_summary_ok tmp_epic.epi.cc
check "--epistasis (logistic) .summary matches the pair report" $?

# --linear dominant with a sex interaction used additive coding for the
# interaction column.  Compare with an explicit covariate model
# (dominant code D, sex S = 1 for males and 0 for females, and D*S).
plink --dummy 300 3 0.02 scalar-pheno --seed 11 --make-bed --out tmp_dom > /dev/null
awk 'BEGIN { OFS = "\t" } { $5 = 1 + (NR % 2); print }' tmp_dom.fam > tmp_dom_sex.fam
mv tmp_dom_sex.fam tmp_dom.fam
plink --bfile tmp_dom --linear dominant sex interaction --out tmp_dom > /dev/null
plink --bfile tmp_dom --recode A --out tmp_dom > /dev/null
dom_ok=0
for j in 0 1 2; do
    awk -v c=$((j + 7)) 'NR == 1 { print "FID IID D S DS"; next }
      $c == "NA" { print $1, $2, -9, -9, -9; next }
      { d = ($c > 0)? 1 : 0; s = ($5 == 1)? 1 : 0; print $1, $2, d, s, d * s }' tmp_dom.raw > tmp_dom_covar.txt
    plink --bfile tmp_dom --linear no-snp --covar tmp_dom_covar.txt --out tmp_dom_ref > /dev/null
    awk -v snp=snp$j '$2 == snp { print $7, $8 }' tmp_dom.assoc.linear > tmp_dom_a.txt
    awk 'NR > 1 { print $3, $4 }' tmp_dom_ref.assoc.linear > tmp_dom_b.txt
    if ! cmp -s tmp_dom_a.txt tmp_dom_b.txt; then
        dom_ok=1
    fi
done
check "--linear dominant sex interaction matches explicit covariates" $dom_ok

# --linear standard-beta did not drop every sample with a missing genotype
# from the phenotype standardization.  The missing calls are placed after the
# first 64 samples, and the result must match a run with them removed.
awk 'BEGIN {
  s = 17
  for (i = 0; i < 100; ++i) {
    s = (s * 69069 + 1) % 4294967296
    r = int(s / 65536) % 3
    s = (s * 69069 + 1) % 4294967296
    q = (int(s / 65536) % 1000) / 100.0 + 1.5 * r
    g = (i >= 70)? "0 0" : ((r == 0)? "A A" : ((r == 1)? "A C" : "C C"))
    print "f" i, "i" i, 0, 0, 1 + (i % 2), q, g
  }
}' > tmp_sb.ped
echo "1 v1 0 100" > tmp_sb.map
awk '$7 == "0" { print $1, $2 }' tmp_sb.ped > tmp_sb_miss.txt
plink --file tmp_sb --linear standard-beta --out tmp_sb1 > /dev/null
plink --file tmp_sb --remove tmp_sb_miss.txt --linear standard-beta --out tmp_sb2 > /dev/null
sb1=$(awk 'NR == 2 { print $6, $7, $8, $9 }' tmp_sb1.assoc.linear)
test -n "$sb1" && test "$sb1" = "$(awk 'NR == 2 { print $6, $7, $8, $9 }' tmp_sb2.assoc.linear)"
check "--linear standard-beta ignores samples with missing genotypes" $?

# --write-set printed the wrong variant IDs for sets stored as a bitfield,
# which starts at a multiple of 128 variants.
plink --dummy 4 600 --seed 1 --make-bed --out tmp_ws > /dev/null
{ echo S1; for i in $(seq 300 2 398); do echo snp$i; done; echo END; } > tmp_ws_in.set
plink --bfile tmp_ws --set tmp_ws_in.set --write-set --out tmp_ws > /dev/null
test "$(grep -v '^$' tmp_ws.set)" = "$(cat tmp_ws_in.set)"
check "--write-set lists a sparse set's own variants" $?

# --show-tags and --flip-scan don't implement --ld-xchr 3, and silently
# treated it as mode 1.  It must be rejected, whatever the flag order.
printf 'X\tx1\t0\t100\tA\tC\nX\tx2\t0\t200\tA\tC\n' > tmp_xc.bim
printf 'f1 i1 0 0 1 2\nf2 i2 0 0 2 1\nf3 i3 0 0 2 2\nf4 i4 0 0 1 1\n' > tmp_xc.fam
printf '\154\033\001\000\000' > tmp_xc.bed
xc_ok=0
test "$(rc_of --bfile tmp_xc --show-tags all --ld-xchr 3)" -eq 5 || xc_ok=1
test "$(rc_of --ld-xchr 3 --bfile tmp_xc --flip-scan)" -eq 5 || xc_ok=1
grep -q 'flip-scan + --ld-xchr 3 is not currently supported' tmp_out.txt || xc_ok=1
check "--show-tags and --flip-scan reject --ld-xchr 3" $xc_ok

# These flags were matched with memcmp() lengths that left out the
# terminating NUL, so any longer flag starting with the same text (here, a
# misspelling) was silently taken as the flag itself.  (--set-r2 is covered
# above.)
printf 'S1\nsnp1\nsnp2\nEND\n' > tmp_cl.set
printf 'snp1 1\nsnp2 0\n' > tmp_cl_tags.txt
cl_ok=0
for args in "--cluster --mccx 2 2" "--set tmp_cl.set --set-maxx 3 --write-set" \
            "--show-tags tmp_cl_tags.txt --tag-mode2x" "--set tmp_cl.set --write-set-r2x"; do
    test "$(rc_of --bfile tmp_data $args)" -eq 5 || cl_ok=1
    grep -q 'Unrecognized flag' tmp_out.txt || cl_ok=1
done
# The correctly spelled flags must still work.
for args in "--cluster --mcc 2 2" "--set tmp_cl.set --set-max 3 --write-set" \
            "--show-tags tmp_cl_tags.txt --tag-mode2"; do
    test "$(rc_of --bfile tmp_data $args)" -eq 0 || cl_ok=1
done
check "misspelled --mcc/--set-max/--tag-mode2/--write-set-r2 are rejected" $cl_ok

# --vif's condition had a misplaced negation, so it was accepted with any
# command instead of only --linear and --epistasis.  (The accepted cases use
# the seeded quantitative-trait data: --epistasis on tmp_data's 8 samples
# fails whenever fewer than two are cases or controls.)
plink --dummy 20 30 0.02 scalar-pheno --seed 1 --make-bed --out tmp_qt > /dev/null
vif_ok=0
test "$(rc_of --bfile tmp_data --assoc --vif 50)" -eq 5 || vif_ok=1
grep -q 'vif must be used with' tmp_out.txt || vif_ok=1
test "$(rc_of --bfile tmp_data --fast-epistasis --vif 50)" -eq 5 || vif_ok=1
test "$(rc_of --bfile tmp_qt --linear --vif 50)" -eq 0 || vif_ok=1
test "$(rc_of --bfile tmp_qt --epistasis --vif 50)" -eq 0 || vif_ok=1
check "--vif is only accepted with --linear/--epistasis" $vif_ok

# --cnv-write freq / --cnv-write-freq checked the --cnv-freq-method2 value
# instead of its flag bits, so they were rejected even with
# --cnv-freq-method2.  CNV filesets are disabled in stable builds, so this
# only runs on development builds.
printf 'f1 i1 0 0 1 2\nf2 i2 0 0 2 1\nf3 i3 0 0 1 2\nf4 i4 0 0 2 1\n' > tmp_cnv.fam
printf 'FID IID CHR BP1 BP2 TYPE SCORE SITES\nf1 i1 1 1000 5000 1 10 5\nf2 i2 1 2000 6000 3 10 5\nf3 i3 1 1500 4000 1 10 5\nf4 i4 2 1000 3000 0 10 5\n' > tmp_cnv.cnv
printf '1\tp1-1000\t0\t1000\n1\tp1-1500\t0\t1500\n1\tp1-2000\t0\t2000\n1\tp1-4000\t0\t4000\n1\tp1-5000\t0\t5000\n1\tp1-6000\t0\t6000\n2\tp2-1000\t0\t1000\n2\tp2-3000\t0\t3000\n' > tmp_cnv.cnv.map
rc_of --cfile tmp_cnv --cnv-write freq --cnv-freq-method2 0.5 > /dev/null
if grep -q 'is either unfinished' tmp_out.txt; then
    echo "skipped: --cnv-write freq is accepted with --cnv-freq-method2 (stable build)"
else
    cnv_ok=0
    test "$(rc_of --cfile tmp_cnv --cnv-write freq --cnv-freq-method2 0.5)" -eq 0 || cnv_ok=1
    test "$(rc_of --cfile tmp_cnv --cnv-write --cnv-write-freq --cnv-freq-method2 0.5)" -eq 0 || cnv_ok=1
    test "$(rc_of --cfile tmp_cnv --cnv-write freq)" -eq 5 || cnv_ok=1
    check "--cnv-write freq is accepted with --cnv-freq-method2" $cnv_ok
fi

exit $failed
