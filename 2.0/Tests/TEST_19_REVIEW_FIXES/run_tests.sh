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

plink --dummy 8 4 --make-bed --out tmp_data > /dev/null

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

exit $failed
