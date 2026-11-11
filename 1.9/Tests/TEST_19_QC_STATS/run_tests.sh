#!/bin/bash

# PLINK 1.9 sample and variant QC statistics and filters, checked against
# plink2: --het, --check-sex, --hardy, --r2/--r, --freq, --maf/--max-maf,
# --nonfounders, --filter-founders, --filter-females and --write-covar.
# (Ported from the old 1.9/tests_old suite, which compared against PLINK
# 1.07.)
#
# Usage: ./run_tests.sh {plink 1.9 dir} {plink2 build dir} {up to 2 plink2 args}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  {plink2 build dir} defaults to
# ../../../2.0/build_dynamic.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}
plink2="${2:-../../../2.0/build_dynamic}/plink2 ${3:-} ${4:-}"
compare_cols="../compare_cols.awk"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}

# 513 samples (513 % 4 == 1), alternating sex, 400 autosomal and 200 chrX
# variants, 3% missing.  The first 40 samples are children of samples 102-181.
# Males' heterozygous chrX calls are set to missing, since plink2 has no
# heterozygous haploid state.
plink --dummy 513 600 0.03 acgt --seed 5 --make-bed --out tmp_raw > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; if (NR > 400) $1 = 23; print }' tmp_raw.bim > tmp_raw.bim.new
mv tmp_raw.bim.new tmp_raw.bim
awk '{ if (NR <= 40) { $3 = "per" (100 + 2 * NR); $4 = "per" (101 + 2 * NR) } $5 = 1 + NR % 2; print }' tmp_raw.fam > tmp_raw.fam.new
mv tmp_raw.fam.new tmp_raw.fam
plink --bfile tmp_raw --set-hh-missing --make-bed --out tmp_data > /dev/null 2>&1

# --het (autosomes only in both).
plink --bfile tmp_data --het --out tmp_het19 > /dev/null
$plink2 --bfile tmp_data --autosome --het --out tmp_het2 > /dev/null
awk -v c1=1,2,3,4,5,6 -v c2=1,2,3,4,5,6 -f $compare_cols tmp_het19.het tmp_het2.het
check "--het matches plink2" $?

# --check-sex, with default and explicit thresholds.  Three males are
# recorded as female, and one as unknown.  (Males have no heterozygous chrX
# calls, so this doesn't create heterozygous haploid calls.)  plink2 prints
# an unknown PEDSEX as NA.
awk '{ if (NR == 50 || NR == 52 || NR == 60) $5 = 2; if (NR == 70) $5 = 0; print }' tmp_data.fam > tmp_sex.fam
sexcheck2() {
    awk '{ if ($3 == "NA") $3 = 0; print }' "$1"
}
plink --bed tmp_data.bed --bim tmp_data.bim --fam tmp_sex.fam --check-sex --out tmp_cs19 > /dev/null 2>&1
$plink2 --bed tmp_data.bed --bim tmp_data.bim --fam tmp_sex.fam --check-sex --out tmp_cs2 > /dev/null 2>&1
awk -v c1=1,2,3,4,5,6 -v c2=1,2,3,4,5,6 -f $compare_cols tmp_cs19.sexcheck <(sexcheck2 tmp_cs2.sexcheck)
check "--check-sex matches plink2" $?
plink --bed tmp_data.bed --bim tmp_data.bim --fam tmp_sex.fam --check-sex 0.3 0.7 --out tmp_cs19t > /dev/null 2>&1
$plink2 --bed tmp_data.bed --bim tmp_data.bim --fam tmp_sex.fam --check-sex max-female-xf=0.3 min-male-xf=0.7 --out tmp_cs2t > /dev/null
awk -v c1=1,2,3,4,5,6 -v c2=1,2,3,4,5,6 -f $compare_cols tmp_cs19t.sexcheck <(sexcheck2 tmp_cs2t.sexcheck) && test "$(grep -c PROBLEM tmp_cs19t.sexcheck)" -eq 4
check "--check-sex with thresholds matches plink2" $?

# --hardy.  Genotype counts are printed as (smaller hom, het, larger hom), since
# the two versions pick A1 differently.  PLINK 1.9's ALL rows must match
# plink2 on everything, and its AFF/UNAFF rows must match plink2 run on just
# the cases or controls.  On chrX, PLINK 1.9 tests only females, while plink2
# also uses males, so only the female counts and heterozygosities are
# compared there.
hwe19() {
    awk -v t="$2" -v x="$3" 'NR > 1 && $3 == t && (($1 == 23) == x) {
      split($6, g, "/"); lo = (g[1] < g[3])? g[1] : g[3]; hi = g[1] + g[3] - lo
      print $2, lo, g[2], hi, $7, $8, $9 }' "$1"
}
hwe2() {
    awk 'NR > 1 { lo = ($5 < $7)? $5 : $7; hi = $5 + $7 - lo
      print $2, lo, $6, hi, (FILENAME ~ /\.x$/)? $10 " " $11 " NA" : $8 " " $9 " " $10 }' "$1"
}
plink --bfile tmp_data --hardy --out tmp_hwe19 > /dev/null
awk '$6 == 2 { print $1, $2 }' tmp_data.fam > tmp_cases.txt
awk '$6 == 1 { print $1, $2 }' tmp_data.fam > tmp_controls.txt
hardy_ok=0
for grp in ALL AFF UNAFF; do
    keep=
    if [[ $grp == AFF ]]; then
        keep="--keep tmp_cases.txt"
    elif [[ $grp == UNAFF ]]; then
        keep="--keep tmp_controls.txt"
    fi
    $plink2 --bfile tmp_data $keep --hardy --out tmp_hwe2_$grp > /dev/null
    awk -v c1=1,2,3,4,5,6,7 -v c2=1,2,3,4,5,6,7 -v skip1=0 -v skip2=0 -f $compare_cols \
        <(hwe19 tmp_hwe19.hwe $grp 0) <(hwe2 tmp_hwe2_$grp.hardy) || hardy_ok=1
    awk -v c1=1,2,3,4,5,6 -v c2=1,2,3,4,5,6 -v skip1=0 -v skip2=0 -f $compare_cols \
        <(hwe19 tmp_hwe19.hwe $grp 1) <(hwe2 tmp_hwe2_$grp.hardy.x) || hardy_ok=1
done
check "--hardy ALL/AFF/UNAFF rows match plink2" $hardy_ok

# --freq, and --maf/--max-maf, with and without --nonfounders.  (The .frq MAF
# column is the A1 frequency.  With --keep-allele-order, A1 is the .bim A1,
# which plink2 reports as ALT.)
freq_ok=0
for nf in "" --nonfounders; do
    plink --bfile tmp_data $nf --keep-allele-order --freq --out tmp_frq19 > /dev/null
    $plink2 --bfile tmp_data $nf --freq --out tmp_frq2 > /dev/null
    awk -v c1=2,3,5,6 -v c2=2,4,6,7 -f $compare_cols tmp_frq19.frq tmp_frq2.afreq || freq_ok=1
    plink --bfile tmp_data $nf --maf 0.45 --max-maf 0.49 --write-snplist --out tmp_mm19 > /dev/null
    $plink2 --bfile tmp_data $nf --maf 0.45 --max-maf 0.49 --write-snplist --out tmp_mm2 > /dev/null
    cmp -s tmp_mm19.snplist tmp_mm2.snplist || freq_ok=1
done
check "--freq and --maf/--max-maf, with and without --nonfounders, match plink2" $freq_ok
plink --bfile tmp_data --maf 0.45 --max-maf 0.49 --write-snplist --out tmp_mmf > /dev/null
plink --bfile tmp_data --nonfounders --maf 0.45 --max-maf 0.49 --write-snplist --out tmp_mmn > /dev/null
! cmp -s tmp_mmf.snplist tmp_mmn.snplist
check "--nonfounders changes which variants pass --maf/--max-maf" $?

# --covar --covar-number --filter-founders --write-covar.  C2 has missing values.
awk 'BEGIN { s = 9; print "FID IID C1 C2 C3 C4" }
  {
    s = (s * 69069 + 1) % 4294967296; a = int(s / 65536) % 100
    s = (s * 69069 + 1) % 4294967296; b = int(s / 65536) % 7
    print $1, $2, a / 10, (b == 0)? -9 : b, NR % 3, a % 5
  }' tmp_data.fam > tmp_covar.txt
plink --bfile tmp_data --covar tmp_covar.txt --covar-number 2,4 --filter-founders --write-covar --out tmp_wc19 > /dev/null
$plink2 --bfile tmp_data --covar tmp_covar.txt --covar-col-nums 4,6 --keep-founders --write-covar --out tmp_wc2 > /dev/null
cmp -s <(sed 's/ $//' tmp_wc19.cov | tail -n +2) <(tr '\t' ' ' < tmp_wc2.cov | sed 's/ NA/ -9/g' | tail -n +2) && test "$(wc -l < tmp_wc19.cov)" -eq 474
check "--covar-number --filter-founders --write-covar matches plink2" $?

# --r2 and --r (default 10-variant window, every pair reported), and --r2
# square, on the autosomes.
plink --bfile tmp_data --chr 1 --r2 --ld-window-r2 0 --out tmp_r2_19 > /dev/null
$plink2 --bfile tmp_data --chr 1 --r2-unphased --ld-window 10 --ld-window-r2 0 --out tmp_r2_2 > /dev/null
awk -v c1=1,2,3,4,5,6,7 -v c2=1,2,3,4,5,6,7 -f $compare_cols tmp_r2_19.ld tmp_r2_2.vcor && test "$(wc -l < tmp_r2_19.ld)" -gt 3000
check "--r2 matches plink2 --r2-unphased" $?
plink --bfile tmp_data --chr 1 --r --ld-window-r2 0 --out tmp_r_19 > /dev/null
$plink2 --bfile tmp_data --chr 1 --r-unphased --ld-window 10 --ld-window-r2 0 --out tmp_r_2 > /dev/null
awk -v c1=3,6,7 -v c2=3,7,9 -f $compare_cols tmp_r_19.ld tmp_r_2.vcor
check "--r matches plink2 --r-unphased" $?
plink --bfile tmp_data --chr 1 --r2 square --out tmp_sq19 > /dev/null
$plink2 --bfile tmp_data --chr 1 --r2-unphased square --out tmp_sq2 > /dev/null
cmp -s tmp_sq19.ld tmp_sq2.unphased.vcor2
check "--r2 square matches plink2 --r2-unphased square" $?

# --filter-females, --filter-controls
plink --bfile tmp_data --filter-females --filter-controls --make-just-fam --out tmp_ff19 > /dev/null
$plink2 --bfile tmp_data --keep-females --keep-if PHENO1 == control --make-just-psam --out tmp_ff2 > /dev/null
cmp -s <(cut -d ' ' -f 1,2 tmp_ff19.fam) <(tail -n +2 tmp_ff2.psam | cut -f 1,2 | tr '\t' ' ')
check "--filter-females --filter-controls matches plink2" $?

exit $failed
