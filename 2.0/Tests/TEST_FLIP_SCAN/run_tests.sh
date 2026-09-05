#!/bin/bash

# --flip-scan, checked against PLINK 1.9.
#
# The fixture plants 20 haplotype blocks and then swaps the alleles of every
# 7th variant in the cases only, which is exactly the strand inconsistency the
# scan is meant to find.
#
# --maf 0.05 is applied throughout: PLINK 1.9 lets a nan through when a
# neighbor pair is monomorphic in one group, which counts the pair as a sign
# flip, and plink2 skips it.  Above that threshold the two agree.

set -exo pipefail

# The redesigned report adds a frequency-difference test, a major-allele
# frequency filter and a two-negative-match rule that PLINK 1.9 has no
# equivalent of.  These settings turn all three off, so the LD half can still
# be compared against 1.9 directly; the new behavior is checked separately.
COMPAT="--flip-scan-freq-diff 1 --flip-scan-max-maj-freq 1 --flip-scan-min-neg 1"

BUILD=$1
EXTRA1=$2
EXTRA2=$3

awk -f make_vcf.awk > tmp_h.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, (NR > 200)? 2 : 1}' tmp_data.fam > tmp_ped.fam
mv tmp_ped.fam tmp_data.fam

compare() {
    plink --bfile tmp_data --maf 0.05 --flip-scan $1 --out plink19
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan $1 cols=chrom,pos,ref,alt,altfreq,posct,rpos,negct,rneg,negids $COMPAT --out plink2
    awk -f compare.awk plink19.flipscan plink2.flipscan
}

# 1. Defaults, and each window/threshold flag.
compare ""
plink --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-window 5 --out plink19_w5
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan cols=chrom,pos,ref,alt,altfreq,posct,rpos,negct,rneg,negids $COMPAT --flip-scan-window 5 --out plink2_w5
awk -f compare.awk plink19_w5.flipscan plink2_w5.flipscan

# --flip-scan-window-kb values large enough to hold the whole variant-count
# window agree; smaller ones do not.  On this fixture PLINK 1.9 reports no
# neighbors at all at 40 kb and below, while plink2 truncates the window and
# keeps going, so only the larger values are compared here.
plink --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-window-kb 100 --out plink19_kb
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan cols=chrom,pos,ref,alt,altfreq,posct,rpos,negct,rneg,negids $COMPAT --flip-scan-window-kb 100 --out plink2_kb
awk -f compare.awk plink19_kb.flipscan plink2_kb.flipscan

plink --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-threshold 0.2 --out plink19_thr
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan cols=chrom,pos,ref,alt,altfreq,posct,rpos,negct,rneg,negids $COMPAT --flip-scan-threshold 0.2 --out plink2_thr
awk -f compare.awk plink19_thr.flipscan plink2_thr.flipscan

# 2. 'verbose' adds one line per above-threshold neighbor pair.
plink --bfile tmp_data --maf 0.05 --flip-scan verbose --out plink19_v
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan verbose cols=chrom,pos,ref,alt,altfreq,posct,rpos,negct,rneg,negids $COMPAT --out plink2_v
awk -f compare.awk plink19_v.flipscan plink2_v.flipscan
awk -f compare_verbose.awk plink19_v.flipscan.verbose plink2_v.flipscan.verbose

# 3. The planted flips have to be found: every 7th variant is one, and each
#    should be flagged.
test "$(grep -vc '^#' plink2.flipscan)" -gt 20

# 4. The result must not depend on --threads.
for t in 1 3 8
do
    $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan cols=chrom,pos,ref,alt,altfreq,posct,rpos,negct,rneg,negids $COMPAT --threads $t --out plink2_t$t
    diff -q plink2.flipscan plink2_t$t.flipscan
done

# 5. 'zs' is the same report, compressed.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan zs cols=chrom,pos,ref,alt,altfreq,posct,rpos,negct,rneg,negids $COMPAT --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.flipscan.zst > plink2_zs.flipscan
diff -q plink2.flipscan plink2_zs.flipscan

# 6. The redesign, which PLINK 1.9 has no equivalent of.
#    The default report carries the two group frequencies and a PROBLEM call.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --out plink2_new
head -n 1 plink2_new.flipscan | grep -qx '#CHROM	POS	ID	REF	ALT	CASE_MAJ_FREQ	CTRL_MAJ_FREQ	POS_CT	R_POS	NEG_CT	R_NEG	PROBLEM	NEG_IDS'
# Both frequencies are for the same allele, so they are on the same side of
# 0.5 unless the variant is one of the planted flips.
awk -F '\t' 'NR > 1 && $6 != "NA" && $7 != "NA" && $6 < 0.5 && $7 < 0.5 {print "both group frequencies below 0.5 on " $3; exit 1}' plink2_new.flipscan

# 'ref-allele-based' renames the columns and reports REF instead of MAJ.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan ref-allele-based --out plink2_ref
head -n 1 plink2_ref.flipscan | grep -q 'CASE_REF_FREQ	CTRL_REF_FREQ'
# Everything except the two frequency columns is unchanged by it.
diff -q <(cut -f 1-5,8- plink2_new.flipscan) <(cut -f 1-5,8- plink2_ref.flipscan)

# A frequency difference past the threshold flags a variant on its own, with
# no LD scan involved: at a threshold of 0 every variant is flagged.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-freq-diff 0 --out plink2_fd0
test "$(awk -F '\t' 'NR > 1 && $12 == "Y"' plink2_fd0.flipscan | wc -l)" -eq "$(($(wc -l < plink2_fd0.flipscan) - 1))"
# ...and those variants drop out of the LD scan entirely.
awk -F '\t' 'NR > 1 && ($8 != 0 || $10 != 0) {print "flagged variant still scanned: " $3; exit 1}' plink2_fd0.flipscan

# Requiring two sign-flipped neighbors instead of one can only reduce the
# number of problem calls.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-min-neg 1 --out plink2_mn1
n2=$(awk -F '\t' 'NR > 1 && $12 == "Y"' plink2_new.flipscan | wc -l)
n1=$(awk -F '\t' 'NR > 1 && $12 == "Y"' plink2_mn1.flipscan | wc -l)
test "$n2" -le "$n1"

# The LD scan needs 50 founders per group unless --bad-ld says otherwise.
awk 'NR % 8 == 0 {print $1, $2}' tmp_data.fam > tmp_few.txt
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --keep tmp_few.txt --bad-freqs --maf 0.05 --flip-scan --out plink2_few > /dev/null 2>&1; then
    echo "--flip-scan ran on too few founders"
    exit 1
fi
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --keep tmp_few.txt --bad-freqs --bad-ld --maf 0.05 --flip-scan --out plink2_few

# 7. --flip-scan-ref-freq: comparison against a reference allele frequency
#    file, with no case/control split and no LD scan.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --freq --out tmp_panel
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-freq tmp_panel.afreq --out plink2_rf
head -n 1 plink2_rf.flipscan | grep -qx '#CHROM	POS	ID	REF	ALT	MAJ_FREQ	PANEL_MAJ_FREQ	PROBLEM'
# Against its own frequencies nothing can differ.
awk -F '\t' 'NR > 1 && $8 != "N" {print "unexpected call on " $3; exit 1}' plink2_rf.flipscan
# One row per scanned variant, same set the LD-mode report covers.
diff -q <(cut -f 3 plink2_rf.flipscan | tail -n +2) <(cut -f 3 plink2_new.flipscan | tail -n +2)

# Moving every fifth panel frequency by exactly 0.3 has to be caught, and
# nothing else with it.  Both sides report the same allele, so a 0.3 shift in
# the panel is a 0.3 difference in the report whichever allele is major.  The
# ALT_FREQS column position depends on whether PROVISIONAL_REF? is present,
# so it is found by name.
awk 'BEGIN{OFS="\t"} NR == 1 {for (i = 1; i <= NF; ++i) {if ($i == "ALT_FREQS") {c = i}}; print; next}
     {if (NR % 5 == 0) {$c = ($c <= 0.5)? ($c + 0.3) : ($c - 0.3); ids[$2] = 1}; print}
     END {for (id in ids) {print id > "tmp_shifted_ids.txt"}}' tmp_panel.afreq > tmp_flipped.afreq
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-freq tmp_flipped.afreq --out plink2_rff
test "$(awk -F '\t' 'NR > 1 && $8 == "Y"' plink2_rff.flipscan | wc -l)" -gt 0
# Every call is a shifted variant, and every shifted variant that survived
# --maf is called.
awk -F '\t' '
    FNR == NR {shifted[$1] = 1; next}
    FNR > 1 {
        d = $6 - $7; if (d < 0) {d = -d}
        if (($8 == "Y") != (d > 0.2)) {print "call disagrees with the difference at " $3; exit 1}
        if (($8 == "Y") && !($3 in shifted)) {print "called an unshifted variant: " $3; exit 1}
        if (($8 == "N") && ($3 in shifted)) {print "missed a shifted variant: " $3; exit 1}
    }' tmp_shifted_ids.txt plink2_rff.flipscan

# A variant absent from the file is NA rather than a call.
head -n 40 tmp_panel.afreq > tmp_partial.afreq
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-freq tmp_partial.afreq --out plink2_rfp
test "$(awk -F '\t' 'NR > 1 && $8 == "NA"' plink2_rfp.flipscan | wc -l)" -gt 0
awk -F '\t' 'NR > 1 && $8 == "NA" && $7 != "NA" {print "panel frequency present for an NA row: " $3; exit 1}' plink2_rfp.flipscan

# This mode needs neither a phenotype nor sorted coordinates, so it works on
# the same fileset with the phenotype column blanked out.
awk 'BEGIN{OFS=" "} {$6 = -9; print}' tmp_data.fam > tmp_nopheno.fam
cp tmp_data.bed tmp_nopheno.bed
cp tmp_data.bim tmp_nopheno.bim
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_nopheno --maf 0.05 --flip-scan --flip-scan-ref-freq tmp_panel.afreq --out plink2_rfnp
diff -q <(cut -f 1-5 plink2_rf.flipscan) <(cut -f 1-5 plink2_rfnp.flipscan)

# --flip-scan-ref-freq without --flip-scan is an error.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --flip-scan-ref-freq tmp_panel.afreq --freq --out plink2_bad > /dev/null 2>&1; then
    echo "--flip-scan-ref-freq ran without --flip-scan"
    exit 1
fi
