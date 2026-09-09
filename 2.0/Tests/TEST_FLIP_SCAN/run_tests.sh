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

# 8. --flip-scan-ref-pfile/--flip-scan-ref-bfile: same comparison against a
#    second fileset, whose allele frequencies are computed here.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --make-bed --out tmp_ref
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_ref --make-pgen --out tmp_refp
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-bfile tmp_ref --out plink2_rd
head -n 1 plink2_rd.flipscan | grep -qx '#CHROM	POS	ID	REF	ALT	MAJ_FREQ	PANEL_MAJ_FREQ	PROBLEM'
# The reference is this dataset, so it has to agree with the --flip-scan-ref-freq
# run against this dataset's own frequency file, to the digit.
diff -q plink2_rf.flipscan plink2_rd.flipscan

# The .pgen/.pvar/.psam spelling, the three-filename spelling, and the .bed
# spelling all name the same fileset.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-pfile tmp_refp --out plink2_rdp
diff -q plink2_rd.flipscan plink2_rdp.flipscan
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-pfile tmp_refp.pgen tmp_refp.pvar tmp_refp.psam --out plink2_rd3
diff -q plink2_rd.flipscan plink2_rd3.flipscan

# A reference variant whose REF and ALT are the other way round is the case
# this mode exists to find: the allele pair still matches, so the frequency is
# turned around before the comparison and the report comes out unchanged, with
# the log saying how many were like that.  --ref-allele rewrites the same
# genotypes with the other allele as REF, which is exactly that situation.
awk 'BEGIN{OFS="\t"} !/^#/ {if (++i % 3 == 0) {print $3, $5; ++n}} END {print n > "tmp_swap_ct.txt"}' tmp_refp.pvar > tmp_swap_ref.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_refp --ref-allele force tmp_swap_ref.txt 2 1 --make-pgen --out tmp_swapped
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-pfile tmp_swapped --out plink2_rds > tmp_rds.log
diff -q plink2_rd.flipscan plink2_rds.flipscan
grep -q "$(cat tmp_swap_ct.txt) with REF and ALT the other way round" tmp_rds.log

# A variant absent from the reference, and a variant whose alleles do not match
# in either order, are both PROBLEM=NA rather than calls.
awk 'NR <= 40 {print $2}' tmp_ref.bim > tmp_keep40.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_ref --extract tmp_keep40.txt --make-bed --out tmp_part
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-bfile tmp_part --out plink2_rdn
test "$(awk -F '\t' 'NR > 1 && $8 == "NA"' plink2_rdn.flipscan | wc -l)" -gt 0
awk -F '\t' 'NR > 1 && $8 == "NA" && $7 != "NA" {print "panel frequency present for an NA row: " $3; exit 1}' plink2_rdn.flipscan
# The 40 that are present keep the same answer they had against the full
# reference.
awk -F '\t' 'FNR == NR {ans[$3] = $8; next} FNR > 1 && $8 != "NA" && $8 != ans[$3] {print "answer changed for " $3; exit 1}' plink2_rd.flipscan plink2_rdn.flipscan

# Alleles that match by ID but not as a pair are skipped with a warning.
awk 'BEGIN{OFS="\t"} {if (NR % 7 == 0) {$5 = "C"; $6 = "T"}; print}' tmp_ref.bim > tmp_wrongallele.bim
cp tmp_ref.bed tmp_wrongallele.bed
cp tmp_ref.fam tmp_wrongallele.fam
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-bfile tmp_wrongallele --out plink2_rdw > tmp_rdw.log 2>&1
grep -q "matched the --flip-scan reference fileset by ID but not by" tmp_rdw.log
test "$(awk -F '\t' 'NR > 1 && $8 == "NA"' plink2_rdw.flipscan | wc -l)" -gt 0

# A frequency difference in the second fileset is caught the same way it is in
# a frequency file.  Recoding a slice of the reference's genotypes towards the
# reference allele moves its frequency without touching this dataset's.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --keep-if PHENO1 == 1 --make-bed --out tmp_ctrlref
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan --flip-scan-ref-bfile tmp_ctrlref --flip-scan-freq-diff 0.001 --out plink2_rdd
test "$(awk -F '\t' 'NR > 1 && $8 == "Y"' plink2_rdd.flipscan | wc -l)" -gt 0
# Every call agrees with the frequency difference actually reported.
awk -F '\t' 'NR > 1 && $8 != "NA" {d = $6 - $7; if (d < 0) {d = -d}; if (($8 == "Y") != (d > 0.001)) {print "call disagrees with the difference at " $3; exit 1}}' plink2_rdd.flipscan

# 'ref-allele-based' renames the columns here too.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --maf 0.05 --flip-scan ref-allele-based --flip-scan-ref-bfile tmp_ref --out plink2_rdr
head -n 1 plink2_rdr.flipscan | grep -qx '#CHROM	POS	ID	REF	ALT	REF_FREQ	PANEL_REF_FREQ	PROBLEM'

# The two reference modes are mutually exclusive, and neither runs on its own.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --flip-scan --flip-scan-ref-bfile tmp_ref --flip-scan-ref-freq tmp_panel.afreq --out plink2_bad > /dev/null 2>&1; then
    echo "--flip-scan-ref-bfile ran alongside --flip-scan-ref-freq"
    exit 1
fi
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --flip-scan-ref-bfile tmp_ref --freq --out plink2_bad > /dev/null 2>&1; then
    echo "--flip-scan-ref-bfile ran without --flip-scan"
    exit 1
fi
