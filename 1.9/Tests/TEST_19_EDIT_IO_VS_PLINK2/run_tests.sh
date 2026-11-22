#!/bin/bash

# PLINK 1.9 data-editing flags and input/output format options, checked
# against the equivalent plink2 flags.
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

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
# Text with tabs turned into spaces.  (plink2 writes tab-delimited .fam/.ped
# files.)
spaces() {
    tr '\t' ' ' < "$1"
}
# $1, $2: filesets that must be identical (.bed, .bim, and .fam up to
# delimiters).
same_fileset() {
    test -s "$1.bed" && cmp -s "$1.bed" "$2.bed" && cmp -s "$1.bim" "$2.bim" && cmp -s <(spaces "$1.fam") <(spaces "$2.fam")
}

# --- Data editing ---

# 213 samples in families of 5 (two founders, then three children), with
# every 11th sample's sex unknown.  300 variants, 10 kb and 0.01 cM apart: 140
# on chr1, 80 on chr2, 40 on chrX, 10 on chrY, 10 on XY and 20 on MT.
# tmp_raw keeps the heterozygous haploid and MT calls; tmp_data has them set
# to missing, since plink2 has no such state.
plink --dummy 213 300 0.02 acgt --seed 43 --make-bed --out tmp_raw > /dev/null
awk '{
    base = NR - (NR - 1) % 5
    pat = 0
    mat = 0
    if ((NR - 1) % 5 >= 2) {
        pat = "i" (base + 1)
        mat = "i" base
    }
    print "f" int((NR - 1) / 5), "i" NR, pat, mat, (NR % 11 == 0)? 0 : 1 + NR % 2, $6
}' tmp_raw.fam > tmp_raw.fam.new
mv tmp_raw.fam.new tmp_raw.fam
awk 'BEGIN { OFS = "\t" } {
    $1 = (NR <= 140)? 1 : (NR <= 220)? 2 : (NR <= 260)? 23 : (NR <= 270)? 24 : (NR <= 280)? 25 : 26
    $3 = NR / 100
    $4 = 10000 * NR
    print
}' tmp_raw.bim > tmp_raw.bim.new
mv tmp_raw.bim.new tmp_raw.bim
plink --bfile tmp_raw --set-hh-missing --keep-allele-order --make-bed --out tmp_data > /dev/null 2>&1

# $1: label, $2: input fileset, $3: PLINK 1.9 flags, $4: plink2 flags, $5:
# nonempty if plink2 must sort the variants (which it can only do while
# writing a .pgen).  The output filesets must be identical, and must differ
# from the input.
same_edit() {
    plink --bfile "$2" --keep-allele-order $3 --make-bed --out tmp_e19 > /dev/null 2>&1
    if [[ -n "${5:-}" ]]; then
        $plink2 --bfile "$2" $4 --make-pgen --sort-vars --out tmp_e2p > /dev/null 2>&1
        $plink2 --pfile tmp_e2p --output-chr 26 --make-bed --out tmp_e2 > /dev/null 2>&1
    else
        $plink2 --bfile "$2" $4 --output-chr 26 --make-bed --out tmp_e2 > /dev/null 2>&1
    fi
    same_fileset tmp_e19 tmp_e2 && ! same_fileset tmp_e19 "$2"
    check "$1" $?
    rm -f tmp_e19.* tmp_e2.* tmp_e2p.*
}
awk '$1 == 1 && NR % 9 == 0 { print $2, 2 }' tmp_data.bim > tmp_uchr.txt
awk 'NR % 4 == 0 { print $2, NR / 50 }' tmp_data.bim > tmp_ucm.txt
awk '$1 <= 2 && NR % 5 == 0 { print $2, 10000 * NR + 15003 }' tmp_data.bim > tmp_umap.txt
awk 'NR % 6 == 0 { print $2, "rs" NR }' tmp_data.bim > tmp_uname.txt
awk 'NR % 7 == 0 { print $1, $2, "i1", "i2" }' tmp_data.fam > tmp_upar.txt
awk 'NR % 8 == 0 { print $1, $2, ($5 == 1)? 2 : 1 }' tmp_data.fam > tmp_usex.txt
same_edit "--update-chr matches plink2" tmp_data "--update-chr tmp_uchr.txt" "--update-chr tmp_uchr.txt" sort
same_edit "--update-cm matches plink2" tmp_data "--update-cm tmp_ucm.txt" "--update-cm tmp_ucm.txt"
same_edit "--update-map matches plink2" tmp_data "--update-map tmp_umap.txt" "--update-map tmp_umap.txt" sort
same_edit "--update-name matches plink2" tmp_data "--update-name tmp_uname.txt" "--update-name tmp_uname.txt"
same_edit "--update-parents matches plink2" tmp_data "--update-parents tmp_upar.txt" "--update-parents tmp_upar.txt"
same_edit "--update-sex matches plink2" tmp_data "--update-sex tmp_usex.txt" "--update-sex tmp_usex.txt"
same_edit "--zero-cms matches plink2" tmp_data "--zero-cms" "--zero-cms"
same_edit "--set-mixed-mt-missing matches plink2" tmp_raw "--set-mixed-mt-missing" "--set-mixed-mt-missing"
awk '{ print $1, $2, "k" (NR % 3) }' tmp_data.fam > tmp_zc_within.txt
awk 'NR % 5 == 0 { print $2, "k1" }' tmp_data.bim > tmp_zc.txt
same_edit "--zero-cluster matches plink2" tmp_data "--within tmp_zc_within.txt --zero-cluster tmp_zc.txt" "--within tmp_zc_within.txt --zero-cluster tmp_zc.txt"

# --real-ref-alleles keeps the .bim allele order, as plink2 always does.
# (Without it, PLINK 1.9 makes A1 the minor allele.)
plink --bfile tmp_data --real-ref-alleles --make-bed --out tmp_rr19 > /dev/null 2>&1
plink --bfile tmp_data --make-bed --out tmp_nr19 > /dev/null 2>&1
$plink2 --bfile tmp_data --output-chr 26 --make-bed --out tmp_rr2 > /dev/null 2>&1
same_fileset tmp_rr19 tmp_rr2 && ! cmp -s tmp_nr19.bim tmp_rr19.bim
check "--real-ref-alleles matches plink2" $?

# --set-missing-var-ids ('$1'/'$2' are the alleles in alphabetical order in
# both), with --new-id-max-allele-len (PLINK 1.9 truncates; plink2 does with
# 'truncate'; every 17th variant has a two-base allele), and
# --missing-var-code.
awk 'BEGIN { OFS = "\t" } { if (NR % 4 == 0) $2 = "."; if (NR % 17 == 0) $5 = $5 "T"; print }' tmp_data.bim > tmp_mid.bim
cp tmp_data.bed tmp_mid.bed
cp tmp_data.fam tmp_mid.fam
same_edit "--set-missing-var-ids matches plink2" tmp_mid "--set-missing-var-ids @:#\$1,\$2" "--set-missing-var-ids @:#\$1,\$2"
same_edit "--set-missing-var-ids --new-id-max-allele-len matches plink2" tmp_mid "--set-missing-var-ids @:#\$1,\$2 --new-id-max-allele-len 1" "--set-missing-var-ids @:#\$1,\$2 --new-id-max-allele-len 1 truncate"
awk 'BEGIN { OFS = "\t" } { if (NR % 4 == 0) $2 = "NA"; print }' tmp_data.bim > tmp_mvc.bim
cp tmp_data.bed tmp_mvc.bed
cp tmp_data.fam tmp_mvc.fam
same_edit "--missing-var-code matches plink2" tmp_mvc "--missing-var-code NA --set-missing-var-ids @_#" "--missing-var-code NA --set-missing-var-ids @_#"

# --split-x and --merge-x, on a copy without XY (which --split-x requires),
# with chrX positions inside and between the b37 pseudoautosomal regions.
# plink2 --split-par places PAR1 before X and PAR2 after it, while PLINK 1.9
# puts all of XY after X, so the split is compared variant by variant.
plink --bfile tmp_data --not-chr XY --keep-allele-order --make-bed --out tmp_xs > /dev/null 2>&1
awk 'BEGIN { OFS = "\t" } $1 == 23 {
    k = ++xk
    $4 = (k <= 10)? 100000 * k : (k <= 30)? 3000000 + 1000000 * k : 154931043 + 10000 * (k - 30)
} { print }' tmp_xs.bim > tmp_xs.bim.new
mv tmp_xs.bim.new tmp_xs.bim
awk '{ print $2, $5 }' tmp_xs.bim > tmp_xs_a1.txt
# Prints "ID IID dosage" triples from a --recode A/--export A file.
raw_triples() {
    awk 'NR == 1 { for (i = 7; i <= NF; ++i) { n = $i; sub(/_[^_]*$/, "", n); id[i] = n } next }
      { for (i = 7; i <= NF; ++i) print id[i], $2, $i }' "$1" | LC_ALL=C sort
}
plink --bfile tmp_xs --keep-allele-order --split-x b37 --make-bed --out tmp_sx19 > /dev/null 2>&1
plink --bfile tmp_sx19 --recode A --recode-allele tmp_xs_a1.txt --out tmp_sx19 > /dev/null 2>&1
$plink2 --bfile tmp_xs --split-par b37 --make-pgen --sort-vars --out tmp_sx2 > /dev/null 2>&1
$plink2 --pfile tmp_sx2 --export A --export-allele tmp_xs_a1.txt --out tmp_sx2 > /dev/null 2>&1
cmp -s <(raw_triples tmp_sx19.raw) <(raw_triples tmp_sx2.raw) && \
    cmp -s <(awk '{ print $2, $1, $4 }' tmp_sx19.bim | LC_ALL=C sort) \
      <(awk '!/^#/ { c = $1; if (c == "X") c = 23; if (c == "PAR1" || c == "PAR2") c = 25; if (c == "Y") c = 24; if (c == "MT") c = 26; print $3, c, $2 }' tmp_sx2.pvar | LC_ALL=C sort) && \
    test "$(awk '$1 == 25' tmp_sx19.bim | wc -l)" -eq 20
check "--split-x matches plink2 --split-par" $?
same_edit "--merge-x matches plink2" tmp_sx19 "--merge-x" "--merge-x" sort

# --- Input and output formats ---

# 37 samples, 40 variants on chr1, chrX, chrY, XY and MT, 5% missing, as a
# text fileset.
plink --dummy 37 40 0.05 acgt --seed 5 --make-bed --out tmp_p > /dev/null
awk 'BEGIN { OFS = "\t" } { $1 = (NR <= 20)? 1 : (NR <= 28)? 23 : (NR <= 32)? 24 : (NR <= 36)? 25 : 26; $4 = 1000 * NR; print }' tmp_p.bim > tmp_p.bim.new
mv tmp_p.bim.new tmp_p.bim
plink --bfile tmp_p --keep-allele-order --recode --out tmp_p > /dev/null 2>&1

# $1: label, $2: PLINK 1.9 input flags, $3: plink2 input flags, $4: extra flags
# for both, $5: optional extended regular expression that PLINK 1.9's .ped
# must match somewhere, so that the flag can be seen to have an effect.  The
# exported .ped/.map files must match.  (plink2 writes chromosome names unless
# told otherwise.)
same_ped() {
    local chr_flag="--output-chr 26"
    if [[ "${4:-}" == *--output-chr* ]]; then
        chr_flag=
    fi
    plink $2 ${4:-} --keep-allele-order --recode --out tmp_r19 > /dev/null 2>&1
    $plink2 $3 ${4:-} $chr_flag --export ped --out tmp_r2 > /dev/null 2>&1
    test -s tmp_r19.ped && cmp -s tmp_r19.ped <(spaces tmp_r2.ped) && cmp -s tmp_r19.map tmp_r2.map && \
        { [[ -z "${5:-}" ]] || grep -Eq -- "$5" tmp_r19.ped; }
    check "$1" $?
    rm -f tmp_r19.ped tmp_r19.map tmp_r2.ped tmp_r2.map
}
cp tmp_p.ped tmp_named1.txt
cp tmp_p.map tmp_named2.txt
same_ped "--ped/--map matches plink2" "--ped tmp_named1.txt --map tmp_named2.txt" "--ped tmp_named1.txt --map tmp_named2.txt"
plink --file tmp_p --keep-allele-order --recode transpose --out tmp_t > /dev/null 2>&1
cp tmp_t.tped tmp_named3.txt
cp tmp_t.tfam tmp_named4.txt
same_ped "--tped/--tfam matches plink2" "--tped tmp_named3.txt --tfam tmp_named4.txt" "--tped tmp_named3.txt --tfam tmp_named4.txt"
plink --file tmp_p --keep-allele-order --recode lgen --out tmp_lg > /dev/null 2>&1
same_ped "--lgen matches plink2" "--lgen tmp_lg.lgen --map tmp_lg.map --fam tmp_lg.fam" "--lgen tmp_lg.lgen --map tmp_lg.map --fam tmp_lg.fam"
plink --file tmp_p --keep-allele-order --recode lgen-ref --out tmp_lr > /dev/null 2>&1
same_ped "--lgen --reference matches plink2" "--lgen tmp_lr.lgen --map tmp_lr.map --fam tmp_lr.fam --reference tmp_lr.ref" "--lgen tmp_lr.lgen --map tmp_lr.map --fam tmp_lr.fam --reference tmp_lr.ref"
# The same genotypes as reference-allele counts.  (Missing calls are left
# out, so they become homozygous reference in both.)
awk 'NR == FNR { ref[$1] = $2; next } $4 != "0" { print $1, $2, $3, ($4 == ref[$3]) + ($5 == ref[$3]) }' tmp_lr.ref tmp_lr.lgen > tmp_ac.lgen
same_ped "--lgen --reference --allele-count matches plink2" "--lgen tmp_ac.lgen --map tmp_lr.map --fam tmp_lr.fam --reference tmp_lr.ref --allele-count" "--lgen tmp_ac.lgen --map tmp_lr.map --fam tmp_lr.fam --reference tmp_lr.ref --allele-count"

# --no-parents/--no-sex/--no-pheno, and --no-fid.  With --no-fid, PLINK 1.9
# sets each FID to the IID, while plink2 sets it to 0, so FIDs are left out.
awk '{ $3 = $4 = $5 = $6 = ""; gsub(/  +/, " "); print }' tmp_p.ped > tmp_nps.ped
cp tmp_p.map tmp_nps.map
same_ped "--no-parents --no-sex --no-pheno matches plink2" "--file tmp_nps" "--pedmap tmp_nps" "--no-parents --no-sex --no-pheno"
awk '{ $1 = ""; sub(/^ /, ""); print }' tmp_p.ped > tmp_nofid.ped
cp tmp_p.map tmp_nofid.map
plink --file tmp_nofid --no-fid --keep-allele-order --recode --out tmp_nf19 > /dev/null 2>&1
$plink2 --pedmap tmp_nofid --no-fid --export ped --out tmp_nf2 > /dev/null 2>&1
cmp -s <(cut -d ' ' -f 2- tmp_nf19.ped) <(spaces tmp_nf2.ped | cut -d ' ' -f 2-) && test "$(cut -d ' ' -f 1 tmp_nf19.ped | head -1)" = per0
check "--no-fid matches plink2 (apart from the FID column)" $?

# Missing-genotype and missing-phenotype codes.
awk '{ for (i = 7; i <= NF; ++i) if ($i == 0) $i = "N"; print }' tmp_p.ped > tmp_mg.ped
cp tmp_p.map tmp_mg.map
same_ped "--missing-genotype/--output-missing-genotype matches plink2" "--file tmp_mg --missing-genotype N" "--pedmap tmp_mg --input-missing-genotype N" "--output-missing-genotype X" " X X( |$)"
awk '{ if (NR % 5 == 0) $6 = -99; print }' tmp_p.ped > tmp_mp.ped
cp tmp_p.map tmp_mp.map
same_ped "--missing-phenotype/--output-missing-phenotype matches plink2" "--file tmp_mp --missing-phenotype -99" "--pedmap tmp_mp --input-missing-phenotype -99" "--output-missing-phenotype -7" "^[^ ]+ [^ ]+ [^ ]+ [^ ]+ [^ ]+ -7 "
awk '{ $6 = (NR % 7 == 0)? -9 : (NR % 3 == 0)? 1 : 0; print }' tmp_p.ped > tmp_one.ped
cp tmp_p.map tmp_one.map
same_ped "--1 matches plink2" "--file tmp_one" "--pedmap tmp_one" "--1" "^[^ ]+ [^ ]+ [^ ]+ [^ ]+ [^ ]+ 2 "

# Allele codes.  (plink2 --alleleACGT requires 1234-coded input.)
same_ped "--allele1234 matches plink2" "--file tmp_p" "--pedmap tmp_p" "--allele1234" "^([^ ]+ ){6}[1-4] "
plink --file tmp_p --allele1234 --keep-allele-order --recode --out tmp_n12 > /dev/null 2>&1
same_ped "--alleleACGT matches plink2" "--file tmp_n12" "--pedmap tmp_n12" "--alleleACGT" "^([^ ]+ ){6}[ACGT] "

for code in M MT 26 chrM chrMT chr26; do
    same_ped "--output-chr $code matches plink2" "--file tmp_p" "--pedmap tmp_p" "--output-chr $code"
done

# VCF sample IDs: "fam<n>:ind <n>" with --id-delim and --vcf-idspace-to, and
# single IDs with --const-fid.  (chr1 only: plink2 won't import chrX from a
# VCF without sex information.)
plink --file tmp_p --chr 1 --keep-allele-order --recode vcf --out tmp_v > /dev/null 2>&1
sed 's/per\([0-9]*\)_per[0-9]*/fam\1:ind \1/g' tmp_v.vcf > tmp_vs.vcf
same_ped "--id-delim --vcf-idspace-to matches plink2 --id-delim --idspace-to" "--vcf tmp_vs.vcf --id-delim : --vcf-idspace-to _" "--vcf tmp_vs.vcf --id-delim : --idspace-to _" "" "^fam0 ind_0 "
sed 's/per\([0-9]*\)_per[0-9]*/s\1/g' tmp_v.vcf > tmp_vc.vcf
same_ped "--const-fid matches plink2" "--vcf tmp_vc.vcf --const-fid F1" "--vcf tmp_vc.vcf --const-fid F1" "" "^F1 s0 "

# Oxford format.  The probabilities put all weight on two adjacent
# genotypes: PLINK 1.9's --hard-call-threshold compares 1 minus the largest
# probability with the threshold, and plink2's compares the dosage's
# distance from the nearest integer, and these agree only then.  The .gen
# file has no chromosome codes (--oxford-single-chr), and the .sample file has
# an extra phenotype with NA missing codes.
plink --bfile tmp_p --chr 1 --keep-allele-order --recode oxford --out tmp_ox > /dev/null 2>&1
awk 'BEGIN { s = 9 } {
    printf "--- %s %s %s %s", $2, $3, $4, $5
    for (i = 6; i <= NF; i += 3) {
        s = (s * 69069 + 1) % 4294967296
        t = int(s / 65536) % 1000 / 1000 * 0.4
        if ($i + $(i + 1) + $(i + 2) == 0) printf " 0 0 0"
        else if ($i == 1) printf " %.3f %.3f 0", 1 - t, t
        else if ($(i + 1) == 1) printf " 0 %.3f %.3f", 1 - t, t
        else printf " 0 %.3f %.3f", t, 1 - t
    }
    printf "\n"
}' tmp_ox.gen > tmp_oxp.gen
awk 'NR == 1 { print $0, "qt"; next } NR == 2 { print $0, "P"; next } { print $0, (NR % 4 == 0)? "NA" : (NR % 7 == 0)? "-9" : NR / 10 }' tmp_ox.sample > tmp_oxp.sample
for thr in 0.1 0.25; do
    same_ped "--oxford-single-chr --hard-call-threshold $thr matches plink2" "--gen tmp_oxp.gen --sample tmp_oxp.sample" "--gen tmp_oxp.gen ref-last --sample tmp_oxp.sample" "--oxford-single-chr 5 --hard-call-threshold $thr"
    plink --gen tmp_oxp.gen --sample tmp_oxp.sample --oxford-single-chr 5 --hard-call-threshold $thr --keep-allele-order --recode --out tmp_hc$thr > /dev/null 2>&1
done
! cmp -s tmp_hc0.1.ped tmp_hc0.25.ped
check "--hard-call-threshold changes the calls" $?
same_ped "--oxford-pheno-name --missing-code matches plink2 --pheno-name --missing-code" "--gen tmp_oxp.gen --sample tmp_oxp.sample --oxford-pheno-name qt" "--gen tmp_oxp.gen ref-last --sample tmp_oxp.sample --pheno-name qt" "--oxford-single-chr 5 --missing-code NA" "^per0 per0 0 0 [0-9] 0\\.3 "

# --keep-autoconv has no plink2 counterpart in .bed form: the files it keeps
# must match a --make-bed conversion.
rm -f tmp_ka.bed tmp_ka.bim tmp_ka.fam
plink --vcf tmp_v.vcf --keep-autoconv --freq --out tmp_ka > /dev/null 2>&1
plink --vcf tmp_v.vcf --make-bed --out tmp_kb > /dev/null 2>&1
same_fileset tmp_ka tmp_kb
check "--keep-autoconv keeps the converted fileset" $?

exit $failed
