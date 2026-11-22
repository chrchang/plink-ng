#!/bin/bash

# PLINK 1.9 sample and variant filters, checked against the equivalent plink2
# filters.
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

# 213 samples in families of 5 (two founders, then three children), with
# every 11th sample's sex unknown and every 13th phenotype missing.  300
# variants, 10 kb apart: 150 on chr1, 80 on chr2, 40 on chrX, 20 on XY and 10
# on MT.  Every 17th variant has a two-base A1 allele, and every 23rd an A2
# allele of N.  Heterozygous haploid and MT calls are set to missing, since
# plink2 has no such state.
plink --dummy 213 300 0.02 acgt --seed 41 --make-bed --out tmp_raw > /dev/null
awk '{
    base = NR - (NR - 1) % 5
    pat = 0
    mat = 0
    if ((NR - 1) % 5 >= 2) {
        pat = "i" (base + 1)
        mat = "i" base
    }
    print "f" int((NR - 1) / 5), "i" NR, pat, mat, (NR % 11 == 0)? 0 : 1 + NR % 2, (NR % 13 == 0)? -9 : $6
}' tmp_raw.fam > tmp_raw.fam.new
mv tmp_raw.fam.new tmp_raw.fam
awk 'BEGIN { OFS = "\t" } {
    $1 = (NR <= 150)? 1 : (NR <= 230)? 2 : (NR <= 270)? 23 : (NR <= 290)? 25 : 26
    $4 = 10000 * NR
    if (NR % 17 == 0) $5 = $5 "T"
    if (NR % 23 == 0) $6 = "N"
    print
}' tmp_raw.bim > tmp_raw.bim.new
mv tmp_raw.bim.new tmp_raw.bim
plink --bfile tmp_raw --set-hh-missing --keep-allele-order --make-bed --out tmp_data > /dev/null 2>&1

# --- Sample filters ---

# $1: label, $2: PLINK 1.9 flags, $3: plink2 flags.  The remaining samples'
# IDs must match.
same_samples() {
    plink --bfile tmp_data $2 --make-just-fam --out tmp_s19 > /dev/null 2>&1
    $plink2 --bfile tmp_data $3 --make-just-psam --out tmp_s2 > /dev/null 2>&1
    test -s tmp_s19.fam && cmp -s <(cut -d ' ' -f 1,2 tmp_s19.fam) <(tail -n +2 tmp_s2.psam | cut -f 1,2 | tr '\t' ' ') && \
        test "$(wc -l < tmp_s19.fam)" -lt 213
    check "$1" $?
    rm -f tmp_s19.fam tmp_s2.psam
}
same_samples "--filter-cases matches plink2 --keep-if PHENO1 == case" "--filter-cases" "--keep-if PHENO1 == case"
same_samples "--filter-males matches plink2 --keep-males" "--filter-males" "--keep-males"
same_samples "--filter-nonfounders matches plink2 --keep-nonfounders" "--filter-nonfounders" "--keep-nonfounders"
same_samples "--prune matches plink2 --require-pheno" "--prune" "--require-pheno"
printf 'f3\nf7\nf40\n' > tmp_fams.txt
same_samples "--keep-fam matches plink2" "--keep-fam tmp_fams.txt" "--keep-fam tmp_fams.txt"
same_samples "--remove-fam matches plink2" "--remove-fam tmp_fams.txt" "--remove-fam tmp_fams.txt"

# --filter/--mfilter and the cluster filters.  plink2 reads the same file as
# categorical phenotypes C1 and C2.
awk '{ print $1, $2, "c" (NR % 4), "g" (NR % 3) }' tmp_data.fam > tmp_cat.txt
awk 'BEGIN { print "FID IID C1 C2" } { print }' tmp_cat.txt > tmp_cat_h.txt
printf 'c1\nc3\n' > tmp_cl.txt
same_samples "--filter matches plink2 --keep-if" "--filter tmp_cat.txt c2" "--pheno tmp_cat_h.txt --keep-if C1 == c2"
same_samples "--filter --mfilter 2 matches plink2 --keep-if" "--filter tmp_cat.txt g1 --mfilter 2" "--pheno tmp_cat_h.txt --keep-if C2 == g1"
same_samples "--keep-clusters matches plink2 --keep-cats" "--within tmp_cat.txt --keep-clusters tmp_cl.txt" "--pheno tmp_cat_h.txt --keep-cat-pheno C1 --keep-cats tmp_cl.txt"
same_samples "--remove-clusters matches plink2 --remove-cats" "--within tmp_cat.txt --remove-clusters tmp_cl.txt" "--pheno tmp_cat_h.txt --remove-cat-pheno C1 --remove-cats tmp_cl.txt"
same_samples "--mwithin 2 --keep-cluster-names matches plink2 --keep-cat-names" "--within tmp_cat.txt --mwithin 2 --keep-cluster-names g0 g2" "--pheno tmp_cat_h.txt --keep-cat-pheno C2 --keep-cat-names g0 g2"
same_samples "--mwithin 2 --remove-cluster-names matches plink2 --remove-cat-names" "--within tmp_cat.txt --mwithin 2 --remove-cluster-names g0" "--pheno tmp_cat_h.txt --remove-cat-pheno C2 --remove-cat-names g0"

# --thin-indiv-count draws randomly, so only the count can be checked.
plink --bfile tmp_data --thin-indiv-count 50 --seed 3 --make-just-fam --out tmp_ti > /dev/null 2>&1
test "$(wc -l < tmp_ti.fam)" -eq 50
check "--thin-indiv-count keeps the given number of samples" $?

# --must-have-sex has no plink2 equivalent: it sets the phenotypes of
# unknown-sex samples to missing in the output.
plink --bfile tmp_data --must-have-sex --make-just-fam --out tmp_mhs > /dev/null 2>&1
awk '{ if ($5 == 0) $6 = -9; print }' tmp_data.fam | cmp -s - tmp_mhs.fam && ! cmp -s tmp_data.fam tmp_mhs.fam
check "--must-have-sex sets unknown-sex phenotypes to missing" $?

# --- Variant filters ---

# $1: label, $2: PLINK 1.9 flags, $3: plink2 flags.  The remaining variants'
# IDs must match, and the filter must remove something.
same_variants() {
    plink --bfile tmp_data $2 --write-snplist --out tmp_v19 > /dev/null 2>&1
    $plink2 --bfile tmp_data $3 --write-snplist --out tmp_v2 > /dev/null 2>&1
    test -s tmp_v19.snplist && cmp -s tmp_v19.snplist tmp_v2.snplist && test "$(wc -l < tmp_v19.snplist)" -lt 300
    check "$1" $?
    rm -f tmp_v19.snplist tmp_v2.snplist
}
same_variants "--not-chr matches plink2" "--not-chr 2,X" "--not-chr 2,X"
same_variants "--autosome-xy matches plink2 --autosome-par" "--autosome-xy" "--autosome-par"
same_variants "--snp --window matches plink2" "--snp snp50 --window 200" "--snp snp50 --window 200"
same_variants "--snps matches plink2" "--snps snp10-snp20,snp40,snp200-snp205" "--snps snp10-snp20,snp40,snp200-snp205"
same_variants "--exclude-snps matches plink2" "--exclude-snps snp10-snp20,snp40" "--exclude-snps snp10-snp20,snp40"
same_variants "--d with --snps matches plink2" "--d : --snps snp10:snp20,snp40" "--d : --snps snp10:snp20,snp40"
same_variants "--from-bp/--to-bp matches plink2" "--chr 1 --from-bp 50000 --to-bp 200000" "--chr 1 --from-bp 50000 --to-bp 200000"
same_variants "--from-kb/--to-kb matches plink2" "--chr 2 --from-kb 1600 --to-kb 2000" "--chr 2 --from-kb 1600 --to-kb 2000"
same_variants "--from-mb/--to-mb matches plink2" "--chr 1 --from-mb 0.3 --to-mb 1.15" "--chr 1 --from-mb 0.3 --to-mb 1.15"
same_variants "--snps-only matches plink2" "--snps-only" "--snps-only"
same_variants "--snps-only just-acgt matches plink2" "--snps-only just-acgt" "--snps-only just-acgt"
# Minor allele counts are over founders.  MT is left out: PLINK 1.9 counts two
# alleles per MT call there, and plink2 one.
same_variants "--max-mac matches plink2 --max-mac --ac-founders" "--not-chr MT --max-mac 70" "--not-chr MT --max-mac 70 --ac-founders"

plink --bfile tmp_data --thin-count 40 --seed 3 --write-snplist --out tmp_tc > /dev/null 2>&1
test "$(wc -l < tmp_tc.snplist)" -eq 40
check "--thin-count keeps the given number of variants" $?

# --maf-succ: frequencies must be (count + 1) / (total + 2), from --freq
# counts.  plink2 --af-pseudocount 1 is the same estimator, but it doesn't
# change plink2's --freq report, so it is compared through --maf: across
# these thresholds, --maf-succ must matter at least once and always match.
plink --bfile tmp_data --chr 1-2 --keep-allele-order --freq counts --out tmp_fc > /dev/null 2>&1
plink --bfile tmp_data --chr 1-2 --keep-allele-order --maf-succ --freq --out tmp_fs > /dev/null 2>&1
paste <(tail -n +2 tmp_fc.frq.counts | tr -s ' ' | sed 's/^ //') <(tail -n +2 tmp_fs.frq | tr -s ' ' | sed 's/^ //') | \
    awk '{ e = ($5 + 1) / ($5 + $6 + 2); d = e - $12; if (d < 0) d = -d; if ($2 != $9 || d > 6e-5) bad = 1 } END { exit (bad || NR != 230) }'
succ_ok=$?
succ_mattered=0
for t in 0.43 0.44 0.45 0.46 0.47; do
    plink --bfile tmp_data --chr 1-2 --maf-succ --maf $t --write-snplist --out tmp_ms19 > /dev/null 2>&1
    plink --bfile tmp_data --chr 1-2 --maf $t --write-snplist --out tmp_mp19 > /dev/null 2>&1
    $plink2 --bfile tmp_data --chr 1-2 --af-pseudocount 1 --maf $t --write-snplist --out tmp_ms2 > /dev/null 2>&1
    cmp -s tmp_ms19.snplist tmp_ms2.snplist || succ_ok=1
    cmp -s tmp_ms19.snplist tmp_mp19.snplist || succ_mattered=1
done
test $succ_ok -eq 0 && test $succ_mattered -eq 1
check "--maf-succ matches (count + 1) / (total + 2) and plink2 --af-pseudocount 1" $?

# --fill-missing-a2 with --keep-allele-order: A2 is what plink2 reads as REF.
plink --bfile tmp_data --keep-allele-order --fill-missing-a2 --make-bed --out tmp_fm19 > /dev/null 2>&1
$plink2 --bfile tmp_data --fill-missing-with-ref --output-chr 26 --make-bed --out tmp_fm2 > /dev/null 2>&1
cmp -s tmp_fm19.bed tmp_fm2.bed && cmp -s tmp_fm19.bim tmp_fm2.bim && ! cmp -s tmp_fm19.bed tmp_data.bed
check "--fill-missing-a2 matches plink2 --fill-missing-with-ref" $?

# --- VCF and quality-score filters ---

# 40 sites, 9 samples.  FILTER cycles through PASS, '.', q10, lowDP;q10 and
# lowDP, and every 7th site (from the 3rd) has only a DS field.
awk 'BEGIN {
    OFS = "\t"
    s = 3
    print "##fileformat=VCFv4.2"
    print "##contig=<ID=1,length=100000000>"
    print "##FILTER=<ID=q10,Description=\"q\">"
    print "##FILTER=<ID=lowDP,Description=\"d\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">"
    print "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"d\">"
    h = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (i = 1; i <= 9; ++i) h = h "\ts" i
    print h
    split("PASS . q10 lowDP;q10 lowDP", filt, " ")
    for (v = 1; v <= 40; ++v) {
        fmt = (v % 7 == 3)? "DS" : "GT"
        line = "1\t" (1000 * v) "\tv" v "\tA\tG\t" (v * 3) "\t" filt[1 + v % 5] "\t.\t" fmt
        for (i = 1; i <= 9; ++i) {
            s = (s * 69069 + 1) % 4294967296
            r = int(s / 65536) % 3
            line = line "\t" ((fmt == "DS")? r : ((r == 0)? "0/0" : (r == 1)? "0/1" : "1/1"))
        }
        print line
    }
}' > tmp_filt.vcf
# The quality scores are the QUAL values, with a header line and a dummy
# column, so the column arguments are exercised.
awk 'BEGIN { print "ID X QUAL" } !/^#/ { print $3, "x", $6 }' tmp_filt.vcf > tmp_qual.txt

# $1: label, $2: PLINK 1.9 flags, $3: plink2 flags.  plink2 converts the VCF
# with --make-pgen: if nothing needs genotypes, it reads the VCF as a variant
# file and ignores --vcf-require-gt.
same_vcf_variants() {
    plink --vcf tmp_filt.vcf $2 --write-snplist --out tmp_vcf19 > /dev/null 2>&1
    $plink2 --vcf tmp_filt.vcf $3 --make-pgen --out tmp_vcf2 > /dev/null 2>&1
    test -s tmp_vcf19.snplist && cmp -s tmp_vcf19.snplist <(grep -v '^#' tmp_vcf2.pvar | cut -f 3) && test "$(wc -l < tmp_vcf19.snplist)" -lt 40
    check "$1" $?
    rm -f tmp_vcf19.snplist tmp_vcf2.pvar
}
same_vcf_variants "--vcf-filter matches plink2 --var-filter" "--vcf-filter" "--var-filter"
same_vcf_variants "--vcf-filter with an exception matches plink2 --var-filter" "--vcf-filter q10" "--var-filter q10"
same_vcf_variants "--vcf-require-gt matches plink2" "--vcf-require-gt" "--vcf-require-gt"
same_vcf_variants "--vcf-filter --vcf-require-gt matches plink2" "--vcf-filter lowDP --vcf-require-gt" "--var-filter lowDP --vcf-require-gt"
same_vcf_variants "--qual-scores --qual-threshold --qual-max-threshold matches plink2 --extract-col-cond" \
    "--qual-scores tmp_qual.txt 3 1 1 --qual-threshold 20 --qual-max-threshold 80" \
    "--extract-col-cond tmp_qual.txt 3 1 1 --extract-col-cond-min 20 --extract-col-cond-max 80"
same_vcf_variants "--qual-scores --qual-threshold matches plink2 --extract-col-cond-min" \
    "--qual-scores tmp_qual.txt 3 1 1 --qual-threshold 50" \
    "--extract-col-cond tmp_qual.txt 3 1 1 --extract-col-cond-min 50"

exit $failed
