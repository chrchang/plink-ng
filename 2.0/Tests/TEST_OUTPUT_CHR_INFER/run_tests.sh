#!/bin/bash

set -exo pipefail

# Writes a small VCF whose chromosome codes are the arguments, in order.  Each
# code gets a ##contig line and one variant; positions are kept outside the
# human PARs.
make_vcf() {
    out=$1
    shift
    {
        printf '##fileformat=VCFv4.2\n'
        for c in "$@"; do
            printf '##contig=<ID=%s,length=100000000>\n' $c
        done
        printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n'
        pos=10000000
        for c in "$@"; do
            printf '%s\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\t0/0\n' $c $pos $pos
            pos=$((pos + 1000))
        done
    } > $out
}

# Prints the ##contig IDs, then the CHROM column.
vcf_chrs() {
    grep '^##contig=<ID=' $1 | sed 's/^##contig=<ID=\([^,>]*\).*/\1/'
    grep -v '^#' $1 | cut -f 1
}

pvar_chrs() {
    grep -v '^#' $1 | cut -f 1
}

check_vcf() {
    vcf_chrs $1 > tmp_expected.txt
    vcf_chrs $2 > tmp_observed.txt
    diff -q tmp_expected.txt tmp_observed.txt
}

# Sex information is needed to import chrX.
printf '#IID\tSEX\ns1\t1\ns2\t2\ns3\t2\n' > tmp_samples.psam

# UCSC-style names, plus an extra contig which must pass through unchanged.
make_vcf tmp_ucsc.vcf chr1 chr2 chr22 chrX chrY chrM chrUn_gl000220
# Ensembl-style names.
make_vcf tmp_ensembl.vcf 1 2 22 X Y MT
make_vcf tmp_plainm.vcf 1 2 X M
make_vcf tmp_numeric.vcf 1 2 23 24 26
make_vcf tmp_ucscmt.vcf chr1 chrX chrMT
make_vcf tmp_autosomes.vcf chr1 chr2

# Without infer, the prefix is dropped: this is what issue #236 is about.
$1/plink2 $2 $3 --vcf tmp_ucsc.vcf --psam tmp_samples.psam --allow-extra-chr --export vcf --out plink2_default
test "$(vcf_chrs plink2_default.vcf | tr '\n' ' ')" = "1 2 22 X Y MT chrUn_gl000220 1 2 22 X Y MT chrUn_gl000220 "

# Direct VCF -> VCF.
for f in ucsc ensembl plainm numeric ucscmt autosomes; do
    $1/plink2 $2 $3 --vcf tmp_$f.vcf --psam tmp_samples.psam --allow-extra-chr --output-chr infer --export vcf --out plink2_$f
    check_vcf tmp_$f.vcf plink2_$f.vcf
done

# VCF -> .pgen -> VCF round trip, with infer at each step.
$1/plink2 $2 $3 --vcf tmp_ucsc.vcf --psam tmp_samples.psam --allow-extra-chr --output-chr infer --make-pgen --out plink2_rt
test "$(pvar_chrs plink2_rt.pvar | tr '\n' ' ')" = "chr1 chr2 chr22 chrX chrY chrM chrUn_gl000220 "
$1/plink2 $2 $3 --pfile plink2_rt --allow-extra-chr --output-chr infer --export vcf --out plink2_rt
check_vcf tmp_ucsc.vcf plink2_rt.vcf

# .bim input and output.
$1/plink2 $2 $3 --vcf tmp_ucsc.vcf --psam tmp_samples.psam --allow-extra-chr --output-chr infer --make-bed --out plink2_bim
test "$(cut -f 1 plink2_bim.bim | tr '\n' ' ')" = "chr1 chr2 chr22 chrX chrY chrM chrUn_gl000220 "
$1/plink2 $2 $3 --bfile plink2_bim --allow-extra-chr --output-chr infer --make-just-pvar --out plink2_bim
test "$(pvar_chrs plink2_bim.pvar | tr '\n' ' ')" = "chr1 chr2 chr22 chrX chrY chrM chrUn_gl000220 "
# A plain .bim stays plain.
$1/plink2 $2 $3 --vcf tmp_ensembl.vcf --psam tmp_samples.psam --make-bed --out plink2_plainbim
$1/plink2 $2 $3 --bfile plink2_plainbim --output-chr infer --make-just-pvar --out plink2_plainbim
test "$(pvar_chrs plink2_plainbim.pvar | tr '\n' ' ')" = "1 2 22 X Y MT "

# BCF input.  (The BCF's own contig names follow --output-chr as well.)
$1/plink2 $2 $3 --vcf tmp_ucsc.vcf --psam tmp_samples.psam --allow-extra-chr --output-chr infer --export bcf --out plink2_bcf
$1/plink2 $2 $3 --bcf plink2_bcf.bcf --psam tmp_samples.psam --allow-extra-chr --output-chr infer --export vcf --out plink2_bcf
check_vcf tmp_ucsc.vcf plink2_bcf.vcf
$1/plink2 $2 $3 --bcf plink2_bcf.bcf --psam tmp_samples.psam --allow-extra-chr --export vcf --out plink2_bcfdefault
test "$(grep -v '^#' plink2_bcfdefault.vcf | cut -f 1 | tr '\n' ' ')" = "1 2 22 X Y MT chrUn_gl000220 "

# Inconsistent 'chr' prefixes: the default (no prefix) is used, while the
# consistent 'M' is still honored.
make_vcf tmp_mixed.vcf chr1 2 chrX chrM
$1/plink2 $2 $3 --vcf tmp_mixed.vcf --psam tmp_samples.psam --output-chr infer --export vcf --out plink2_mixed
test "$(grep -v '^#' plink2_mixed.vcf | cut -f 1 | tr '\n' ' ')" = "1 2 X M "

# A numeric X next to a lettered M: the default (letters) is used for X, while
# the M spelling is still honored.
make_vcf tmp_mixed2.vcf chr1 chr23 chrM
$1/plink2 $2 $3 --vcf tmp_mixed2.vcf --psam tmp_samples.psam --output-chr infer --export vcf --out plink2_mixed2
test "$(grep -v '^#' plink2_mixed2.vcf | cut -f 1 | tr '\n' ' ')" = "chr1 chrX chrM "

# No standard chromosome codes at all: default.
make_vcf tmp_extra.vcf chrUn_gl000220 contigA
$1/plink2 $2 $3 --vcf tmp_extra.vcf --psam tmp_samples.psam --allow-extra-chr --output-chr infer --export vcf --out plink2_extra
check_vcf tmp_extra.vcf plink2_extra.vcf

# Other output files follow the inferred codes too.
$1/plink2 $2 $3 --vcf tmp_ucsc.vcf --psam tmp_samples.psam --allow-extra-chr --output-chr infer --freq --out plink2_freq
test "$(grep -v '^#' plink2_freq.afreq | cut -f 1 | tr '\n' ' ')" = "chr1 chr2 chr22 chrX chrY chrM chrUn_gl000220 "

fails() {
    "$@" && false || true
}

fails $1/plink2 $2 $3 --vcf tmp_ensembl.vcf --psam tmp_samples.psam --output-chr inferred --export vcf --out plink2_bad
