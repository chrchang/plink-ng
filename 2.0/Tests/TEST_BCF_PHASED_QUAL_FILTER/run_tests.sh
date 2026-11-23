#!/bin/bash

# --bcf with --vcf-min-gq/--vcf-min-dp must handle a phased call in the first
# sample.  BcfScanGt() used to look up that sample's first allele at
# [sample_idx * 2 - 1] relative to the second-allele pointer; with an
# unsigned sample_idx of 0 this indexed ~4 GiB past the GT vector instead of
# one element before it, and crashed on 64-bit builds.  Both the int8 and the
# int16 (>63 ALT alleles) GT paths are checked, with a sample count below and
# above the SIMD pre-scan width.  --bcf results must match --vcf on the same
# text.

set -exo pipefail

ALTS70=$(awk 'BEGIN {b = "ACGT"; for (i = 1; i < 70; ++i) {a = "C"; for (j = i; j; j = int(j / 4)) a = a substr(b, j % 4 + 1, 1); printf "%s%s", (i > 1)? "," : "", a}}')

# Columns after FORMAT for the first two samples; make_vcf pads the rest.
cat > body.txt << EOF2
1	100	v1	A	C	.	.	.	GT:GQ:DP	0|1:50:30	0/0:50:30
1	200	v2	A	C	.	.	.	GT:GQ:DP	0|1:5:2	0|1:50:30
1	300	v3	A	C	.	.	.	GT:GQ:DP	0|0:50:30	1|0:50:30
1	400	v4	A	C	.	.	.	GT:GQ:DP	0|1:5:2	0/1:50:30
1	500	v5	A	$ALTS70	.	.	.	GT:GQ:DP	0|69:50:30	0/0:50:30
1	600	v6	A	$ALTS70	.	.	.	GT:GQ:DP	1|69:5:2	2|0:50:30
EOF2

make_vcf() {
    local n=$1
    {
        printf '##fileformat=VCFv4.3\n##contig=<ID=1,length=1000000>\n'
        printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
        printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
        for i in $(seq 1 $n); do printf '\ts%d' $i; done
        printf '\n'
        awk -v n=$n 'BEGIN {OFS = FS = "\t"} {s = $0; for (i = 3; i <= n; ++i) s = s "\t0/0:50:30"; print s}' body.txt
    } > tmp_s$n.vcf
}

for n in 2 40; do
    make_vcf $n
    python3 ../utils/vcf_to_bcf.py tmp_s$n.vcf tmp_s$n.bcf
    for filt in "" "--vcf-min-gq 10" "--vcf-min-dp 10" "--vcf-min-gq 10 --vcf-min-dp 10"; do
        $1/plink2 $2 $3 --vcf tmp_s$n.vcf $filt --export vcf --out tmp_want
        $1/plink2 $2 $3 --bcf tmp_s$n.bcf $filt --export vcf --out tmp_got
        diff <(grep -v '^##' tmp_want.vcf) <(grep -v '^##' tmp_got.vcf)
    done
    # Spot-check the filtered calls: s1's failing phased hets become missing,
    # and s1's passing phased calls keep their phase.
    $1/plink2 $2 $3 --bcf tmp_s$n.bcf --vcf-min-gq 10 --export vcf --out tmp_got
    test "$(grep -v '^#' tmp_got.vcf | cut -f 10 | tr '\n' ' ')" = "0|1 .|. 0|0 .|. 0|69 .|. "
    test "$(grep -v '^#' tmp_got.vcf | cut -f 11 | tr '\n' ' ')" = "0|0 0|1 1|0 0/1 0/0 2|0 "
done
