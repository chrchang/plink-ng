#!/bin/bash

# --score variance-standardize on a variant whose scored allele is fixed.
#
# A variant with zero variance contributes nothing to a variance-standardized
# score.  When the scored allele is absent, that already worked; when it is
# fixed instead (e.g. the REF allele of a variant with no ALT call), --score
# used to stop with an error claiming the allele frequency was zero while
# some dosages were not.  Frequencies that genuinely contradict the genotypes
# must still be rejected, in either direction.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# v2 has no ALT call.  v4 is multiallelic with only ALT2 calls, so REF and
# ALT1 are absent and ALT2 is fixed.
cat > tmp_data.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4	s5	s6
1	100	v1	A	G	.	.	.	GT	0/1	0/0	1/1	0/1	0/0	0/1
1	200	v2	C	T	.	.	.	GT	0/0	0/0	0/0	./.	0/0	0/0
1	300	v3	G	A	.	.	.	GT	0/0	1/1	0/1	0/1	1/1	0/0
1	400	v4	A	C,G	.	.	.	GT	2/2	2/2	2/2	2/2	./.	2/2
EOF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

# $1: score file lines.  $2: output prefix.  Remaining args: extra flags.
run_score() {
    printf "$1" > $2.txt
    local out=$2
    shift 2
    $plink2 --pfile tmp_data --score $out.txt variance-standardize cols=+scoresums "$@" --out $out > /dev/null
}

# The zero-variance alleles, fixed or absent, must leave the score sums of v1
# and v3 alone.
run_score 'v1\tA\t0.5\nv3\tA\t-1.25\n' tmp_base --bad-freqs
for extra in 'v2\tC\t2\n' 'v2\tT\t3\n' 'v4\tG\t4\n' 'v4\tA\t5\n' 'v4\tC\t6\n'; do
    run_score "v1\tA\t0.5\nv3\tA\t-1.25\n$extra" tmp_extra --bad-freqs
    awk 'FNR == 1 { next }
         NR == FNR { s[$1] = $NF; next }
         { if (s[$1] != $NF) { print "sample " $1 ": " s[$1] " vs " $NF; exit 1 }; ++n }
         END { if (n != 6) { print "expected 6 samples, got " n + 0; exit 1 } }' tmp_base.sscore tmp_extra.sscore
done

# Frequencies that contradict the genotypes are still an error: v2's ALT
# frequency is claimed to be 1, so neither of its alleles can be scored.
printf '#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\tv1\tA\tG\t0.4\t12\n1\tv2\tC\tT\t1\t10\n1\tv3\tG\tA\t0.5\t12\n' > tmp_bad.afreq
for allele in C T; do
    if run_score "v2\t$allele\t1\n" tmp_bad --read-freq tmp_bad.afreq; then
        echo "v2:$allele was scored despite contradictory frequencies"
        exit 1
    fi
    grep -q 'variance-standardize failure' tmp_bad.log
done
