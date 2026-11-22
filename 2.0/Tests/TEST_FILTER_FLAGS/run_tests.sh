#!/bin/bash

# Sample/variant filter and allele-order flags with no other test coverage:
# --keep-nosex, --remove-females, --remove-males, --remove-nosex,
# --exclude-palindromic-snps, --min-alleles, --from/--to, --require-info,
# --require-no-info, --require-covar, --thin, --thin-count,
# --thin-indiv-count, --maj-ref, and --alt-allele.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# v2/v3/v4 are palindromic (v4 in lowercase); v5 is multiallelic with an A/T
# allele pair; v6 has a missing ALT; v7 is an indel.  DP=. counts as a
# missing INFO key, and DB is a flag key.
cat > tmp_data.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
1	100	v1	A	G	.	.	DP=5;DB	GT	1/1	1/1	0/1	1/1	0/0	1/1
1	200	v2	A	T	.	.	DP=.	GT	0/1	0/0	0/0	0/0	0/1	0/0
1	300	v3	C	G	.	.	DB	GT	0/0	0/1	0/0	0/0	0/0	0/0
1	400	v4	g	c	.	.	.	GT	0/0	0/0	0/1	0/0	0/0	0/0
1	500	v5	A	C,T	.	.	DP=7	GT	2/2	2/1	2/2	0/2	1/1	2/2
1	600	v6	T	.	.	.	DP=3	GT	0/0	0/0	0/0	0/0	0/0	0/0
1	700	v7	AT	A	.	.	DB	GT	0/1	0/0	0/0	0/0	0/0	0/0
1	800	v8	C	T	.	.	DP=9;DB	GT	0/0	0/1	0/0	0/0	0/0	0/0
VCF
printf '#IID\tSEX\tQ\nS1\t1\t1.5\nS2\t2\tNA\nS3\t0\t2.0\nS4\t1\t0.5\nS5\t2\t3.0\nS6\tNA\t1.0\n' > tmp_data.psam
printf '#IID\tC1\tC2\nS1\t1\t2\nS2\t3\tNA\nS3\tNA\t1\nS4\t4\t5\nS5\t2\t2\nS6\t0\t0\n' > tmp_covar.txt
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

# $1: expected space-separated IIDs; remaining args: plink2 flags.
check_samples() {
    local expected=$1
    shift
    $plink2 --pfile tmp_data "$@" --make-just-psam --out tmp_s > /dev/null
    test "$(tail -n +2 tmp_s.psam | cut -f 1 | tr '\n' ' ')" = "$expected "
}

# $1: expected space-separated variant IDs; remaining args: plink2 flags.
check_variants() {
    local expected=$1
    shift
    $plink2 --pfile tmp_data "$@" --write-snplist --out tmp_v > /dev/null
    test "$(tr '\n' ' ' < tmp_v.snplist)" = "$expected "
}

# $1: .gcount file, $2: variant ID, remaining args: column names.  Prints
# the requested fields, tab-separated.
gcount_fields() {
    local fname=$1
    local id=$2
    shift 2
    awk -F '\t' -v id="$id" -v cols="$*" 'BEGIN {OFS = "\t"; n = split(cols, want, " ")}
        NR == 1 {for (i = 1; i <= NF; i++) colidx[$i] = i; next}
        $2 == id {line = $colidx[want[1]]; for (j = 2; j <= n; j++) line = line OFS $colidx[want[j]]; print line}' "$fname"
}

# S3 has sex code 0 and S6 has NA.
check_samples "S3 S6" --keep-nosex
check_samples "S1 S3 S4 S6" --remove-females
check_samples "S2 S3 S5 S6" --remove-males
check_samples "S1 S2 S4 S5" --remove-nosex
check_samples "S1 S4" --remove-females --remove-nosex

# Only biallelic A/T and C/G SNPs are palindromic, regardless of case.
check_variants "v1 v5 v6 v7 v8" --exclude-palindromic-snps
# A lone missing ALT counts as one allele.
check_variants "v1 v2 v3 v4 v5 v7 v8" --min-alleles 2
check_variants "v5" --min-alleles 3
check_variants "v3 v4 v5 v6" --from v3 --to v6

# --require-info needs every key; --require-no-info excludes on any key.
check_variants "v1 v5 v6 v8" --require-info DP
check_variants "v1 v8" --require-info DP DB
check_variants "v2 v4 v5 v6" --require-no-info DB
check_variants "v2 v4" --require-no-info DP DB

# With no arguments, --require-covar requires every covariate.
check_samples "S1 S4 S5 S6" --covar tmp_covar.txt --require-covar
check_samples "S1 S3 S4 S5 S6" --covar tmp_covar.txt --require-covar C2
check_samples "S1 S4 S5 S6" --covar tmp_covar.txt --require-pheno Q --require-covar C1

# Random thinning: exact counts, and reproducible under a fixed seed.
for seed in 1 2; do
    $plink2 --pfile tmp_data --thin-count 3 --seed $seed --write-snplist --out tmp_tc$seed > /dev/null
    test "$(wc -l < tmp_tc$seed.snplist)" -eq 3
    $plink2 --pfile tmp_data --thin-count 3 --seed $seed --write-snplist --out tmp_tc${seed}b > /dev/null
    diff -q tmp_tc$seed.snplist tmp_tc${seed}b.snplist
    $plink2 --pfile tmp_data --thin-indiv-count 2 --seed $seed --make-just-psam --out tmp_tic$seed > /dev/null
    test "$(tail -n +2 tmp_tic$seed.psam | wc -l)" -eq 2
    $plink2 --pfile tmp_data --thin 0.5 --seed $seed --write-snplist --out tmp_t$seed > /dev/null
    $plink2 --pfile tmp_data --thin 0.5 --seed $seed --write-snplist --out tmp_t${seed}b > /dev/null
    diff -q tmp_t$seed.snplist tmp_t${seed}b.snplist
    # --thin keeps a subset, in the original order.
    grep -vxF -f tmp_t$seed.snplist <(cut -f 3 tmp_data.pvar | grep '^v') > tmp_t${seed}_removed.txt || true
    test "$(cat tmp_t$seed.snplist tmp_t${seed}_removed.txt | wc -l)" -eq 8
done
if $plink2 --pfile tmp_data --thin 1 --write-snplist --out tmp_t1x > /dev/null 2>&1; then
    exit 1
fi

# --maj-ref: VCF REF alleles are 'known', so without 'force' this is a
# no-op.  With 'force', v1 swaps, and v5's ALTs are sorted by decreasing
# frequency (T: 8, C: 3, A: 1).
$plink2 --pfile tmp_data --maj-ref --make-pgen --out tmp_mr0 > /dev/null
diff -q tmp_data.pvar tmp_mr0.pvar
$plink2 --pfile tmp_data --maj-ref force --make-pgen --out tmp_mrf > /dev/null
grep -v '^##' tmp_mrf.pvar | cut -f 3-5 | diff - <(printf 'ID\tREF\tALT\nv1\tG\tA\nv2\tA\tT\nv3\tC\tG\nv4\tg\tc\nv5\tT\tC,A\nv6\tT\t.\nv7\tAT\tA\nv8\tC\tT\n')
$plink2 --pfile tmp_mrf --geno-counts --out tmp_mrf > /dev/null
# Recoded under T,C,A: 0/0 0/1 0/0 0/2 1/1 0/0, and both swaps are now
# provisional.
test "$(gcount_fields tmp_mrf.gcount v5 PROVISIONAL_REF? HOM_REF_CT HET_REF_ALT_CTS TWO_ALT_GENO_CTS)" = "$(printf 'Y\t3\t1,1\t1,0,0')"
grep -c ';PR$' tmp_mrf.pvar | grep -qx 2

# Without 'force', only provisional-REF variants are touched: mark v1 (but
# not v5) provisional.
awk 'BEGIN {OFS="\t"} /^#CHROM/ {print "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">"} !/^#/ && $3 == "v1" {$6 = $6 ";PR"} {print}' tmp_data.pvar > tmp_pr.pvar
cp tmp_data.pgen tmp_pr.pgen
cp tmp_data.psam tmp_pr.psam
$plink2 --pfile tmp_pr --maj-ref --make-pgen --out tmp_mr > /dev/null
grep -v '^##' tmp_mr.pvar | cut -f 3-6 | head -n 2 | diff - <(printf 'ID\tREF\tALT\tINFO\nv1\tG\tA\tDP=5;DB;PR\n')
grep -q "$(printf '^1\t500\tv5\tA\tC,T\t')" tmp_mr.pvar

# --alt-allele: v1 and v7 already match; v5's ALTs are reordered (REF stays
# A, so no 'force' needed); v9 isn't in the dataset; '.' is skipped.
printf 'v1\tG\nv5\tT,C\nv7\tA\nv9\tC\nv2\t.\n' > tmp_alt1.txt
$plink2 --pfile tmp_data --alt-allele tmp_alt1.txt 2 1 --make-pgen --out tmp_aa > /dev/null
grep -v '^##' tmp_aa.pvar | cut -f 3-5 | diff - <(printf 'ID\tREF\tALT\nv1\tA\tG\nv2\tA\tT\nv3\tC\tG\nv4\tg\tc\nv5\tA\tT,C\nv6\tT\t.\nv7\tAT\tA\nv8\tC\tT\n')
# Genotypes follow: v5 was 2/2 2/1 2/2 0/2 1/1 2/2 under A,C,T.
$plink2 --pfile tmp_aa --geno-counts --out tmp_aa > /dev/null
test "$(gcount_fields tmp_aa.gcount v5 HOM_REF_CT HET_REF_ALT_CTS TWO_ALT_GENO_CTS)" = "$(printf '0\t1,0\t3,1,1')"

# Changing a known REF requires 'force', and marks the REF provisional.
printf 'v8\tC\n' > tmp_alt2.txt
if $plink2 --pfile tmp_data --alt-allele tmp_alt2.txt 2 1 --make-pgen --out tmp_ab > /dev/null 2>&1; then
    exit 1
fi
$plink2 --pfile tmp_data --alt-allele force tmp_alt2.txt 2 1 --make-pgen --out tmp_ab > /dev/null
grep -q "$(printf '^1\t800\tv8\tT\tC\tDP=9;DB;PR$')" tmp_ab.pvar
$plink2 --pfile tmp_ab --geno-counts --out tmp_ab > /dev/null
test "$(gcount_fields tmp_ab.gcount v8 REF ALT PROVISIONAL_REF? HOM_REF_CT HET_REF_ALT_CTS TWO_ALT_GENO_CTS)" = "$(printf 'T\tC\tY\t0\t1\t5')"
