#!/bin/bash

# Phenotype/covariate transform flags with no other test coverage:
# --variance-standardize, --quantile-normalize, --pheno-quantile-normalize,
# and --covar-quantile-normalize (alongside --covar-variance-standardize).
# Results are compared against an independent reimplementation.

set -exo pipefail

plink2="$1/plink2 $2 $3"

{
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\n'
    printf '1\t100\tv1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\t1/1\t0/1\t0/0\t0/1\n'
} > tmp_data.vcf
# Q1 has a tie (3) and an NA; B is binary; CAT is categorical.  C2 has a tie
# at its minimum.
printf '#IID\tSEX\tQ1\tQ2\tB\nS1\t1\t1\t10\t1\nS2\t2\t3\t20\t2\nS3\t1\t2\tNA\t1\nS4\t2\t10\t40\t2\nS5\t1\t3\t30\t1\nS6\t2\tNA\t50\t2\n' > tmp_data.psam
printf '#IID\tC1\tC2\tCAT\nS1\t5\t1\tx\nS2\t7\t1\ty\nS3\t9\t2\tx\nS4\t11\t3\ty\nS5\tNA\t5\tx\nS6\t15\t8\tz\n' > tmp_covar.txt
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --covar tmp_covar.txt --make-just-psam --write-covar --out tmp_orig > /dev/null

# $1: label; $2: phenotype specs; $3: covariate specs; remaining args:
# transform flags.
check() {
    local label=$1
    local pheno_specs=$2
    local covar_specs=$3
    shift 3
    $plink2 --pfile tmp_data --covar tmp_covar.txt "$@" --make-just-psam --write-covar --out tmp_$label > /dev/null
    python3 check_transform.py tmp_orig.psam tmp_$label.psam $pheno_specs
    python3 check_transform.py tmp_orig.cov tmp_$label.cov $covar_specs
}

# With no arguments, everything quantitative is transformed.
check vs_all "Q1=vs Q2=vs" "C1=vs C2=vs" --variance-standardize
check vs_some "Q1=vs" "C2=vs" --variance-standardize Q1 C2
check cvs "" "C1=vs C2=vs" --covar-variance-standardize
check qn_all "Q1=qn Q2=qn" "C1=qn C2=qn" --quantile-normalize
check qn_some "Q2=qn" "C1=qn" --quantile-normalize Q2 C1
check pqn "Q1=qn Q2=qn" "" --pheno-quantile-normalize
check pqn_some "Q1=qn" "" --pheno-quantile-normalize Q1
check cqn "" "C1=qn C2=qn" --covar-quantile-normalize
check cqn_some "" "C1=qn" --covar-quantile-normalize C1
# Mixed: quantile normalization happens before variance standardization.
check qn_vs "Q1=qn+vs Q2=qn+vs" "C1=qn+vs C2=qn+vs" --quantile-normalize --variance-standardize
check pqn_cvs "Q1=qn" "C1=vs C2=vs" --pheno-quantile-normalize Q1 --covar-variance-standardize
check pqn_cqn "Q1=qn Q2=qn" "C2=qn" --pheno-quantile-normalize --covar-quantile-normalize C2

# Spot-check one value by hand: Q1's tied 3s have mean rank 3.5 of 5, and
# Phi^-1(0.6) = 0.253347.
test "$(awk -F '\t' '$1 == "S2" {print $3}' tmp_pqn.psam)" = "0.253347"

# Conflicting flag combinations.
for flags in "--quantile-normalize --pheno-quantile-normalize" "--quantile-normalize --covar-quantile-normalize" "--variance-standardize --covar-variance-standardize"; do
    if $plink2 --pfile tmp_data --covar tmp_covar.txt $flags --make-just-psam --write-covar --out tmp_bad > /dev/null 2>&1; then
        exit 1
    fi
done
