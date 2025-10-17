## Data Conditions - General for All Tests
# Test Options (Aside from Core Test Data):
# General Conditions 1 - Sample Size:
# - Large data with non-64 multiple and odd (size of the original test data): --thin-indiv-count $((64*500+17))
# - Large data with multiples of 64: --thin-indiv-count $((64*500))
# - Small data: --thin-indiv-count $((1000))
# General Conditions 2 - Threads:
# - Single thread: --threads 1
# - Multiple threads: --threads 4

# Define bfiles
bfiles_=(
  "1kgp3_50k_nomiss_Av_nonintdose" 
  "1kgp3_50k_yesmiss_Av_nonintdose"
)

# Define general conditions
general_conditions_1=(
  "--thin-indiv-count $((1000))" 
  "--thin-indiv-count $((64*500+17))" 
  "--thin-indiv-count $((64*500))"
)
general_conditions_2=(
  "--threads 1" 
  "--threads 4"
)

# GWAS-Specific Conditions
pheno_=(
  "y" 
  "ybool"
)
covar_=(
  0 
  1
)

# Nested Loops for Combinations of Conditions
for bfile in "${bfiles_[@]}"; do
    for d1 in "${general_conditions_1[@]}"; do
        for d2 in "${general_conditions_2[@]}"; do
            for pheno in "${pheno_[@]}"; do
                for covar in "${covar_[@]}"; do
                    echo "Running analysis with:"
                    echo "  BFILE: $bfile"
                    echo "  Condition 1: $d1"
                    echo "  Condition 2: $d2"
                    echo "  Phenotype: $pheno"
                    echo "  Covariate: $covar"
                    # Add your GWAS command here, e.g.:
                    # plink2 --bfile $bfile --pheno $pheno $d1 $d2 --covar-number $covar ...
                done
            done
        done
    done
done
