#!/bin/bash
# ==========================================================
# CONFIGURATION FILE for PLINK2 vs R GWAS GLM Testing
# ==========================================================
# This config file loads batch parameters from CI_reference_files
# and sets local paths for your testing environment.
# ==========================================================

# -------------------
#  PATH CONFIGURATION
# -------------------
# Base folder containing PLINK2 .pfile and phenotype files
datapath="test_data/"

# Optional subpath (if files live in subdirectories)
sbpath=""

# Output folder for PLINK2 results
outdir="./results"

# Path to Python comparison script
compare_script="./2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py"

# Correlation threshold for pass/fail (must be between 0 and 1)
# Tests will FAIL if any correlation is below this threshold
correlation_threshold=0.8

# -------------------
#  BATCH SELECTION
# -------------------
# Select which batch to run (1-12)
# Can be overridden by GLM_BATCH_NUM environment variable (used in CI)
BATCH_NUM=${GLM_BATCH_NUM:-1}

echo "Using Batch Number: $BATCH_NUM"

# -------------------
#  LOAD BATCH PARAMETERS
# -------------------
# Load test configurations from the selected batch file
BATCH_FILE="./2.0/build_dynamic/CI-SCRIPTS/CI_reference_files/BATCH_$(printf %02d $BATCH_NUM)_PARAMS.txt"

if [[ ! -f "$BATCH_FILE" ]]; then
    echo "❌ Batch file not found: $BATCH_FILE"
    echo "   Run 'python3 generate_batch_configs.py' first to generate batch files."
    exit 1
fi

# Read parameters from batch file
params=()
while IFS= read -r line; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    params+=("$line")
done < "$BATCH_FILE"

echo "✓ Loaded ${#params[@]} tests from Batch $BATCH_NUM"

# -------------------
#  PARAMETER FORMAT
# -------------------
# Each parameter line contains:
#   dataset_name,model,n,covariate,threads
#
# Parameters:
#   dataset_name : root filename prefix (no extension)
#                  e.g. 1kgp3_50k_nomiss_Av_nonintdose
#                       1kgp3_50k_yesmiss_Av_nonintdose
#
#   model        : one of
#                    linear   → phenotype = y
#                    logistic → phenotype = ybool
#                    firth    → phenotype = ybool, hybrid logistic
#
#   n            : subset size (matches .keep file) (1000, 32000, 32017)
#
#   covariate    : covariate column in phenotype file (e.g. "COV_1")
#                  leave blank for no covariates
#
#   threads      : number of threads for PLINK2 (e.g. 1, 4, 8)
# -------------------