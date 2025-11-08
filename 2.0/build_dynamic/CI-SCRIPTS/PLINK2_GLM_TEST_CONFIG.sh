#!/bin/bash
# ==========================================================
# CONFIGURATION FILE for PLINK2 vs R GWAS GLM Testing Pipeline
# ==========================================================
# Defines all datasets, models, covariates, subset sizes, and threads.
# Each line in `params` defines one full test run.
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
compare_script="2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py"

# -------------------
#  PARAMETER SCHEMA
# -------------------
# Each line in params[] follows:
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
#   n            : subset size (matches .keep file)
#
#   covariate    : covariate column in phenotype file (e.g. "COV_1")
#                  leave blank for no covariates
#
#   threads      : number of threads for PLINK2 (e.g. 1, 4, 8)
#
# Lines starting with "#" are ignored.
# -------------------

params=(
  # --- No-missing dataset tests ---
  "1kgp3_50k_nomiss_Av_nonintdose,linear,1000,,1"
  "1kgp3_50k_nomiss_Av_nonintdose,logistic,1000,,4"
  "1kgp3_50k_nomiss_Av_nonintdose,firth,32017,COV_1,4"

  # --- Yes-missing dataset tests ---
  "1kgp3_50k_yesmiss_Av_nonintdose,linear,1000,,1"
  "1kgp3_50k_yesmiss_Av_nonintdose,logistic,32000,,4"
  "1kgp3_50k_yesmiss_Av_nonintdose,firth,32000,COV_1,4"

  # Example (disabled):
  # "1kgp3_50k_nomiss_Av_nonintdose,firth,1000,,4"
)
