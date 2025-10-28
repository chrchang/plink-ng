#!/bin/bash

set -e


echo "Setting up Python environment..."
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
pip install pandas numpy 


# ---------- Run PLINK2 tests ----------
echo "Running PLINK2 GLM tests..."
mkdir -p ./results/

# # quick test to make sure plink2 is working
# plink2 --pfile test_data/1kgp3_50k_yesmiss_Av_nonintdose \
#        --genotyping-rate dosage \
#        --out ./derivatives/tmp

# # ---------- Generate R results ----------
## Should be done already and part of test data. 
## Code kept for reference.


# Rscript CI-SCRIPTS/R_gwas_tests.R

#---------- Run PLINK2 GLM tests ----------

## Run GLM
## Logistic - no firth, all samples, and one covariate

gtype="yesmiss"

datapath="test_data/" 
pfile="1kgp3_50k_${gtype}_Av_nonintdose_recode_varIDs"
phenotype="y"
nofirth=""
phenofile="1kgp3_50k_${gtype}_Av_nonintdose_combined_phenocov.csv"
cov="COV_1"
d1="" #--thin-indiv-count $((1000))"
d2="--threads 4"

plink2 --glm $nofirth hide-covar \
  --pfile "$datapath$pfile" \
  --allow-extra-chr \
  --pheno "$datapath$phenofile" \
  --pheno-name "$phenotype" \
  --covar "$datapath$phenofile" \
  --covar-name $cov\
  $d1 $d2 \
  --out "./results/${phenofile}_${phenotype}_${nofirth}_glm"

## Compare to R results

# python 2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py \
#   "results/${phenofile}_${phenotype}_${nofirth}_glm.ybool.glm.logistic" \
#   "test_data/1kgp3_50k_${gtype}_Av_nonintdose_ybool_COV_1_glm_logistic_logistic.csv"

python 2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py \
  "results/${phenofile}_${phenotype}_${nofirth}_glm.y.glm.linear" \
  "test_data/1kgp3_50k_${gtype}_Av_nonintdose_y_COV_1_glm_linear_linear.csv"



# # ---------- List files ----------
# echo "Test data files:"
# ls -l ./test_data
# echo "Derivatives files:"
# ls -l ./derivatives


echo "✅ PLINK2 tests complete!"


## If run locally, cleanup
rm -rf venv results