#!/bin/bash
set -e
# # ---------- Build PLINK 2.0 ----------
# echo "Building PLINK 2.0 dynamically..."
# cd 2.0/build_dynamic
# make clean
# make
# cp plink2 "../../bin/plink2"

# # ---------- Add built/downloaded binaries to PATH ----------
# export PATH="../../bin/:$PATH"
# echo "Added $BIN_DIR to PATH"


## ---------- DL Data----------
## Remove existing venv if present
# rm -rf venv

# python3 -m venv venv
# source venv/bin/activate
# python3 -m pip install --upgrade pip
# pip install gdown

# echo "Downloading test data..."
# GDRIVE_FILE_ID="17x0g1SSzmkjEhuKapV192Ym3ejhB_FkG"
# gdown "https://drive.google.com/uc?id=$GDRIVE_FILE_ID" -O test_data.zip
# unzip -q test_data.zip -d test_data

# # ---------- Generate R results ----------


# fname="test_data/1kgp3_50k_nomiss_Av_nonintdose"

# plink2 \
#     --pfile "$fname" \
#     --make-pgen \
#     --set-all-var-ids @:#\$r,\$a \
#     --out "${fname}_recode_varIDs"

# plink2 \
#   --pfile "${fname}_recode_varIDs" \
#   --export A \
#   --out "${fname}_recode_varIDs_A"

# Rscript CI-SCRIPTS/R_gwas_tests.R

# # ---------- Run PLINK2 GLM tests ----------
echo "Running PLINK2 GLM tests..."

## GLM Tests

## Logistic - no firth, all samples, and one covariate

# datapath="test_data/" 
# pfile="1kgp3_50k_nomiss_Av_nonintdose_recode_varIDs"
# phenotype="ybool"
# nofirth="no-firth"
# phenofile="1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv"
# cov="COV_1"
# d1="" #--thin-indiv-count $((1000))"
# d2="--threads 4"

# plink2 --glm $nofirth hide-covar \
#   --pfile "$datapath$pfile" \
#   --allow-extra-chr \
#   --pheno "$datapath$phenofile" \
#   --pheno-name "$phenotype" \
#   --covar "$datapath$phenofile" \
#   --covar-name $cov\
#   $d1 $d2 \
#   --out "./results/${phenofile}_${phenotype}_${nofirth}_glm"

## Logistic - Firth only, all samples, and one covariate

# datapath="test_data/" 
# pfile="1kgp3_50k_nomiss_Av_nonintdose_recode_varIDs"
# phenotype="ybool"
# nofirth="firth"
# phenofile="1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv"
# cov="COV_1"
# d1="" #--thin-indiv-count $((1000))"
# d2="--threads 4"

# plink2 --glm $nofirth hide-covar \
#   --pfile "$datapath$pfile" \
#   --allow-extra-chr \
#   --pheno "$datapath$phenofile" \
#   --pheno-name "$phenotype" \
#   --covar "$datapath$phenofile" \
#   --covar-name $cov\
#   $d1 $d2 \
#   --out "./results/${phenofile}_${phenotype}_${nofirth}_glm"

## Linear, all samples, and one covariate

# datapath="test_data/" 
# pfile="1kgp3_50k_nomiss_Av_nonintdose_recode_varIDs"
# phenotype="y"
# nofirth=""
# phenofile="1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv"
# cov="COV_1"
# d1="" #--thin-indiv-count $((1000))"
# d2="--threads 4"

# plink2 --glm $nofirth hide-covar \
#   --pfile "$datapath$pfile" \
#   --allow-extra-chr \
#   --pheno "$datapath$phenofile" \
#   --pheno-name "$phenotype" \
#   --covar "$datapath$phenofile" \
#   --covar-name $cov\
#   $d1 $d2 \
#   --out "./results/${phenofile}_${phenotype}_${nofirth}_glm"



# # ---------- Compare PLINK2 with R  ----------    

## Load both R and PLINK2 results in python and make comparisons

# plink_results="results/1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv_ybool_no-firth_glm.ybool.glm.logistic"
# R_results="results/1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv_ybool_COV_1_glm_logistic_logistic.csv"
# plink_results="results/1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv_y__glm.y.glm.linear"
# R_results="results/1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv_y_COV_1_glm_linear_linear.csv"


plink_results="results/1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv_ybool_firth_glm.ybool.glm.firth"
R_results="results/1kgp3_50k_nomiss_Av_nonintdose_ybool_COV_1_glm_firth_firth.csv"

python3 2.0/build_dynamic/CI-SCRIPTS/compare_plink2_R_glm.py $plink_results $R_results
