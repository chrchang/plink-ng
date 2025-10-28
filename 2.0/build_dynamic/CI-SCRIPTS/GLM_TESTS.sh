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
datapath="./test_data/"


## Config
cov="COV_1"
d1="" # --thin-indiv-count $((1000))
d2="--threads 4"

gtype_=("nomiss" "yesmiss")
phenotypes=("y" "ybool")

# Temporary buffer for all Python comparison outputs
summary_output=$(mktemp)

for gtype in "${gtype_[@]}"; do
  pfile="1kgp3_50k_${gtype}_Av_nonintdose_recode_varIDs"
  phenofile="1kgp3_50k_${gtype}_Av_nonintdose_combined_phenocov.csv"

  for phenotype in "${phenotypes[@]}"; do

    if [ "$phenotype" == "ybool" ]; then
      for nofirth in "" "no-firth"; do
        if [ "$nofirth" == "no-firth" ]; then
          model0="logistic"
          model="logistic"
        else
          model0="logistic.hybrid"
          model="firth"
        fi

        outprefix="./results/${phenofile}_${phenotype}_${nofirth}_glm"
        glm_file="${outprefix}.${phenotype}.glm.${model0}"
        ref_file="test_data/1kgp3_50k_${gtype}_Av_nonintdose_${phenotype}_${cov}_glm_${model}_${model}.csv"

        echo "Running: $phenofile / $phenotype / $model0"
        plink2 --glm $nofirth hide-covar \
          --pfile "${datapath}${pfile}" \
          --allow-extra-chr \
          --pheno "${datapath}${phenofile}" \
          --pheno-name "$phenotype" \
          --covar "${datapath}${phenofile}" \
          --covar-name "$cov" \
          $d1 $d2 \
          --out "$outprefix" >/dev/null 2>&1

        # Append comparison output to buffer
        echo -e "\n### $phenofile / $phenotype / $model0 ###" >>"$summary_output"
        python 2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py \
          "$glm_file" "$ref_file" >>"$summary_output" 2>&1
      done
    else
      nofirth=""
      model="linear"

      outprefix="./results/${phenofile}_${phenotype}_${nofirth}_glm"
      glm_file="${outprefix}.${phenotype}.glm.${model}"
      ref_file="test_data/1kgp3_50k_${gtype}_Av_nonintdose_${phenotype}_${cov}_glm_${model}_${model}.csv"

      echo "Running: $phenofile / $phenotype / $model"
      plink2 --glm hide-covar \
        --pfile "${datapath}${pfile}" \
        --allow-extra-chr \
        --pheno "${datapath}${phenofile}" \
        --pheno-name "$phenotype" \
        --covar "${datapath}${phenofile}" \
        --covar-name "$cov" \
        $d1 $d2 \
        --out "$outprefix" >/dev/null 2>&1

      echo -e "\n### $phenofile / $phenotype / $model ###" >>"$summary_output"
      python 2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py \
        "$glm_file" "$ref_file" >>"$summary_output" 2>&1
    fi

  done
done

# Print the compiled summary
echo
echo "============================================"
echo "🔍 GLM Comparison Summary"
echo "============================================"
cat "$summary_output"
echo "============================================"
rm -f "$summary_output"


# # ---------- List files ----------
# echo "Test data files:"
# ls -l ./test_data
# echo "Derivatives files:"
# ls -l ./derivatives


echo "✅ PLINK2 tests complete!"


## If run locally, cleanup
rm -rf venv results