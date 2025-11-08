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



#---------- Run PLINK2 GLM tests ----------

./CI-SCRIPTS/LOCAL_PLINK2_GLM_TESTING.sh


## If run locally, cleanup
rm -rf venv results