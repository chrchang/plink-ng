set -e


echo "Setting up Python environment..."
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
pip install gdown

# # ---------- Download test data ----------
echo "Downloading test data..."
GDRIVE_FILE_ID="17x0g1SSzmkjEhuKapV192Ym3ejhB_FkG"
gdown "https://drive.google.com/uc?id=$GDRIVE_FILE_ID" -O test_data.zip
unzip -q test_data.zip -d test_data

# ---------- Run PLINK2 tests ----------
echo "Running PLINK2 GLM tests..."
mkdir -p ./derivatives/

# quick test to make sure plink2 is working
plink2 --pfile test_data/1kgp3_50k_yesmiss_Av_nonintdose \
       --genotyping-rate dosage \
       --out ./derivatives/tmp

## GLM Tests
# datapath="test_data/" 
# pfile="1kgp3_50k_nomiss_Av_nonintdose"
# phenotype="ybool"
# phenofile="1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv"
# d1="" #--thin-indiv-count $((1000))"
# d2="--threads 4"

# plink2 --glm \
#   --pfile "$datapath$pfile" \
#   --allow-extra-chr \
#   --pheno "$datapath$phenofile" \
#   --pheno-name "$phenotype" \
#   --covar "$datapath$phenofile" \
#   --covar-name "COV_1" "COV_2" "COV_3"\
#   $d1 $d2 \
#   --out "./derivatives/${phenofile}_${phenotype}_glm"


# # ---------- List files ----------
# echo "Test data files:"
# ls -l ./test_data
# echo "Derivatives files:"
# ls -l ./derivatives

# ---------- Cleanup ----------
rm -rf test_data derivatives
make -C 2.0/build_dynamic/ clean 

echo "✅ PLINK2 tests complete!"