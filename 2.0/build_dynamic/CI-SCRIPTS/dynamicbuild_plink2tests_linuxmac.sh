#!/bin/bash
set -e

echo "Setting up Python environment..."
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
pip install gdown


echo "Building Plink2..."
cd 2.0/build_dynamic
make clean
make
sudo cp plink2 /usr/local/bin/plink2

# echo "Downloading test data..."
# GDRIVE_FILE_ID="1duspKrxtYdf_dJe5jX4Y_GMkzmy6i8dj"
# gdown "https://drive.google.com/uc?id=$GDRIVE_FILE_ID" -O test_data.zip
# unzip test_data.zip -d test_data
# mv test_data/yesmiss_testdata/* test_data/
# rm -rf test_data/yesmiss_testdata

# echo "Running Plink2 tests..."
# mkdir -p ./derivatives/
# plink2 --pfile test_data/1kgp3_50k_yesmiss_Av_nonintdose --genotyping-rate dosage --out ./derivatives/tmp

# echo "Listing files..."
# ls ./test_data/*
# ls ./derivatives/*
