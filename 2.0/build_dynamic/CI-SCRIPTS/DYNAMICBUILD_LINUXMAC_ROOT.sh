#!/bin/bash
set -e

echo "Setting up Python environment..."
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
pip install gdown


echo "Building Plink2..."

## goto build directory
cd 2.0/build_dynamic
make clean
make
sudo cp plink2 /usr/local/bin/plink2

## return to root directory

cd ../../



