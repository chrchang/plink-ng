#!/bin/bash
set -e

echo "Building Plink2..."

## goto build directory
cd 2.0/build_dynamic
make clean
make
sudo cp plink2 /usr/local/bin/plink2

## return to root directory

cd ../../



