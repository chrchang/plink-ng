#!/bin/bash

set -exo pipefail

$1/plink2 $2 $3 --dummy 33 65537 0.1 acgt dosage-freq=0.1 --out tmp_data

$1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-dosage=DS --out tmp_data2
$1/plink2 $2 $3 --vcf tmp_data2.vcf dosage=DS --out tmp_data2
diff -q tmp_data.pgen tmp_data2.pgen

$1/plink2 $2 $3 --pfile tmp_data --export vcf vcf-dosage=GP --out tmp_data3
$1/plink2 $2 $3 --vcf tmp_data3.vcf dosage=GP --out tmp_data3
diff -q tmp_data.pgen tmp_data3.pgen

$1/plink2 $2 $3 --pfile tmp_data --export bgen-1.1 --out tmp_data4
cat tmp_data4.sample | cut -d ' ' -f 1-3 > tmp_data4.min.sample
$1/plink2 $2 $3 --bgen tmp_data4.bgen ref-last --sample tmp_data4.min.sample --out tmp_data4
diff -q tmp_data.pgen tmp_data4.pgen

$1/plink2 $2 $3 --pfile tmp_data --export bgen-1.2 --out tmp_data5
$1/plink2 $2 $3 --bgen tmp_data5.bgen ref-last --out tmp_data5
diff -q tmp_data.pgen tmp_data5.pgen

$1/plink2 $2 $3 --pfile tmp_data --export oxford --out tmp_data6
$1/plink2 $2 $3 --data tmp_data6 ref-last --out tmp_data6
diff -q tmp_data.pgen tmp_data6.pgen

$1/plink2 $2 $3 --pfile tmp_data --export A-transpose --out tmp_data7
$1/plink2 $2 $3 --import-dosage tmp_data7.traw id-delim=_ skip0=1 skip1=2 chr-col-num=1 pos-col-num=4 ref-first --psam tmp_data.psam --out tmp_data7
diff -q tmp_data.pgen tmp_data7.pgen

# $1/plink2 $2 $3 --seed 1 --dummy 33 65537 0.1 acgt phase-freq=0.1 dosage-freq=0.1 --out tmp_data8
$1/plink2 $2 $3 --dummy 33 65537 0.1 acgt phase-freq=0.1 dosage-freq=0.1 --out tmp_data8

$1/plink2 $2 $3 --pfile tmp_data8 --export vcf vcf-dosage=HDS --out tmp_data8
$1/plink2 $2 $3 --vcf tmp_data8.vcf dosage=HDS --make-pgen --export vcf vcf-dosage=HDS --out tmp_data9
# Can't compare .pgens for equality for now, because --dummy (and
# MakePgenThread) don't yet try their best to omit inferrable HDS values.
# Probably want to clean this up before beta testing begins.
diff -q tmp_data8.vcf tmp_data9.vcf

$1/plink2 $2 $3 --pfile tmp_data8 --export bcf vcf-dosage=HDS-force --out tmp_data8
$1/plink2 $2 $3 --bcf tmp_data8.bcf dosage=HDS --out tmp_data10
# VCF-import and BCF-import are both smart about inferrable HDS.
diff -q tmp_data9.pgen tmp_data10.pgen
