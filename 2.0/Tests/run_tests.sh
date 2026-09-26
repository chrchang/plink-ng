#!/bin/bash

# Usage: ./run_tests.sh {plink2 + pgen_compress build dir}
#   {up to 2 args, e.g. --randmem, "--threads 1"}
# Requires plink to be in the system PATH.

set -exo pipefail

if [[ $# -eq 0 ]]; then
    d=../../build_dynamic
else
    # doesn't always work, but should be good enough
    DIR=$1
    if [ "${DIR:0:1}" = "/" ]; then
        d=$DIR
    else
        d=../$DIR
    fi
fi

cd TEST_EXTRACT_CHR
./run_tests.sh $d $2 $3 > TEST_EXTRACT_CHR.log
cd ..
echo "TEST_EXTRACT_CHR passed."

cd TEST_MAF_FILTER
./run_tests.sh $d $2 $3 > TEST_MAF_FILTER.log
cd ..
echo "TEST_MAF_FILTER passed."

cd TEST_PGEN_FREQ
./run_tests.sh $d $2 $3 > TEST_PGEN_FREQ.log
cd ..
echo "TEST_PGEN_FREQ passed."

cd TEST_PHASED_VCF
./run_tests.sh $d $2 $3 > TEST_PHASED_VCF.log
cd ..
echo "TEST_PHASED_VCF passed."

cd TEST_SAMPLE_SUBSET
./run_tests.sh $d $2 $3 > TEST_SAMPLE_SUBSET.log
cd ..
echo "TEST_SAMPLE_SUBSET passed."

cd TEST_DOSAGE_ROUND_TRIP
./run_tests.sh $d $2 $3 > TEST_DOSAGE_ROUND_TRIP.log
cd ..
echo "TEST_DOSAGE_ROUND_TRIP passed."

cd TEST_EXPORT_OPTIONS
./run_tests.sh $d $2 $3 > TEST_EXPORT_OPTIONS.log
cd ..
echo "TEST_EXPORT_OPTIONS passed."

cd TEST_ONE_WAY_EXPORT
./run_tests.sh $d $2 $3 > TEST_ONE_WAY_EXPORT.log
cd ..
echo "TEST_ONE_WAY_EXPORT passed."

cd TEST_PMERGE
./run_tests.sh $d $2 $3 > TEST_PMERGE.log
cd ..
echo "TEST_PMERGE passed."

cd TEST_EXPORT_LEGACY
./run_tests.sh $d $2 $3 > TEST_EXPORT_LEGACY.log
cd ..
echo "TEST_EXPORT_LEGACY passed."

cd TEST_HOMOZYG
./run_tests.sh $d $2 $3 > TEST_HOMOZYG.log
cd ..
echo "TEST_HOMOZYG passed."

cd TEST_23FILE
./run_tests.sh $d $2 $3 > TEST_23FILE.log
cd ..
echo "TEST_23FILE passed."

cd TEST_HET_IBC
./run_tests.sh $d $2 $3 > TEST_HET_IBC.log
cd ..
echo "TEST_HET_IBC passed."

cd TEST_TWOLOCUS
./run_tests.sh $d $2 $3 > TEST_TWOLOCUS.log
cd ..
echo "TEST_TWOLOCUS passed."

cd TEST_DISTANCE
./run_tests.sh $d $2 $3 > TEST_DISTANCE.log
cd ..
echo "TEST_DISTANCE passed."

cd TEST_TEST_MISSING
./run_tests.sh $d $2 $3 > TEST_TEST_MISSING.log
cd ..
echo "TEST_TEST_MISSING passed."

cd TEST_FLIP_SCAN
./run_tests.sh $d $2 $3 > TEST_FLIP_SCAN.log
cd ..
echo "TEST_FLIP_SCAN passed."

cd TEST_MGF
./run_tests.sh $d $2 $3 > TEST_MGF.log
cd ..
echo "TEST_MGF passed."

cd TEST_KING_RT
./run_tests.sh $d $2 $3 > TEST_KING_RT.log
cd ..
echo "TEST_KING_RT passed."

cd TEST_SHOW_TAGS
./run_tests.sh $d $2 $3 > TEST_SHOW_TAGS.log
cd ..
echo "TEST_SHOW_TAGS passed."

cd TEST_META_ANALYSIS
./run_tests.sh $d $2 $3 > TEST_META_ANALYSIS.log
cd ..
echo "TEST_META_ANALYSIS passed."

cd TEST_PGEN_MALFORMED
./run_tests.sh $d $2 $3 > TEST_PGEN_MALFORMED.log
cd ..
echo "TEST_PGEN_MALFORMED passed."

cd TEST_MAKE_PHENO
./run_tests.sh $d $2 $3 > TEST_MAKE_PHENO.log
cd ..
echo "TEST_MAKE_PHENO passed."

cd TEST_MAKE_PERM_PHENO
./run_tests.sh $d $2 $3 > TEST_MAKE_PERM_PHENO.log
cd ..
echo "TEST_MAKE_PERM_PHENO passed."

cd TEST_Q_SCORE_RANGE
./run_tests.sh $d $2 $3 > TEST_Q_SCORE_RANGE.log
cd ..
echo "TEST_Q_SCORE_RANGE passed."

cd TEST_BLOCKS
./run_tests.sh $d $2 $3 > TEST_BLOCKS.log
cd ..
echo "TEST_BLOCKS passed."

cd TEST_PHENO_UTILS
./run_tests.sh $d $2 $3 > TEST_PHENO_UTILS.log
cd ..
echo "TEST_PHENO_UTILS passed."

cd TEST_GRM_MAF
./run_tests.sh $d $2 $3 > TEST_GRM_MAF.log
cd ..
echo "TEST_GRM_MAF passed."

cd TEST_EPISTASIS
./run_tests.sh $d $2 $3 > TEST_EPISTASIS.log
cd ..
echo "TEST_EPISTASIS passed."

cd TEST_TUCC
./run_tests.sh $d $2 $3 > TEST_TUCC.log
cd ..
echo "TEST_TUCC passed."

cd TEST_WRITE_VAR_RANGES
./run_tests.sh $d $2 $3 > TEST_WRITE_VAR_RANGES.log
cd ..
echo "TEST_WRITE_VAR_RANGES passed."

cd TEST_LD_SCORE_ANNOT
./run_tests.sh $d $2 $3 > TEST_LD_SCORE_ANNOT.log
cd ..
echo "TEST_LD_SCORE_ANNOT passed."

cd TEST_MISSING_FUSION
./run_tests.sh $d $2 $3 > TEST_MISSING_FUSION.log
cd ..
echo "TEST_MISSING_FUSION passed."

cd TEST_ATTRIB
./run_tests.sh $d $2 $3 > TEST_ATTRIB.log
cd ..
echo "TEST_ATTRIB passed."

cd TEST_NEIGHBOUR
./run_tests.sh $d $2 $3 > TEST_NEIGHBOUR.log
cd ..
echo "TEST_NEIGHBOUR passed."

cd TEST_SETS
./run_tests.sh $d $2 $3 > TEST_SETS.log
cd ..
echo "TEST_SETS passed."

cd TEST_TEST_MISHAP
./run_tests.sh $d $2 $3 > TEST_TEST_MISHAP.log
cd ..
echo "TEST_TEST_MISHAP passed."

cd TEST_EPISTASIS_QT
./run_tests.sh $d $2 $3 > TEST_EPISTASIS_QT.log
cd ..
echo "TEST_EPISTASIS_QT passed."

cd TEST_OBLIG_MISSING
./run_tests.sh $d $2 $3 > TEST_OBLIG_MISSING.log
cd ..
echo "TEST_OBLIG_MISSING passed."

cd TEST_GLM_MISSING
./run_tests.sh $d $2 $3 > TEST_GLM_MISSING.log
cd ..
echo "TEST_GLM_MISSING passed."

cd TEST_PGEN_EMPTY_DOSAGE_LIST
./run_tests.sh $d $2 $3 > TEST_PGEN_EMPTY_DOSAGE_LIST.log
cd ..
echo "TEST_PGEN_EMPTY_DOSAGE_LIST passed."

cd TEST_PGEN_CORRUPT_BODY
./run_tests.sh $d $2 $3 > TEST_PGEN_CORRUPT_BODY.log
cd ..
echo "TEST_PGEN_CORRUPT_BODY passed."

cd TEST_OXFORD_SAMPLE
./run_tests.sh $d $2 $3 > TEST_OXFORD_SAMPLE.log
cd ..
echo "TEST_OXFORD_SAMPLE passed."

cd TEST_KING_TABLE_MALFORMED
./run_tests.sh $d $2 $3 > TEST_KING_TABLE_MALFORMED.log
cd ..
echo "TEST_KING_TABLE_MALFORMED passed."

cd TEST_VCF_MALFORMED
./run_tests.sh $d $2 $3 > TEST_VCF_MALFORMED.log
cd ..
echo "TEST_VCF_MALFORMED passed."

cd TEST_BGEN_MALFORMED
./run_tests.sh $d $2 $3 > TEST_BGEN_MALFORMED.log
cd ..
echo "TEST_BGEN_MALFORMED passed."

cd TEST_PVAR_ALLELE_CT
./run_tests.sh $d $2 $3 > TEST_PVAR_ALLELE_CT.log
cd ..
echo "TEST_PVAR_ALLELE_CT passed."

cd TEST_ADJUST_FILE_PROVREF
./run_tests.sh $d $2 $3 > TEST_ADJUST_FILE_PROVREF.log
cd ..
echo "TEST_ADJUST_FILE_PROVREF passed."

cd TEST_SCORE_VS_FIXED_ALLELE
./run_tests.sh $d $2 $3 > TEST_SCORE_VS_FIXED_ALLELE.log
cd ..
echo "TEST_SCORE_VS_FIXED_ALLELE passed."

cd TEST_PCA_ALLELE_WTS_PROJECT
./run_tests.sh $d $2 $3 > TEST_PCA_ALLELE_WTS_PROJECT.log
cd ..
echo "TEST_PCA_ALLELE_WTS_PROJECT passed."

cd TEST_KING_TABLE_IBS1
./run_tests.sh $d $2 $3 > TEST_KING_TABLE_IBS1.log
cd ..
echo "TEST_KING_TABLE_IBS1 passed."

cd TEST_CHRY_DOSAGE_FREQ
./run_tests.sh $d $2 $3 > TEST_CHRY_DOSAGE_FREQ.log
cd ..
echo "TEST_CHRY_DOSAGE_FREQ passed."

cd TEST_SCORE_CENTER_MISSING
./run_tests.sh $d $2 $3 > TEST_SCORE_CENTER_MISSING.log
cd ..
echo "TEST_SCORE_CENTER_MISSING passed."

cd TEST_GLM_LOGISTIC_MISSING
./run_tests.sh $d $2 $3 > TEST_GLM_LOGISTIC_MISSING.log
cd ..
echo "TEST_GLM_LOGISTIC_MISSING passed."

cd TEST_BCF_MALFORMED
./run_tests.sh $d $2 $3 > TEST_BCF_MALFORMED.log
cd ..
echo "TEST_BCF_MALFORMED passed."

cd TEST_VCF_MULTIALLELIC_DS_FIELDS
./run_tests.sh $d $2 $3 > TEST_VCF_MULTIALLELIC_DS_FIELDS.log
cd ..
echo "TEST_VCF_MULTIALLELIC_DS_FIELDS passed."

cd TEST_VCF_MULTIALLELIC_DS_X
./run_tests.sh $d $2 $3 > TEST_VCF_MULTIALLELIC_DS_X.log
cd ..
echo "TEST_VCF_MULTIALLELIC_DS_X passed."

cd TEST_BCF_HAPLOID_HDS
./run_tests.sh $d $2 $3 > TEST_BCF_HAPLOID_HDS.log
cd ..
echo "TEST_BCF_HAPLOID_HDS passed."

cd TEST_BCF_UNPHASED_DOSAGE_HDS
./run_tests.sh $d $2 $3 > TEST_BCF_UNPHASED_DOSAGE_HDS.log
cd ..
echo "TEST_BCF_UNPHASED_DOSAGE_HDS passed."

cd TEST_BCF_MULTIALLELIC_HDS
./run_tests.sh $d $2 $3 > TEST_BCF_MULTIALLELIC_HDS.log
cd ..
echo "TEST_BCF_MULTIALLELIC_HDS passed."

cd TEST_BGEN_REF_FIRST_PHASE
./run_tests.sh $d $2 $3 > TEST_BGEN_REF_FIRST_PHASE.log
cd ..
echo "TEST_BGEN_REF_FIRST_PHASE passed."

cd TEST_SPLIT_PHASED
./run_tests.sh $d $2 $3 > TEST_SPLIT_PHASED.log
cd ..
echo "TEST_SPLIT_PHASED passed."

cd TEST_PGEN_DIFFLIST_COMMON
./run_tests.sh $d $2 $3 > TEST_PGEN_DIFFLIST_COMMON.log
cd ..
echo "TEST_PGEN_DIFFLIST_COMMON passed."

echo "All tests passed."
