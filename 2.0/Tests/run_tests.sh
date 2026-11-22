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

cd TEST_PMERGE_MULTIPASS
./run_tests.sh $d $2 $3 > TEST_PMERGE_MULTIPASS.log
cd ..
echo "TEST_PMERGE_MULTIPASS passed."

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

cd TEST_COMMON_FLAGS
./run_tests.sh $d $2 $3 > TEST_COMMON_FLAGS.log
cd ..
echo "TEST_COMMON_FLAGS passed."

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

cd TEST_GENE_REPORT
./run_tests.sh $d $2 $3 > TEST_GENE_REPORT.log
cd ..
echo "TEST_GENE_REPORT passed."

cd TEST_INFO_TO_COLS
./run_tests.sh $d $2 $3 > TEST_INFO_TO_COLS.log
cd ..
echo "TEST_INFO_TO_COLS passed."

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

cd TEST_SCORE_DOMREC
./run_tests.sh $d $2 $3 > TEST_SCORE_DOMREC.log
cd ..
echo "TEST_SCORE_DOMREC passed."

cd TEST_GLM_MULTINOMIAL
./run_tests.sh $d $2 $3 > TEST_GLM_MULTINOMIAL.log
cd ..
echo "TEST_GLM_MULTINOMIAL passed."

cd TEST_GLM_LINEAR_BATCH
./run_tests.sh $d $2 $3 > TEST_GLM_LINEAR_BATCH.log
cd ..
echo "TEST_GLM_LINEAR_BATCH passed."

cd TEST_OUTPUT_CHR_INFER
./run_tests.sh $d $2 $3 > TEST_OUTPUT_CHR_INFER.log
cd ..
echo "TEST_OUTPUT_CHR_INFER passed."

cd TEST_HDS_FORCE_INTEGRAL
./run_tests.sh $d $2 $3 > TEST_HDS_FORCE_INTEGRAL.log
cd ..
echo "TEST_HDS_FORCE_INTEGRAL passed."

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

cd TEST_GEN_IMPORT_MT
./run_tests.sh $d $2 $3 > TEST_GEN_IMPORT_MT.log
cd ..
echo "TEST_GEN_IMPORT_MT passed."

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

cd TEST_ADJUST_EXTREME_P
./run_tests.sh $d $2 $3 > TEST_ADJUST_EXTREME_P.log
cd ..
echo "TEST_ADJUST_EXTREME_P passed."

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

cd TEST_MULTIALLELIC_PHASED_R
./run_tests.sh $d $2 $3 > TEST_MULTIALLELIC_PHASED_R.log
cd ..
echo "TEST_MULTIALLELIC_PHASED_R passed."

cd TEST_MALFORMED_INPUTS
./run_tests.sh $d $2 $3 > TEST_MALFORMED_INPUTS.log
cd ..
echo "TEST_MALFORMED_INPUTS passed."

cd TEST_LGEN
./run_tests.sh $d $2 $3 > TEST_LGEN.log
cd ..
echo "TEST_LGEN passed."

cd TEST_VARIANT_JOIN
./run_tests.sh $d $2 $3 > TEST_VARIANT_JOIN.log
cd ..
echo "TEST_VARIANT_JOIN passed."

cd TEST_NORMALIZE_INDEL_JOIN
./run_tests.sh $d $2 $3 > TEST_NORMALIZE_INDEL_JOIN.log
cd ..
echo "TEST_NORMALIZE_INDEL_JOIN passed."

cd TEST_SPLIT_PAR_ALL_PAR1
./run_tests.sh $d $2 $3 > TEST_SPLIT_PAR_ALL_PAR1.log
cd ..
echo "TEST_SPLIT_PAR_ALL_PAR1 passed."

cd TEST_VALIDATE_DUP_ALLELE
./run_tests.sh $d $2 $3 > TEST_VALIDATE_DUP_ALLELE.log
cd ..
echo "TEST_VALIDATE_DUP_ALLELE passed."

cd TEST_ALLELE_UNIQUENESS
./run_tests.sh $d $2 $3 > TEST_ALLELE_UNIQUENESS.log
cd ..
echo "TEST_ALLELE_UNIQUENESS passed."

cd TEST_LEGACY_FLAG_ALIASES
./run_tests.sh $d $2 $3 > TEST_LEGACY_FLAG_ALIASES.log
cd ..
echo "TEST_LEGACY_FLAG_ALIASES passed."

cd TEST_KEEP_IF_BINARY
./run_tests.sh $d $2 $3 > TEST_KEEP_IF_BINARY.log
cd ..
echo "TEST_KEEP_IF_BINARY passed."

cd TEST_INFO_FILTER_STR
./run_tests.sh $d $2 $3 > TEST_INFO_FILTER_STR.log
cd ..
echo "TEST_INFO_FILTER_STR passed."

cd TEST_MIND_GENO_Y_NOSEX
./run_tests.sh $d $2 $3 > TEST_MIND_GENO_Y_NOSEX.log
cd ..
echo "TEST_MIND_GENO_Y_NOSEX passed."

cd TEST_SNPS_RANGE_DUP
./run_tests.sh $d $2 $3 > TEST_SNPS_RANGE_DUP.log
cd ..
echo "TEST_SNPS_RANGE_DUP passed."

cd TEST_Q_SCORE_RANGE_SPARSE
./run_tests.sh $d $2 $3 > TEST_Q_SCORE_RANGE_SPARSE.log
cd ..
echo "TEST_Q_SCORE_RANGE_SPARSE passed."

cd TEST_IMPORT_DOSAGE_CERTAINTY
./run_tests.sh $d $2 $3 > TEST_IMPORT_DOSAGE_CERTAINTY.log
cd ..
echo "TEST_IMPORT_DOSAGE_CERTAINTY passed."

cd TEST_PHENO_SVD
./run_tests.sh $d $2 $3 > TEST_PHENO_SVD.log
cd ..
echo "TEST_PHENO_SVD passed."

cd TEST_GLM_JOINT_DF
./run_tests.sh $d $2 $3 > TEST_GLM_JOINT_DF.log
cd ..
echo "TEST_GLM_JOINT_DF passed."

cd TEST_GLM_PERM
./run_tests.sh $d $2 $3 > TEST_GLM_PERM.log
cd ..
echo "TEST_GLM_PERM passed."

cd TEST_LD_CONSOLE_CHRX_FEMALES
./run_tests.sh $d $2 $3 > TEST_LD_CONSOLE_CHRX_FEMALES.log
cd ..
echo "TEST_LD_CONSOLE_CHRX_FEMALES passed."

cd TEST_SAMPLE_DIFF_ALL_MISSING
./run_tests.sh $d $2 $3 > TEST_SAMPLE_DIFF_ALL_MISSING.log
cd ..
echo "TEST_SAMPLE_DIFF_ALL_MISSING passed."

cd TEST_INDIV_SORT_PHASE_SEX
./run_tests.sh $d $2 $3 > TEST_INDIV_SORT_PHASE_SEX.log
cd ..
echo "TEST_INDIV_SORT_PHASE_SEX passed."

cd TEST_PROVREF_FLAGS_SUBSET
./run_tests.sh $d $2 $3 > TEST_PROVREF_FLAGS_SUBSET.log
cd ..
echo "TEST_PROVREF_FLAGS_SUBSET passed."

cd TEST_FLIP_SUBSET_NONSNP
./run_tests.sh $d $2 $3 > TEST_FLIP_SUBSET_NONSNP.log
cd ..
echo "TEST_FLIP_SUBSET_NONSNP passed."

cd TEST_ZERO_CLUSTER_MULTIALLELIC
./run_tests.sh $d $2 $3 > TEST_ZERO_CLUSTER_MULTIALLELIC.log
cd ..
echo "TEST_ZERO_CLUSTER_MULTIALLELIC passed."

cd TEST_PMERGE_EMPTY_FILESET
./run_tests.sh $d $2 $3 > TEST_PMERGE_EMPTY_FILESET.log
cd ..
echo "TEST_PMERGE_EMPTY_FILESET passed."

cd TEST_EXPORT_ID_DELIMS
./run_tests.sh $d $2 $3 > TEST_EXPORT_ID_DELIMS.log
cd ..
echo "TEST_EXPORT_ID_DELIMS passed."

cd TEST_EXPORT_OXFORD_NOSEX_Y
./run_tests.sh $d $2 $3 > TEST_EXPORT_OXFORD_NOSEX_Y.log
cd ..
echo "TEST_EXPORT_OXFORD_NOSEX_Y passed."

cd TEST_GLM_MULTIALLELIC_CONST_ALLELE
./run_tests.sh $d $2 $3 > TEST_GLM_MULTIALLELIC_CONST_ALLELE.log
cd ..
echo "TEST_GLM_MULTIALLELIC_CONST_ALLELE passed."

cd TEST_VCF_IDSPACE_TO
./run_tests.sh $d $2 $3 > TEST_VCF_IDSPACE_TO.log
cd ..
echo "TEST_VCF_IDSPACE_TO passed."

cd TEST_ADJUST_LOG10_MINP
./run_tests.sh $d $2 $3 > TEST_ADJUST_LOG10_MINP.log
cd ..
echo "TEST_ADJUST_LOG10_MINP passed."

cd TEST_SET_MISSING_DOSAGE
./run_tests.sh $d $2 $3 > TEST_SET_MISSING_DOSAGE.log
cd ..
echo "TEST_SET_MISSING_DOSAGE passed."

cd TEST_GLM_LOCAL_COVAR
./run_tests.sh $d $2 $3 > TEST_GLM_LOCAL_COVAR.log
cd ..
echo "TEST_GLM_LOCAL_COVAR passed."

cd TEST_GLM_CHRY_COVAR_SUBSETS
./run_tests.sh $d $2 $3 > TEST_GLM_CHRY_COVAR_SUBSETS.log
cd ..
echo "TEST_GLM_CHRY_COVAR_SUBSETS passed."

cd TEST_KING_SPARSE_NSNP
./run_tests.sh $d $2 $3 > TEST_KING_SPARSE_NSNP.log
cd ..
echo "TEST_KING_SPARSE_NSNP passed."

cd TEST_CLUMP_FORCE_A1_MULTIALLELIC
./run_tests.sh $d $2 $3 > TEST_CLUMP_FORCE_A1_MULTIALLELIC.log
cd ..
echo "TEST_CLUMP_FORCE_A1_MULTIALLELIC passed."

cd TEST_SCORE_CHRY_NOSEX_XCHR0
./run_tests.sh $d $2 $3 > TEST_SCORE_CHRY_NOSEX_XCHR0.log
cd ..
echo "TEST_SCORE_CHRY_NOSEX_XCHR0 passed."

cd TEST_FST_JACKKNIFE_FULL_BLOCK
./run_tests.sh $d $2 $3 > TEST_FST_JACKKNIFE_FULL_BLOCK.log
cd ..
echo "TEST_FST_JACKKNIFE_FULL_BLOCK passed."

cd TEST_GWAS_SSF_RSID_RETRY
./run_tests.sh $d $2 $3 > TEST_GWAS_SSF_RSID_RETRY.log
cd ..
echo "TEST_GWAS_SSF_RSID_RETRY passed."

cd TEST_PMERGE_PSAM_CONFLICTS
./run_tests.sh $d $2 $3 > TEST_PMERGE_PSAM_CONFLICTS.log
cd ..
echo "TEST_PMERGE_PSAM_CONFLICTS passed."

cd TEST_NORMALIZE_RTRIM
./run_tests.sh $d $2 $3 > TEST_NORMALIZE_RTRIM.log
cd ..
echo "TEST_NORMALIZE_RTRIM passed."

cd TEST_FA_CONTIG_LENGTHS
./run_tests.sh $d $2 $3 > TEST_FA_CONTIG_LENGTHS.log
cd ..
echo "TEST_FA_CONTIG_LENGTHS passed."

cd TEST_SPLIT_CAT_PHENO_OMIT
./run_tests.sh $d $2 $3 > TEST_SPLIT_CAT_PHENO_OMIT.log
cd ..
echo "TEST_SPLIT_CAT_PHENO_OMIT passed."

cd TEST_PHENO_COL_NUMS_WIDE
./run_tests.sh $d $2 $3 > TEST_PHENO_COL_NUMS_WIDE.log
cd ..
echo "TEST_PHENO_COL_NUMS_WIDE passed."

cd TEST_RECOVER_VAR_IDS_ROUNDTRIP
./run_tests.sh $d $2 $3 > TEST_RECOVER_VAR_IDS_ROUNDTRIP.log
cd ..
echo "TEST_RECOVER_VAR_IDS_ROUNDTRIP passed."

cd TEST_GENO_DOSAGE_CHRX_SUBSET
./run_tests.sh $d $2 $3 > TEST_GENO_DOSAGE_CHRX_SUBSET.log
cd ..
echo "TEST_GENO_DOSAGE_CHRX_SUBSET passed."

cd TEST_HET_SMALL_SAMPLE
./run_tests.sh $d $2 $3 > TEST_HET_SMALL_SAMPLE.log
cd ..
echo "TEST_HET_SMALL_SAMPLE passed."

cd TEST_MENDEL_CHRY_NOBS
./run_tests.sh $d $2 $3 > TEST_MENDEL_CHRY_NOBS.log
cd ..
echo "TEST_MENDEL_CHRY_NOBS passed."

cd TEST_BGEN13_SHORT_CHR_THREADS
./run_tests.sh $d $2 $3 > TEST_BGEN13_SHORT_CHR_THREADS.log
cd ..
echo "TEST_BGEN13_SHORT_CHR_THREADS passed."

cd TEST_FILTER_FLAGS
./run_tests.sh $d $2 $3 > TEST_FILTER_FLAGS.log
cd ..
echo "TEST_FILTER_FLAGS passed."

cd TEST_VCF_IMPORT_FLAGS
./run_tests.sh $d $2 $3 > TEST_VCF_IMPORT_FLAGS.log
cd ..
echo "TEST_VCF_IMPORT_FLAGS passed."

cd TEST_PMERGE_MODES
./run_tests.sh $d $2 $3 > TEST_PMERGE_MODES.log
cd ..
echo "TEST_PMERGE_MODES passed."

cd TEST_COL_MATCH_FILTERS
./run_tests.sh $d $2 $3 > TEST_COL_MATCH_FILTERS.log
cd ..
echo "TEST_COL_MATCH_FILTERS passed."

cd TEST_SPECIES_FLAGS
./run_tests.sh $d $2 $3 > TEST_SPECIES_FLAGS.log
cd ..
echo "TEST_SPECIES_FLAGS passed."

cd TEST_KING_MATRIX
./run_tests.sh $d $2 $3 > TEST_KING_MATRIX.log
cd ..
echo "TEST_KING_MATRIX passed."

cd TEST_PHENO_TRANSFORMS
./run_tests.sh $d $2 $3 > TEST_PHENO_TRANSFORMS.log
cd ..
echo "TEST_PHENO_TRANSFORMS passed."

cd TEST_ADJUST_FILE_FIELDS
./run_tests.sh $d $2 $3 > TEST_ADJUST_FILE_FIELDS.log
cd ..
echo "TEST_ADJUST_FILE_FIELDS passed."

cd TEST_CLUMP_FIELDS
./run_tests.sh $d $2 $3 > TEST_CLUMP_FIELDS.log
cd ..
echo "TEST_CLUMP_FIELDS passed."

cd TEST_PGEN_INFO_COUNTS
./run_tests.sh $d $2 $3 > TEST_PGEN_INFO_COUNTS.log
cd ..
echo "TEST_PGEN_INFO_COUNTS passed."

cd TEST_GLM_OUTPUT_FLAGS
./run_tests.sh $d $2 $3 > TEST_GLM_OUTPUT_FLAGS.log
cd ..
echo "TEST_GLM_OUTPUT_FLAGS passed."

cd TEST_SAMPLE_ID_PHENO_FLAGS
./run_tests.sh $d $2 $3 > TEST_SAMPLE_ID_PHENO_FLAGS.log
cd ..
echo "TEST_SAMPLE_ID_PHENO_FLAGS passed."

cd TEST_VARIANT_FILTER_OPTS
./run_tests.sh $d $2 $3 > TEST_VARIANT_FILTER_OPTS.log
cd ..
echo "TEST_VARIANT_FILTER_OPTS passed."

cd TEST_LD_OPTIONS
./run_tests.sh $d $2 $3 > TEST_LD_OPTIONS.log
cd ..
echo "TEST_LD_OPTIONS passed."

cd TEST_IMPORT_OUTPUT_OPTS
./run_tests.sh $d $2 $3 > TEST_IMPORT_OUTPUT_OPTS.log
cd ..
echo "TEST_IMPORT_OUTPUT_OPTS passed."

echo "All tests passed."
