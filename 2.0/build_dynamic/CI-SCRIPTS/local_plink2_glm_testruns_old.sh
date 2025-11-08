#!/bin/bash
set -e

# =================== CONFIG ===================
datapath="test_data/"
sbpath=""         # optional subdirectory
outdir="./results"
threads=4
mkdir -p "$outdir"

# ------------------- Datasets -------------------
froot_=(
  "1kgp3_50k_nomiss_Av_nonintdose"
  "1kgp3_50k_yesmiss_Av_nonintdose"
)

# ------------------- Covariates -------------------
cov_options=( "COV_1" "" )  # empty string means no covariate

# ------------------- Subset sizes -------------------
thin=(1000 32000 32017)

# =================== FUNCTIONS ===================

run_plink() {
    local fname="$1"
    local model="$2"   # "linear", "logistic", "firth"
    local n="$3"
    local cov="$4"

    # Determine phenotype based on model
    local phenotype=""
    if [[ "$model" == "linear" ]]; then
        phenotype="y"
    else
        phenotype="ybool"
    fi

    local pfile="${fname}_recode_varIDs"
    local phenofile="${fname}_combined_phenocov.csv"
    local keep_file="${datapath}${fname}_recode_varIDs_subset_${n}.keep"

    # Covariate args
    local cov_arg=""
    local allow_no_covars=""
    if [[ -n "$cov" ]]; then
        cov_arg="--covar ${datapath}${phenofile} --covar-name $cov"
    else
        allow_no_covars="allow-no-covars"
    fi

    # Determine PLINK2 --glm flags
    local glm_flags="hide-covar"
    if [[ "$model" == "logistic" ]]; then
        glm_flags="no-firth hide-covar"
    fi

    # Construct output prefix
    local cov_suffix=$( [[ -z "$cov" ]] && echo "noCov" || echo "$cov" )
    local outprefix="${outdir}/${fname}_${phenotype}_${cov_suffix}_glm${model:+_$model}_keep${n}"

    echo "📌 Running PLINK2 GLM: fname=$fname, model=$model, phenotype=$phenotype, n=$n, cov=${cov:-none}"

    plink2 --glm $glm_flags \
        $allow_no_covars \
        --pfile "${datapath}${pfile}" \
        --allow-extra-chr \
        --pheno "${datapath}${phenofile}" \
        --pheno-name "$phenotype" \
        --keep "$keep_file" \
        --threads "$threads" \
        $cov_arg \
        --out "$outprefix"

    # Run comparison to R reference if applicable
    compare_to_r "$phenotype" "$fname" "$n" "$cov_suffix" "$model" "$outprefix"

    echo "✅ Completed: $outprefix"
}

compare_to_r() {
    local phenotype="$1"
    local fname="$2"
    local n="$3"
    local cov_suffix="$4"
    local model="$5"
    local outprefix="$6"

    local glm_file=""
    local ref_file=""

    case "$model" in
        linear)
            glm_file="${outprefix}.${phenotype}.glm.linear"
            ref_file="${datapath}${fname}_${phenotype}_${cov_suffix}_glm_linear_keep${n}_linear.csv"
            ;;
        logistic)
            glm_file="${outprefix}.${phenotype}.glm.logistic"
            ref_file="${datapath}${fname}_${phenotype}_${cov_suffix}_glm_logistic_keep${n}_logistic.csv"
            ;;
        firth)
            glm_file="${outprefix}.${phenotype}.glm.logistic.hybrid"
            ref_file="${datapath}${fname}_${phenotype}_${cov_suffix}_glm_firth_keep${n}_firth.csv"
            ;;
        *)
            echo "ℹ️  Skipping comparison for unknown model $model"
            return
            ;;
    esac

    if [[ -f "$glm_file" && -f "$ref_file" ]]; then
        echo "🔍 Comparing PLINK2 GLM output to R reference"
        echo "   PLINK2 file: $glm_file"
        echo "   R reference: $ref_file"
        python3 2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py "$glm_file" "$ref_file"
    else
        echo "⚠️  Skipping comparison: one or both files do not exist"
        echo "   PLINK2 file: $glm_file"
        echo "   R reference: $ref_file"
    fi
}

# =================== SINGLE-PARAMETER TEST ===================
echo "=================== SINGLE-PARAMETER TEST ==================="
test_fname="1kgp3_50k_nomiss_Av_nonintdose"
test_model="firth"   # <-- specify the model here: linear, logistic, firth
test_n=32000
test_cov="COV_1"          # "" = no covariate, or "COV_1"

run_plink "$test_fname" "$test_model" "$test_n" "$test_cov"
echo "=================== SINGLE-PARAMETER TEST COMPLETE ==================="

# =================== FULL LOOP ===================
: <<'COMMENT'
for fname in "${froot_[@]}"; do
    for n in "${thin[@]}"; do
        for cov in "${cov_options[@]}"; do
            # Specify the model for this run
            for model in "linear" "logistic" "firth"; do
                run_plink "$fname" "$model" "$n" "$cov"
            done
        done
    done
done
COMMENT
