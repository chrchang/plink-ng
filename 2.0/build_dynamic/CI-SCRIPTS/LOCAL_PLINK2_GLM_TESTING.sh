#!/bin/bash
set -e

# ==========================================================
#  PLINK2 vs R GLM Comparison Pipeline
# ==========================================================
# Loads parameters from PLINK2_GLM_TEST_CONFIG.sh
# and executes corresponding runs.
# ==========================================================

# Load configuration
CONFIG_FILE="./2.0/build_dynamic/CI-SCRIPTS/PLINK2_GLM_TEST_CONFIG.sh"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "❌ Config file not found: $CONFIG_FILE"
  exit 1
fi
source "$CONFIG_FILE"

mkdir -p "$outdir"

# ---------------------------
# Color helpers
# ---------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
BOLD='\033[1m'
RESET='\033[0m'

# -------------------
# Helper: Run PLINK2
# -------------------
run_plink() {
    local fname="$1"
    local model="$2"
    local n="$3"
    local cov="$4"
    local threads="$5"

    # Phenotype determined by model
    local phenotype=""
    if [[ "$model" == "linear" ]]; then
        phenotype="y"
    else
        phenotype="ybool"
    fi

    # Determine missing genotype type from fname
    local missing_geno="Unknown"
    if [[ "$fname" == *"yesmiss"* ]]; then
        missing_geno="Yes"
    elif [[ "$fname" == *"nomiss"* ]]; then
        missing_geno="No"
    fi

    local pfile="${fname}_recode_varIDs"
    local phenofile="${fname}_combined_phenocov.csv"
    local keep_file="${datapath}${fname}_recode_varIDs_subset_${n}.keep"

    local cov_arg=""
    local allow_no_covars=""
    if [[ -n "$cov" ]]; then
        cov_arg="--covar ${datapath}${phenofile} --covar-name $cov"
    else
        allow_no_covars="allow-no-covars"
    fi

    local glm_flags="hide-covar"
    if [[ "$model" == "logistic" ]]; then
        glm_flags="no-firth hide-covar"
    fi

    local cov_suffix
    if [[ -z "$cov" ]]; then
        cov_suffix="noCov"
    else
        cov_suffix="$cov"
    fi

    local outprefix="${outdir}/${fname}_${phenotype}_${cov_suffix}_glm_${model}_keep${n}"

    echo -e "${BLUE}---------------------------------------------------${RESET}"
    echo -e "${BOLD}📌 Running PLINK2 GLM${RESET}"
    echo "   Dataset:          $fname"
    echo "   Missing genotypes:${missing_geno}"
    echo "   Model:            $model"
    echo "   Phenotype:        $phenotype"
    echo "   Subset size:      $n"
    echo "   Covariate:        ${cov:-none}"
    echo "   Threads:          $threads"
    echo "   Output:           $outprefix"
    echo -e "${BLUE}---------------------------------------------------${RESET}"

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

    compare_to_r "$phenotype" "$fname" "$n" "$cov_suffix" "$model" "$outprefix" "$missing_geno"
}

# -------------------
# Helper: Compare PLINK2 vs R Results
# -------------------
compare_to_r() {
    local phenotype="$1"
    local fname="$2"
    local n="$3"
    local cov_suffix="$4"
    local model="$5"
    local outprefix="$6"
    local missing_geno="$7"

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
    esac

    echo -e "\n${YELLOW}============================================================${RESET}"
    echo -e "${BOLD}🔎  CORRELATION CHECK — PLINK2 vs R RESULTS${RESET}"
    echo -e "${YELLOW}============================================================${RESET}"
    echo -e " ${BOLD}Phenotype:${RESET}         $phenotype"
    echo -e " ${BOLD}Model:${RESET}             $model"
    echo -e " ${BOLD}Subset size:${RESET}       $n"
    echo -e " ${BOLD}Covariate:${RESET}        $cov_suffix"
    echo -e " ${BOLD}Missing genotypes:${RESET} $missing_geno"
    echo -e " ${BOLD}PLINK2 file:${RESET}      $glm_file"
    echo -e " ${BOLD}R ref file:${RESET}       $ref_file"

    if [[ -f "$glm_file" && -f "$ref_file" ]]; then
        echo -e "${GREEN}✅ Files found. Running Python comparison...${RESET}"

        # Capture Python output
        comparison_output=$(python3 "$compare_script" "$glm_file" "$ref_file" 2>&1)

        echo -e "\n${BOLD}${YELLOW}================== CORRELATION RESULTS ==================${RESET}"
        # Highlight numeric values (correlations) in yellow, rest green
        echo "$comparison_output" | while IFS= read -r line; do
            if [[ "$line" =~ [0-9]+\.[0-9]+ ]]; then
                echo -e "${YELLOW}${line}${RESET}"
            else
                echo -e "${GREEN}${line}${RESET}"
            fi
        done
        echo -e "${BOLD}${YELLOW}=========================================================${RESET}\n"

    else
        echo -e "${RED}⚠️  Skipping comparison: one or both files missing.${RESET}"
        [[ ! -f "$glm_file" ]] && echo -e "   Missing: ${glm_file}"
        [[ ! -f "$ref_file" ]] && echo -e "   Missing: ${ref_file}"
    fi
}

# ==========================================================
# ====================== MAIN LOOP =========================
# ==========================================================
echo -e "${BLUE}=================== PLINK2 vs R TESTS ===================${RESET}"

# ------------------------------
# Single permutation test example
# ------------------------------
# Example: Uncomment to run a single test
# run_plink "1kgp3_50k_nomiss_Av_nonintdose" "firth" 1000 "COV_1" 4

# ------------------------------
# Loop through full parameter list
# ------------------------------
for param in "${params[@]}"; do
   [[ -z "$param" || "$param" =~ ^# ]] && continue
   IFS=',' read -r fname model n cov threads <<< "$param"
   run_plink "$fname" "$model" "$n" "$cov" "$threads"
done

echo -e "${GREEN}=================== ALL TESTS COMPLETE ===================${RESET}"
