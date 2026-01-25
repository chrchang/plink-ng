#!/bin/bash
# Don't use set -e - we want to control error handling explicitly
set -o pipefail

# Force unbuffered output for GitHub Actions
export PYTHONUNBUFFERED=1

# Error trap to show where failures occur
trap 'echo "❌ ERROR at line $LINENO in ${BASH_SOURCE[0]}" >&2; exit 1' ERR

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

    # Check if input files exist before running PLINK2
    if [[ ! -f "${datapath}${pfile}.pgen" ]]; then
        echo -e "${RED}ERROR: PLINK2 pgen file not found: ${datapath}${pfile}.pgen${RESET}" >&2
        exit 1
    fi
    if [[ ! -f "${datapath}${phenofile}" ]]; then
        echo -e "${RED}ERROR: Phenotype file not found: ${datapath}${phenofile}${RESET}" >&2
        exit 1
    fi
    if [[ ! -f "$keep_file" ]]; then
        echo -e "${RED}ERROR: Keep file not found: $keep_file${RESET}" >&2
        exit 1
    fi

    # Run PLINK2 and capture errors
    if ! plink2 --glm $glm_flags \
        $allow_no_covars \
        --pfile "${datapath}${pfile}" \
        --allow-extra-chr \
        --pheno "${datapath}${phenofile}" \
        --pheno-name "$phenotype" \
        --keep "$keep_file" \
        --threads "$threads" \
        $cov_arg \
        --out "$outprefix"; then
        echo -e "${RED}ERROR: PLINK2 command failed with exit code $?${RESET}" >&2
        exit 1
    fi

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

        echo -e "\n${BOLD}${BLUE}╔════════════════════════════════════════════════════════════╗${RESET}"
        echo -e "${BOLD}${BLUE}║           CORRELATION RESULTS vs THRESHOLD                 ║${RESET}"
        echo -e "${BOLD}${BLUE}║           Required: ≥ ${correlation_threshold}                                  ║${RESET}"
        echo -e "${BOLD}${BLUE}╚════════════════════════════════════════════════════════════╝${RESET}"
        
        # Check for correlations below threshold
        correlation_failed=false
        while IFS= read -r line; do
            # Extract correlation values from lines containing numbers
            if [[ "$line" =~ ([0-9]+\.[0-9]+) ]]; then
                corr_value="${BASH_REMATCH[1]}"
                
                # Use bc for floating point comparison
                if (( $(echo "$corr_value < $correlation_threshold" | bc -l) )); then
                    correlation_failed=true
                    echo -e "${RED}  ❌ FAIL: $line${RESET}"
                    echo -e "${RED}     └─ Correlation $corr_value < $correlation_threshold (BELOW THRESHOLD)${RESET}"
                else
                    echo -e "${GREEN}  ✓ PASS: $line${RESET}"
                    echo -e "${GREEN}     └─ Correlation $corr_value ≥ $correlation_threshold${RESET}"
                fi
            else
                # Print non-numeric lines normally
                echo -e "${YELLOW}  $line${RESET}"
            fi
        done <<< "$comparison_output"
        
        echo -e "${BOLD}${BLUE}════════════════════════════════════════════════════════════${RESET}\n"
        
        # Exit if correlation threshold not met
        if [[ "$correlation_failed" == "true" ]]; then
            echo -e "\n${RED}╔════════════════════════════════════════════════════════════╗${RESET}" >&2
            echo -e "${RED}║                    ❌ TEST FAILED ❌                       ║${RESET}" >&2
            echo -e "${RED}║  One or more correlations below ${correlation_threshold} threshold          ║${RESET}" >&2
            echo -e "${RED}╚════════════════════════════════════════════════════════════╝${RESET}" >&2
            echo -e "${RED}${BOLD}Failed Test Details:${RESET}" >&2
            echo -e "${RED}  Test:       $test_counter of $total_tests${RESET}" >&2
            echo -e "${RED}  Dataset:    $fname${RESET}" >&2
            echo -e "${RED}  Model:      $model${RESET}" >&2
            echo -e "${RED}  Subset:     $n samples${RESET}" >&2
            echo -e "${RED}  Covariate:  $cov_suffix${RESET}" >&2
            echo -e "${RED}  Missing:    $missing_geno${RESET}" >&2
            echo -e "" >&2
            echo -e "${YELLOW}${BOLD}Progress Before Failure:${RESET}" >&2
            echo -e "${YELLOW}  Completed: $((test_counter - 1)) of $total_tests tests${RESET}" >&2
            echo -e "${YELLOW}  Remaining: $((total_tests - test_counter)) tests (not run)${RESET}" >&2
            echo -e "" >&2
            echo -e "${RED}${BOLD}Required Action: Fix correlation issues before proceeding${RESET}\n" >&2
            exit 1
        else
            echo -e "${GREEN}${BOLD}✅ ALL CORRELATIONS PASSED (≥ $correlation_threshold)${RESET}\n"
        fi

    else
        echo -e "${RED}⚠️  Skipping comparison: one or both files missing.${RESET}" >&2
        [[ ! -f "$glm_file" ]] && echo -e "   Missing: ${glm_file}" >&2
        [[ ! -f "$ref_file" ]] && echo -e "   Missing: ${ref_file}" >&2
        echo -e "" >&2
        echo -e "${YELLOW}${BOLD}Progress Before Failure:${RESET}" >&2
        echo -e "${YELLOW}  Completed: $((test_counter - 1)) of $total_tests tests${RESET}" >&2
        echo -e "${YELLOW}  Remaining: $((total_tests - test_counter)) tests (not run)${RESET}" >&2
        echo -e "" >&2
        echo -e "${RED}This is a FAILURE condition - missing reference files${RESET}" >&2
        exit 1
    fi
}

# ==========================================================
# ====================== MAIN LOOP =========================
# ==========================================================
echo -e "${BLUE}=================== PLINK2 vs R TESTS ===================${RESET}"
echo -e "${BOLD}Batch:${RESET} $BATCH_NUM of 12"
echo -e "${BOLD}Number of tests:${RESET} ${#params[@]}"
echo -e "${BOLD}Correlation threshold:${RESET} $correlation_threshold"
echo -e "${BLUE}=========================================================${RESET}\n"

# Start timer
start_time=$(date +%s)

# ------------------------------
# Initialize test tracking
# ------------------------------
test_counter=0
declare -a passed_tests
declare -a failed_tests
total_tests=${#params[@]}

echo -e "${YELLOW}DEBUG: params array has ${#params[@]} elements${RESET}"
echo -e "${YELLOW}DEBUG: First param: ${params[0]}${RESET}"
echo -e "${YELLOW}DEBUG: All params:${RESET}"
for p in "${params[@]}"; do
    echo -e "${YELLOW}  - $p${RESET}"
done
echo ""

# ------------------------------
# Loop through parameter list from config
# ------------------------------
for param in "${params[@]}"; do
   [[ -z "$param" || "$param" =~ ^# ]] && continue
   
   ((test_counter++))
   
   echo -e "\n${GREEN}═══════════════════════════════════════════════════════════${RESET}"
   echo -e "${BOLD}${GREEN}TEST $test_counter of $total_tests${RESET}"
   echo -e "${GREEN}═══════════════════════════════════════════════════════════${RESET}"
   echo -e "${YELLOW}DEBUG: Parsing parameter line: $param${RESET}"
   
   IFS=',' read -r fname model n cov threads <<< "$param"
   
   echo -e "${YELLOW}DEBUG: Parsed values:${RESET}"
   echo -e "${YELLOW}  fname=$fname${RESET}"
   echo -e "${YELLOW}  model=$model${RESET}"
   echo -e "${YELLOW}  n=$n${RESET}"
   echo -e "${YELLOW}  cov=$cov${RESET}"
   echo -e "${YELLOW}  threads=$threads${RESET}"
   
   # Track current test info
   current_test="Test $test_counter: $fname, $model, n=$n, cov=${cov:-none}, threads=$threads"
   
   echo -e "${YELLOW}DEBUG: About to call run_plink...${RESET}"
   
   # Run the test (will exit if correlation fails)
   run_plink "$fname" "$model" "$n" "$cov" "$threads"
   
   echo -e "${YELLOW}DEBUG: run_plink completed successfully${RESET}"
   
   # If we reach here, test passed
   passed_tests+=("$current_test")
done

echo -e "\n${GREEN}╔════════════════════════════════════════════════════════════╗${RESET}"
echo -e "${GREEN}║                  🎉 ALL TESTS PASSED 🎉                    ║${RESET}"
echo -e "${GREEN}╚════════════════════════════════════════════════════════════╝${RESET}"

# Calculate elapsed time
end_time=$(date +%s)
elapsed=$((end_time - start_time))
elapsed_min=$((elapsed / 60))
elapsed_sec=$((elapsed % 60))

echo -e "${BOLD}${GREEN}Summary:${RESET}"
echo -e "${GREEN}  ✓ Batch:       $BATCH_NUM of 12${RESET}"
echo -e "${GREEN}  ✓ Completed:   $test_counter / $total_tests tests${RESET}"
echo -e "${GREEN}  ✓ Threshold:   All correlations ≥ $correlation_threshold${RESET}"
echo -e "${GREEN}  ✓ Time:        ${elapsed_min}m ${elapsed_sec}s${RESET}"

echo -e "\n${BOLD}${GREEN}Parameters Tested:${RESET}"
for test_info in "${passed_tests[@]}"; do
    echo -e "${GREEN}  • ${test_info}${RESET}"
done

echo -e "${GREEN}=========================================================${RESET}\n"