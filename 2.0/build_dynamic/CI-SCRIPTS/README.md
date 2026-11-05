# PLINK2 GLM Testing Framework

Continuous integration framework for validating PLINK2 GLM implementations against R reference outputs.

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Test Data](#test-data)
  - [Generation (PLINK2SYN)](#generation-plink2syn)
  - [Caching Workflow](#caching-workflow)
  - [Expected Files](#expected-files)
- [CI/CD Pipeline](#cicd-pipeline)
  - [Automated Testing](#automated-testing)
  - [Test Coverage](#test-coverage)
- [Configuration](#configuration)
  - [Batch System](#batch-system)
  - [Parameter Format](#parameter-format)
  - [Test Conditions](#test-conditions)
  - [Config File](#config-file)
- [Scripts Reference](#scripts-reference)
- [Directory Structure](#directory-structure)
- [Troubleshooting](#troubleshooting)

## Overview

**Purpose:** Ensure PLINK2 GLM statistical accuracy by comparing results with ground-truth R outputs on synthetic datasets.

**Test Coverage:** 108 test configurations across 12 batches (9 tests each)

**CI Schedule:** Weekly via [.github/workflows/auto_checkcommits_tests_linuxmac.yml](/.github/workflows/auto_checkcommits_tests_linuxmac.yml) on Linux/macOS

**Complete Cycle:** ~12 weeks to run all batches

**Validation Criteria:** Correlations (beta, p-values) must be ≥0.8, and median differences (relative for beta and SE, absolute for p-values) must be ≤5%

## Quick Start

```bash
# Run default batch (Batch 1)
./2.0/build_dynamic/CI-SCRIPTS/GLM_TESTS.sh

# Run specific batch (1-12)
export GLM_BATCH_NUM=5
./2.0/build_dynamic/CI-SCRIPTS/GLM_TESTS.sh

# Run all batches sequentially
for i in {1..12}; do
  export GLM_BATCH_NUM=$i
  ./2.0/build_dynamic/CI-SCRIPTS/GLM_TESTS.sh
done
```

**Prerequisites:**
- Test data in `test_data/` directory (see [Test Data](#test-data))
- Python 3.x with pandas and numpy (`pip install pandas numpy`)
- System dependencies: LAPACK, OpenBLAS (Linux), zlib

## How It Works

```
Test Data (PLINK2SYN) → Cache (Google Drive) → CI Workflow → PLINK2 GLM → Validate vs R → Pass/Fail
```

### Pipeline Steps

1. **Test Generation:** Synthetic datasets created via [PLINK2SYN](https://github.com/physioforecast/PLINK2SYN) with known ground truth
   - Based on 1000 Genomes Project Phase 3 data
   - ~32k synthetic samples with controlled genetic effects
   - Both quantitative and binary phenotypes
   - R reference outputs computed for validation

2. **Data Caching:** Test data stored in Google Drive, cached by [.github/workflows/cache_test_data.yml](/.github/workflows/cache_test_data.yml)
   - Downloads once, reused across all workflow runs
   - Cache key versioning for updates

3. **CI Execution:** Weekly workflow builds PLINK2, runs GLM tests, validates against R reference outputs
   - Checks for new commits before running
   - Rotates through batches 1-12 sequentially
   - Runs on both Linux and macOS

4. **Validation:** Correlation check (≥0.8) for beta coefficients and p-values, and median difference check (≤5%) for beta coefficients, standard errors, and p-values
   - Compares PLINK2 output with R reference
   - Fails fast if correlation below threshold
   - Detailed output in `results/` directory

## Test Data

### Generation (PLINK2SYN)

Test datasets generated via [PLINK2SYN](https://github.com/physioforecast/PLINK2SYN), a genomics research pipeline for creating synthetic genotype-phenotype datasets.

**PLINK2SYN Pipeline:**
- **Purpose:** Generates synthetic GWAS datasets with known ground truth for validating PLINK2 statistical methods
- **Source Data:** 1000 Genomes Project Phase 3 (~50,000 biallelic SNPs)
- **Synthetic Samples:** ~32,000 augmented samples created by adding controlled Gaussian noise (N(0, 0.1)) to real genotypes
- **Phenotypes:**
  - Quantitative (`y`) - continuous phenotype
  - Binary (`ybool`) - case/control phenotype
  - Configurable heritability with ~500 causal variants (selected at probability 0.01)
- **Covariates:** Optional covariates (e.g., `COV_1`) for testing covariate adjustment
- **Reference Outputs:** Ground truth R GLM results for validation against PLINK2
- **Scenarios:**
  - Datasets with missing genotypes (`yesmiss`)
  - Datasets without missing genotypes (`nomiss`)
  - Multiple sample sizes (1000, 32000, 32017)
  - Single-threaded and multi-threaded execution

This pipeline ensures PLINK2 GLM implementations can be validated against known statistical outcomes from R, providing confidence in the accuracy of GWAS results.

### Caching Workflow

Test data is managed via [.github/workflows/cache_test_data.yml](/.github/workflows/cache_test_data.yml):

**Storage & Caching:**
- **Storage:** Google Drive (File ID: `1lW5s8kwMeOwHPrbRoREtO5pwgKKK_8FR`)
- **Cache Key:** `test-data-v1` - shared across all workflow runs
- **Trigger:** Manual (workflow_dispatch) via GitHub Actions UI
- **Process:**
  1. Check if cache already exists (skip download if cached)
  2. Download `test_data.zip` from Google Drive using `gdown`
  3. Extract to `test_data/` directory
  4. Save to GitHub Actions cache
  5. Subsequent workflow runs restore from cache (no re-download)

**When to Run:**
- First-time setup of CI test data cache
- After updating test datasets (increment cache key to `test-data-v2`, etc.)
- If cache becomes corrupted or invalidated
- When test data version needs to be synchronized across workflows

**Updating Test Data:**
1. Regenerate test data using [PLINK2SYN](https://github.com/physioforecast/PLINK2SYN) pipeline with desired parameters
2. Create `test_data.zip` archive with all required files
3. Upload to Google Drive and update file ID in workflow
4. Increment cache key in both workflows (e.g., `test-data-v1` → `test-data-v2`)
5. Run cache_test_data.yml workflow to populate new cache

### Expected Files

Located in `test_data/` directory:

**Genotype Files (PLINK2 format):**
- `{dataset}_recode_varIDs.pgen` - PLINK2 genotype file (binary)
- `{dataset}_recode_varIDs.pvar` - Variant information (chromosome, position, alleles)
- `{dataset}_recode_varIDs.psam` - Sample information (IDs)

**Phenotype & Covariate Files:**
- `{dataset}_combined_phenocov.csv` - Phenotypes (`y`, `ybool`) and covariates (`COV_1`, etc.)

**Sample Subset Files:**
- `{dataset}_recode_varIDs_subset_{n}.keep` - Sample IDs for subsets
  - n=1000 (small subset)
  - n=32000 (large subset, 64×500)
  - n=32017 (non-64-multiple subset, 64×500+17)

**R Reference Files (Ground Truth):**
- `{dataset}_{phenotype}_{covariate}_glm_{model}_keep{n}_{model}.csv`
- Contains beta coefficients, standard errors, p-values from R GLM
- Used for correlation-based validation of PLINK2 results

**Example Datasets:**
- `1kgp3_50k_nomiss_Av_nonintdose` - Dataset without missing genotypes
- `1kgp3_50k_yesmiss_Av_nonintdose` - Dataset with missing genotypes

**Local Testing:**
For local development, contact maintainer for test data access or download from Google Drive using the file ID.

## CI/CD Pipeline

### Automated Testing

[.github/workflows/auto_checkcommits_tests_linuxmac.yml](/.github/workflows/auto_checkcommits_tests_linuxmac.yml)

**Trigger Schedule:**
- **Automatic:** Every Sunday at midnight (UTC) via cron schedule
- **Manual:** Via workflow_dispatch in GitHub Actions UI

**Execution Conditions:**
- Only runs if new commits exist on master branch since last successful workflow run
- Skips execution if no changes detected (saves CI resources)

**Platform Matrix:**
- `ubuntu-latest` - Tests on Linux
- `macos-latest` - Tests on macOS
- Both platforms run in parallel

**Workflow Steps:**

1. **Check Commits**
   - Compare latest master commit with last successful workflow commit
   - Output: `new_commits=true/false`
   - Skip remaining steps if no new commits

2. **Setup Environment**
   - Checkout repository with full git history
   - Install system dependencies:
     - **Linux:** libopenblas-dev, liblapack-dev, liblapacke-dev, build-essential, zlib1g-dev, unzip, bc
     - **macOS:** Uses native system libraries

3. **Determine Batch Number**
   - Derived from the date: `BATCH = (ISO week % 12) + 1`
   - Cycles through batches 1→2→...→12→1
   - Set as environment variable `GLM_BATCH_NUM`

4. **Restore Test Data**
   - Restore from GitHub Actions cache (key: `test-data-v1`)
   - Populated by [cache_test_data.yml](/.github/workflows/cache_test_data.yml) workflow
   - No download needed if cache hit

5. **Build PLINK2**
   - Execute `DYNAMICBUILD_LINUXMAC_ROOT.sh`
   - Compiles PLINK2 from source with system libraries
   - Ensures latest code changes are tested

6. **Run GLM Tests**
   - Execute `GLM_TESTS.sh` with current batch number
   - Sets up Python venv, installs dependencies
   - Runs `LOCAL_PLINK2_GLM_TESTING.sh`
   - Tests all configurations in current batch (~9 tests)
   - Validates each result against R reference
   - **Fails fast** if any correlation below threshold

7. **Cleanup**
   - Execute `CLEANUP.sh`
   - Removes test data, results, build artifacts
   - Saves disk space on CI runners

**Batch Rotation Logic:**
- Keeps no state, so a failed run or an evicted cache cannot skip a batch
- Ensures all 108 tests run over ~12 weeks
- Linux and macOS run the same batch

### Test Coverage

**Weekly Execution:**
- 1 batch × 2 platforms = 18 test executions per week
- ~9 tests per batch × 2 platforms = ~18 PLINK2 GLM runs
- Validates beta coefficients, standard errors, p-values

**Complete Coverage:**
- All 108 tests covered every ~12 weeks
- Both Linux and macOS platforms
- Multiple scenarios: missing/complete data, covariates, sample sizes, threading

**Validation Metrics:**
- Beta coefficient correlation ≥ 0.8
- Standard error correlation ≥ 0.8
- P-value correlation ≥ 0.8
- Configurable threshold in `PLINK2_GLM_TEST_CONFIG.sh`

## Configuration

### Batch System

Tests are organized in `CI_reference_files/` directory:

**Batch Files:**
- `BATCH_01_PARAMS.txt` through `BATCH_12_PARAMS.txt` (9 tests each)
- `ALL_BATCHES_MASTER_LIST.txt` (all 108 tests combined)

**Batch Distribution:**
- Tests divided evenly across 12 batches for parallel CI execution
- Each batch takes approximately similar runtime
- Batch selection via `BATCH_NUM` environment variable or config file

**Generating Batches:**
- Use `glm_test_batches.py` to regenerate batch files from master list
- Useful when adding new test configurations
- Ensures even distribution across batches

### Parameter Format

Test parameters stored in CSV format: `dataset_name,model,n,covariate,threads`

| Field | Options | Description |
|-------|---------|-------------|
| `dataset_name` | `1kgp3_50k_{nomiss\|yesmiss}_Av_nonintdose` | Dataset with/without missing genotypes |
| `model` | `linear`, `logistic`, `firth` | GLM model type (see below) |
| `n` | `1000`, `32000`, `32017` | Sample subset size |
| `covariate` | `COV_1` or blank | Covariate name (optional, blank for no covariates) |
| `threads` | `1`, `4`, `8` | PLINK2 thread count for parallelization |

**Model Types:**
- `linear` - Linear regression (uses phenotype `y`, continuous)
- `logistic` - Logistic regression (uses phenotype `ybool`, binary, no Firth correction)
- `firth` - Firth logistic regression (uses phenotype `ybool`, hybrid Firth correction)

**Example Configurations:**
```
1kgp3_50k_yesmiss_Av_nonintdose,linear,1000,COV_1,1
1kgp3_50k_nomiss_Av_nonintdose,logistic,32000,,4
1kgp3_50k_yesmiss_Av_nonintdose,firth,32017,COV_1,8
```

### Test Conditions

**Datasets:**
- `1kgp3_50k_nomiss_Av_nonintdose` - No missing genotypes (complete data)
- `1kgp3_50k_yesmiss_Av_nonintdose` - With missing genotypes (tests genotype QC)

**Sample Sizes:**
- **1000** - Small subset for quick validation
- **32000** - Large subset, multiple of 64 (64 × 500)
- **32017** - Large non-64-multiple subset (64 × 500 + 17), tests edge cases

**Threading Configurations:**
- **1 thread** - Single-threaded execution (tests sequential code path)
- **4 threads** - Multi-threaded (tests parallel code path)
- **8 threads** - Higher parallelization (tests scalability)

**Phenotypes:**
- **y** - Quantitative phenotype (for linear models)
- **ybool** - Binary phenotype (for logistic/Firth models)

**Covariates:**
- **COV_1** - Single covariate adjustment
- **None** - No covariate adjustment (blank field)

### Config File

Edit `PLINK2_GLM_TEST_CONFIG.sh` to customize test behavior:

**Path Configuration:**
```bash
datapath="test_data/"              # Test data directory
outdir="./results"                 # Results output directory
sbpath=""                          # Optional subdirectory path
compare_script="./2.0/build_dynamic/CI-SCRIPTS/COMPARE_GLM_PLINK2_R.py"
```

**Batch Selection:**
```bash
BATCH_NUM=${GLM_BATCH_NUM:-1}      # Batch number (1-12), defaults to 1
                                   # Override with GLM_BATCH_NUM env var (used by CI)
```

**Validation Threshold:**
```bash
correlation_threshold=0.8          # Minimum correlation for pass (0.0-1.0)
                                   # Tests fail if any correlation < threshold
```

**Advanced Configuration:**
- Modify `BATCH_FILE` path to use custom batch files
- Adjust correlation threshold for stricter/looser validation
- Change data paths for different test data locations

## Scripts Reference

| Script | Purpose | Key Features |
|--------|---------|--------------|
| `GLM_TESTS.sh` | CI entry point | Sets up Python venv, installs dependencies (pandas, numpy), creates results directory, executes test runner |
| `LOCAL_PLINK2_GLM_TESTING.sh` | Core test runner | Loads config, runs PLINK2 GLM for each test case, compares with R reference, validates correlations, fail-fast behavior |
| `PLINK2_GLM_TEST_CONFIG.sh` | Configuration | Defines paths, batch selection, correlation threshold, loads batch parameters from files |
| `COMPARE_GLM_PLINK2_R.py` | Validation script | Compares PLINK2 vs R outputs, calculates correlations (beta, p-values) and median differences (beta, SE, p-values), outputs pass/fail status |
| `glm_test_batches.py` | Batch generator | Generates batch files from master list, ensures even distribution across batches |
| `CLEANUP.sh` | Cleanup utility | Removes test_data, results directories, cleans build artifacts with `make clean` |
| `DYNAMICBUILD_LINUXMAC_ROOT.sh` | Build script | Compiles PLINK2 from source using system libraries (LAPACK, OpenBLAS) |

**Execution Flow:**
```
GLM_TESTS.sh → LOCAL_PLINK2_GLM_TESTING.sh → PLINK2_GLM_TEST_CONFIG.sh (config)
                                           → PLINK2 (run GLM)
                                           → COMPARE_GLM_PLINK2_R.py (validate)
```

## Directory Structure

```
CI-SCRIPTS/
├── README.md                                  # This documentation
│
├── Entry Point Scripts
│   ├── GLM_TESTS.sh                          # Main CI entry (sets up env)
│   └── LOCAL_PLINK2_GLM_TESTING.sh           # Core test runner (runs tests)
│
├── Configuration
│   └── PLINK2_GLM_TEST_CONFIG.sh             # Paths, batch, threshold config
│
├── Validation & Analysis
│   ├── COMPARE_GLM_PLINK2_R.py               # PLINK2 vs R comparison
│   └── glm_test_batches.py                   # Batch file generator
│
├── Utilities
│   ├── CLEANUP.sh                            # Cleanup test data/results
│   └── DYNAMICBUILD_LINUXMAC_ROOT.sh         # Build PLINK2
│
├── CI_reference_files/                       # Batch configurations
│   ├── BATCH_01_PARAMS.txt                   # Batch 1 (9 tests)
│   ├── BATCH_02_PARAMS.txt                   # Batch 2 (9 tests)
│   ├── ...                                   # Batches 3-11
│   ├── BATCH_12_PARAMS.txt                   # Batch 12 (9 tests)
│   └── ALL_BATCHES_MASTER_LIST.txt           # All 108 tests
│
└── results/                                  # Test outputs (gitignored)
    └── *.glm.{linear,logistic,logistic.hybrid}
```

## Troubleshooting

### Missing Test Data

**CI Environment:**
- **Symptom:** Workflow fails with "Missing: test_data/..." errors
- **Solution:**
  1. Run [cache_test_data.yml](/.github/workflows/cache_test_data.yml) workflow manually via workflow_dispatch
  2. Verify cache key matches (`test-data-v1`) in both workflows
  3. Check Google Drive file is accessible with file ID `1lW5s8kwMeOwHPrbRoREtO5pwgKKK_8FR`

**Local Environment:**
- **Symptom:** Cannot find test data files when running tests locally
- **Solution:**
  1. Contact repository maintainer for test data access
  2. Download from Google Drive using the file ID
  3. Extract to `test_data/` directory in repository root
  4. Verify all expected files are present (genotypes, phenotypes, reference outputs)

### Correlation Failures

**Symptoms:**
- Test fails with "❌ TEST FAILED - Correlation below threshold" message
- Output shows correlation values < 0.8

**Diagnosis:**
1. Check `results/` directory for detailed PLINK2 vs R comparison output
2. Review which metric failed (beta, SE, or p-value)
3. Examine PLINK2 .glm output file for unexpected values

**Common Causes:**
- **PLINK2 build issue:** Rebuild PLINK2 and rerun tests
- **Numerical precision:** Check for very small p-values or large betas
- **Data mismatch:** Verify R reference files match current test data version
- **Threading differences:** Compare single-threaded vs multi-threaded results

**Solutions:**
- Verify PLINK2 binary is built correctly from latest source
- Check that R reference files are from same PLINK2SYN run as test data
- Review test data integrity (md5sum check)
- Adjust correlation threshold if small numerical differences expected

### Python/Dependency Errors

**Missing Dependencies:**
```bash
# Install required Python packages
pip install pandas numpy

# Or with virtual environment
python3 -m venv venv
source venv/bin/activate
pip install pandas numpy
```

**Import Errors:**
- Verify Python 3.x is being used (`python3 --version`)
- Check virtual environment is activated if using one
- Ensure pip installed packages for correct Python version

**Script Execution Errors:**
- Make sure scripts have execute permissions: `chmod +x *.sh`
- Verify correct working directory (repository root)
- Check that all paths in config file are correct

### Workflow Not Triggering

**Automatic Trigger (Weekly):**
- **Symptom:** Workflow doesn't run on Sunday midnight
- **Check:**
  1. Verify new commits exist on master branch since last successful run
  2. Review workflow run history in GitHub Actions
  3. Check if previous run is still in progress
- **Solution:** Workflow skips if no new commits (expected behavior)

**Manual Trigger:**
- Navigate to Actions → "PLINK2 AUTO RUN TESTS - LINUX & MAC - CHECK COMMIT HX FIRST"
- Click "Run workflow" → Select branch → "Run workflow"

### Build Failures

**LAPACK/OpenBLAS Issues (Linux):**
```bash
# Install required libraries
sudo apt-get update
sudo apt-get install -y libopenblas-dev liblapack-dev liblapacke-dev
```

**macOS Build Issues:**
- Ensure Xcode command line tools installed: `xcode-select --install`
- Check for compatible system libraries
- Review build output in DYNAMICBUILD_LINUXMAC_ROOT.sh

**Compilation Errors:**
- Pull latest source code changes
- Clean build artifacts: `make -C 2.0/build_dynamic/ clean`
- Rebuild from scratch

