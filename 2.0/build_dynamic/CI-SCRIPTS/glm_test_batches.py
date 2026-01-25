#!/usr/bin/env python3
"""
Generate all permutations of GWAS test configurations,
shuffle them, and create separate batch parameter files.

Output: Creates batch parameter files in CI_reference_files/ subdirectory
"""

import random
import os
from itertools import product

# Configuration parameters for test permutations
datasets = [
    "1kgp3_50k_nomiss_Av_nonintdose",
    "1kgp3_50k_yesmiss_Av_nonintdose"
]

models = ["linear", "logistic", "firth"]
sample_sizes = ["1000", "32000", "32017"]
covariates = ["", "COV_1"]  # empty string means no covariate
threads = ["1", "4", "8"]

# Hardcoded output directory
BATCH_DIR = "CI_reference_files"

def generate_all_permutations():
    """Generate all 108 test configurations."""
    all_tests = []
    for dataset, model, n, cov, thread in product(datasets, models, sample_sizes, covariates, threads):
        test_config = f"{dataset},{model},{n},{cov},{thread}"
        all_tests.append(test_config)
    return all_tests

def write_batch_params(tests, batch_size=9):
    """Create separate parameter files for each batch."""
    
    # Create output directory if it doesn't exist
    os.makedirs(BATCH_DIR, exist_ok=True)
    
    # Shuffle tests
    shuffled_tests = tests.copy()
    random.shuffle(shuffled_tests)
    
    num_batches = (len(shuffled_tests) + batch_size - 1) // batch_size
    
    print(f"\nCreating {num_batches} batch parameter files in '{BATCH_DIR}/'...")
    
    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(shuffled_tests))
        batch = shuffled_tests[start_idx:end_idx]
        
        # Create filename with zero-padded batch number
        filename = f"BATCH_{batch_num + 1:02d}_PARAMS.txt"
        filepath = os.path.join(BATCH_DIR, filename)
        
        with open(filepath, 'w') as f:
            # Write header comment
            f.write(f"# Batch {batch_num + 1} of {num_batches}\n")
            f.write(f"# Tests {start_idx + 1}-{end_idx} of {len(shuffled_tests)}\n")
            f.write(f"# Contains {len(batch)} test configurations\n")
            f.write("#\n")
            f.write("# Format: dataset_name,model,n,covariate,threads\n")
            f.write("#\n\n")
            
            # Write parameters (one per line)
            for test in batch:
                f.write(f"{test}\n")
        
        print(f"  ✓ Created {filename} ({len(batch)} tests)")
    
    return shuffled_tests, num_batches

def write_master_list(shuffled_tests, batch_size=9):
    """Write a master list showing all tests in all batches."""
    filename = "ALL_BATCHES_MASTER_LIST.txt"
    filepath = os.path.join(BATCH_DIR, filename)
    
    with open(filepath, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("PLINK2 vs R GWAS GLM Testing Pipeline\n")
        f.write("MASTER BATCH LIST - All 108 Permutations\n")
        f.write("=" * 60 + "\n\n")
        
        num_batches = (len(shuffled_tests) + batch_size - 1) // batch_size
        
        for batch_num in range(num_batches):
            start_idx = batch_num * batch_size
            end_idx = min(start_idx + batch_size, len(shuffled_tests))
            batch = shuffled_tests[start_idx:end_idx]
            
            f.write(f"BATCH {batch_num + 1} (File: BATCH_{batch_num + 1:02d}_PARAMS.txt)\n")
            f.write("-" * 60 + "\n")
            for i, test in enumerate(batch, start=start_idx + 1):
                f.write(f"{i:3d}. {test}\n")
            f.write("\n")
        
        f.write("=" * 60 + "\n")
        f.write(f"Total: {len(shuffled_tests)} tests across {num_batches} batches\n")
        f.write("=" * 60 + "\n")
    
    print(f"  ✓ Created {filename}")

def main():
    # Set random seed for reproducibility
    random.seed(42)
    
    # Generate all test configurations
    all_tests = generate_all_permutations()
    
    print(f"Generated {len(all_tests)} test configurations")
    
    # Create batch parameter files
    shuffled_tests, num_batches = write_batch_params(all_tests, batch_size=9)
    
    # Write master list
    write_master_list(shuffled_tests, batch_size=9)
    
    print(f"\n{'='*60}")
    print(f"✓ Complete!")
    print(f"{'='*60}")
    print(f"Total tests: {len(all_tests)}")
    print(f"Batches created: {num_batches} (9 tests each)")
    print(f"\nFiles created in '{BATCH_DIR}/' directory:")
    print(f"  - BATCH_01_PARAMS.txt through BATCH_{num_batches:02d}_PARAMS.txt")
    print(f"  - ALL_BATCHES_MASTER_LIST.txt")
    print(f"\nThese parameter files can be loaded by your testing config script.")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()