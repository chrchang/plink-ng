import pandas as pd
import numpy as np
import sys
import re

def COMPARE_GLM_PLINK2_R(plink_file, R_file):
    """
    Compare PLINK 2.0 and R glm regression results (linear or logistic).
    
    For logistic models: compares log(OR) vs R beta.
    For linear models: compares beta vs R beta.

    Args:
        plink_file (str): Path to PLINK 2.0 output (--glm) file (tab-delimited).
        R_file (str): Path to R glm output file (comma-delimited).

    Returns:
        tuple: (correlation_beta, correlation_p, differences_ok)
    """
    # --- Read input files ---
    df_plink = pd.read_csv(plink_file, sep="\t", dtype=str)
    df_R = pd.read_csv(R_file, sep=",", dtype=str)

    # --- Convert numeric columns ---
    for col in ["OR", "BETA", "SE", "LOG(OR)_SE", "Z_STAT", "P"]:
        if col in df_plink.columns:
            df_plink[col] = pd.to_numeric(df_plink[col], errors="coerce")

    for col in ["beta", "se", "p"]:
        if col in df_R.columns:
            df_R[col] = pd.to_numeric(df_R[col], errors="coerce")

    # --- Create SNP key for merge ---
    df_plink["SNP_KEY"] = df_plink["ID"].astype(str)
    df_R["SNP_KEY"] = df_R["SNP"].str.replace(r"_[ACGT]+$", "", regex=True)

    # --- Merge on SNP ID ---
    merged = pd.merge(df_plink, df_R, on="SNP_KEY", how="inner", suffixes=("_plink", "_R"))
    n_merge = merged.shape[0]
    print(f"✅ Merged {n_merge} SNPs between PLINK and R output.")

    # --- Detect model type ---
    if "OR" in merged.columns and merged["OR"].notna().any():
        model_type = "logistic"
    elif "BETA" in merged.columns and merged["BETA"].notna().any():
        model_type = "linear"
    else:
        raise ValueError("Could not detect model type from PLINK file (missing OR/BETA columns).")

    print(f"📊 Detected model type: {model_type.upper()}")

    # --- Compute comparable betas ---
    if model_type == "logistic":
        merged["beta_plink"] = -np.log(merged["OR"])
    else:
        merged["beta_plink"] = -merged["BETA"]

    # --- Compute correlations ---
    corr_beta = merged[["beta_plink", "beta"]].corr().iloc[0,1]
    corr_p = merged[["P", "p"]].corr().iloc[0,1]

    # --- Correlation ignores scale, so also compare the values themselves ---
    # Medians, since a few variants may legitimately differ near separation.
    se_col = "SE" if "SE" in merged.columns else "LOG(OR)_SE"
    diff_beta = (merged["beta_plink"] - merged["beta"]).abs().median() / merged["beta"].abs().median()
    diff_se = (merged[se_col] - merged["se"]).abs().median() / merged["se"].abs().median()
    diff_p = (merged["P"] - merged["p"]).abs().median()
    # shortcut: 5% was not calibrated against CI data, tighten it once the weekly job has reported values
    differences_ok = bool(diff_beta <= 0.05 and diff_se <= 0.05 and diff_p <= 0.05)

    # --- Print quick summary ---
    print(f"\n🔹 Correlation between PLINK and R betas/log(OR): {corr_beta:.4f}")
    print(f"🔹 Correlation between p-values: {corr_p:.4f}")

    # No decimal point: the shell caller reads any decimal in this output as a correlation.
    print(f"\n🔹 Median differences: beta {diff_beta:.0%}, SE {diff_se:.0%}, P {diff_p:.0%} (limit 5%)")
    if not differences_ok:
        print("❌ PLINK and R values differ by more than the limit.")

    # --- Return ---
    return corr_beta, corr_p, differences_ok


# --- Command-line execution ---
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <plink_file> <R_file>")
        sys.exit(1)

    plink_file = sys.argv[1]
    R_file = sys.argv[2]

    corr_beta, corr_p, differences_ok = COMPARE_GLM_PLINK2_R(plink_file, R_file)
    if not differences_ok:
        sys.exit(1)
