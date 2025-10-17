import pandas as pd
import matplotlib.pyplot as plt

# === Step 1. Read input files ===

x_file = "test_data/1kgp3_50k_nomiss_Av_nonintdose_combined_known_g_effects.csv"         # file with variant IDs or x-values
data_file = "derivatives/1kgp3_50k_nomiss_Av_nonintdose_combined_phenocov.csv_ybool_glm.ybool.glm.logistic.hybrid"      # file with variant data

x_df = pd.read_csv(x_file)      # assumes header: variant_id
data_df = pd.read_csv(data_file, sep="\t")

# === Step 2. Merge on variant ID ===
merged = pd.merge(x_df, data_df, left_on='variant IDs', right_on='ID', how='inner')

# print(merged.head(100))

# === Step 3. Plot ===
plt.figure(figsize=(8,5))
plt.scatter(merged['x_values'].values, merged['OR'].values, label='BETA', s=4)

# # Optional: if you want -log10(P)
# import numpy as np
# if 'P' in merged.columns:
#     plt.scatter(merged['variant IDs'], -np.log10(merged['P']), label='-log10(P)', color='red')

# plt.xlabel('Variant ID')
# plt.ylabel('Effect Size / -log10(P)')
# plt.title('Variant Association Plot')
# plt.xticks(rotation=90)
# plt.legend()
# plt.tight_layout()
plt.show()