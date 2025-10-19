set -exo pipefail

## Data conditions - general for all
# test options (aside from core test data):
# --thin-indiv-count $((64*1000+17))
# --thin-indiv-count $((64*1000))
# --thin-indiv-count $((1000))
# --threads 1
# --threads 4


bfiles_=('1kgp3_50k_nomiss_Av_nonintdose_biallelic') # '1kgp3_50k_yesmiss_Av_nonintdose' - need to fix yesmiss
general_conditions_1=('--thin-indiv-count '$((1000))) #'--thin-indiv-count '$((64*1000+17)) '--thin-indiv-count '$((64*1000)) )
#general_conditions_1=('--thin-indiv-count 10')
general_conditions_2=('--threads 4') # '--threads 4')

pheno_=('y' 'ybool')

for bfile in "${bfiles_[@]}"; do
    for pheno in "${pheno_[@]}"; do
        for d1 in "${general_conditions_1[@]}"; do
            for d2 in "${general_conditions_2[@]}"; do
                if [[ "$pheno" == "ybool" ]]; then
                    # Use --logistic for binary phenotypes
                    plink --logistic --bfile "$bfile" --allow-extra-chr \
                    --pheno "/test_data/1kgp3_50k_nomiss_Av_phenocov.csv" --pheno-name "$pheno" \
                    --covar "/test_data/1kgp3_50k_nomiss_Av_phenocov.csv" --covar-name 'COV_1' \
                    $d1 $d2 --out ref
                else
                    # Use --linear for continuous phenotypes
                    plink --linear --bfile "$bfile" --allow-extra-chr \
                    --pheno "/test_data/1kgp3_50k_nomiss_Av_phenocov.csv" --pheno-name "$pheno" \
                    --covar "/test_data/1kgp3_50k_nomiss_Av_phenocov.csv" --covar-name 'COV_1' \
                    $d1 $d2 --out ref
                fi
                    plink2 --glm --bfile "$bfile" --allow-extra-chr \
                    --pheno "/test_data/1kgp3_50k_nomiss_Av_phenocov.csv" --pheno-name "$pheno" \
                    --covar "/test_data/1kgp3_50k_nomiss_Av_phenocov.csv" --covar-name 'COV_1' \
                    $d1 $d2 --out test
            done
        done
    done
done
