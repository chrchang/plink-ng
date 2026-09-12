#!/bin/bash

# --grm-maf, and the frequency floor --pca/GRM construction now insists on.

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

# 400 samples, so the instability point is 0.25 / sqrt(400) = 0.0125.  Half the
# variants sit below it and half well above.
{ printf '500 rare 0.0005 0.01 1 1\n'; printf '500 common 0.1 0.5 1 1\n'; } > tmp_mix.sim
plink --simulate tmp_mix.sim --simulate-ncases 200 --simulate-ncontrols 200 --out tmp_mix > /dev/null

printf '500 common 0.1 0.5 1 1\n' > tmp_common.sim
plink --simulate tmp_common.sim --simulate-ncases 200 --simulate-ncontrols 200 --out tmp_common > /dev/null

# 1. A monomorphic variant gets its own message, and 'yes-really' does not
#    bypass it: there is no minor allele to standardize by.  Simulating in the
#    0.0005-0.01 range at 400 samples reliably produces some.
fails $1/plink2 $2 $3 --bfile tmp_mix --pca 4 --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --out plink2_bad > tmp_err.txt 2>&1 || true
grep -q 'carry no minor allele' tmp_err.txt
fails $1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0 yes-really --out plink2_bad

# 1b. With those excluded, what remains is the instability error, and its
#     message names the lowest frequency present.
fails $1/plink2 $2 $3 --bfile tmp_mix --mac 1 --pca 4 --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_mix --mac 1 --pca 4 --out plink2_bad > tmp_err2.txt 2>&1 || true
grep -q 'lowest remaining here is' tmp_err2.txt

# 2. A common-only fileset needs no flag.
$1/plink2 $2 $3 --bfile tmp_common --pca 4 --out plink2_common
test -s plink2_common.eigenvec

# 3. --grm-maf makes it run, and reports what it dropped.
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.05 --out plink2_maf
test -s plink2_maf.eigenvec
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.05 --out plink2_maf > tmp_out.txt 2>&1
grep -q 'variants remaining for --pca/GRM construction' tmp_out.txt

# 4. --grm-min-af is the same flag.
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-min-af 0.05 --out plink2_alias
diff -q plink2_maf.eigenvec plink2_alias.eigenvec
diff -q plink2_maf.eigenval plink2_alias.eigenval

# 5. It filters for --pca/GRM alone: --freq still sees every variant.
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.05 --freq --out plink2_both
test "$(grep -vc '^#' plink2_both.afreq)" -eq 1000

#    and the eigenvectors match the run without --freq.
diff -q plink2_maf.eigenvec plink2_both.eigenvec

# 6. Restricting the global variant set the same way gives the same PCs, which
#    is what "applied just to --pca" has to mean.
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --maf 0.05 --out plink2_globalmaf
diff -q plink2_maf.eigenvec plink2_globalmaf.eigenvec
diff -q plink2_maf.eigenval plink2_globalmaf.eigenval

# 7. A threshold below the instability point needs 'yes-really'.
fails $1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.001 --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.001 yes-really --out plink2_yr
test -s plink2_yr.eigenvec

# 8. The mode parameter --maf takes is accepted, and nonsense is not.
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.05 minor --out plink2_mode
$1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.05 nonmajor yes-really --out plink2_mode2
diff -q plink2_mode.eigenvec plink2_mode2.eigenvec
fails $1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 0.05 nonsense --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf abc --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf 1.5 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_mix --pca 4 --grm-maf -1 --out plink2_bad

# 9. GRM construction is gated the same way, whether or not --pca asked for it.
fails $1/plink2 $2 $3 --bfile tmp_mix --make-rel --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_mix --make-grm-list --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_mix --make-rel --grm-maf 0.05 --out plink2_rel
test -s plink2_rel.rel

# 10. A threshold that removes everything is an error rather than a crash.
fails $1/plink2 $2 $3 --bfile tmp_common --pca 4 --grm-maf 0.99 --out plink2_bad
