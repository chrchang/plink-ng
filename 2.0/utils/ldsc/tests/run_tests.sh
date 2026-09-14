#!/bin/bash

# ldsc, against an independent implementation of the same estimator
# (oracle.py) and against the model the fixtures were simulated under.
#
# Usage: ./run_tests.sh <directory containing the ldsc binary>

set -exo pipefail

L=$1/ldsc

python3 make_data.py 1

# 1. Heritability, default settings: free intercept, two-step estimator.
$L --h2 trait1.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
   --M 6000 --out t_h2
python3 check.py t_h2.h2 oracle.py h2 --ref-ld ldscores.ldscore \
   --w-ld w_ld.ldscore --sumstats trait1.sumstats --M 6000

# 2. Intercept constrained to 1, which also turns the two-step estimator off.
$L --h2 trait1.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
   --M 6000 --no-intercept --out t_h2_noint
python3 check.py t_h2_noint.h2 oracle.py h2 --ref-ld ldscores.ldscore \
   --w-ld w_ld.ldscore --sumstats trait1.sumstats --M 6000 --intercept-h2 1
grep -q "Intercept: constrained to 1" t_h2_noint.log

# 3. A different two-step cutoff, and a different block count.
$L --h2 trait1.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
   --M 6000 --two-step 10 --n-blocks 50 --out t_h2_ts
python3 check.py t_h2_ts.h2 oracle.py h2 --ref-ld ldscores.ldscore \
   --w-ld w_ld.ldscore --sumstats trait1.sumstats --M 6000 --two-step 10 \
   --n-blocks 50

# 4. Genetic correlation.
$L --rg trait1.sumstats,trait2.sumstats --ref-ld ldscores.ldscore \
   --w-ld w_ld.ldscore --M 6000 --out t_rg
python3 check.py t_rg.rg oracle.py rg --ref-ld ldscores.ldscore \
   --w-ld w_ld.ldscore --sumstats trait1.sumstats,trait2.sumstats --M 6000

# 5. The estimates have to land near the parameters the data was simulated
#    under; matching the oracle is not enough on its own.
python3 - << 'PYEOF'
import csv
row = next(csv.DictReader(open('t_rg.rg'), delimiter='\t'))
checks = [('h2_p1', 0.30, float(row['h2_p1_se'])),
          ('h2_obs', 0.18, float(row['h2_obs_se'])),
          ('rg', 0.10 / (0.30 * 0.18) ** 0.5, float(row['se']))]
for name, truth, se in checks:
    est = float(row[name])
    if abs(est - truth) > 4 * se:
        raise SystemExit('%s: %g is more than 4 SE (%g) from the simulated %g'
                         % (name, est, se, truth))
print('h2 and rg are within 4 SE of the simulated values')
PYEOF

# 6. Allele orientation: swapping A1/A2 (with Z negated) and flipping strand
#    must not change anything.
$L --rg trait1.sumstats,trait2_flipped.sumstats --ref-ld ldscores.ldscore \
   --w-ld w_ld.ldscore --M 6000 --out t_rg_flip
diff -q <(tail -n +2 t_rg.rg | cut -f3-) <(tail -n +2 t_rg_flip.rg | cut -f3-)

# 7. --ref-ld-chr/--w-ld-chr, with M read from the .l2.M_5_50 files, has to
#    reproduce the single-file run.
$L --h2 trait1.sumstats --ref-ld-chr chr_ref. --w-ld-chr chr_w. --out t_h2_chr
grep -q "Read M = 6000 from the .l2.M_5_50 files" t_h2_chr.log
diff -q t_h2.h2 t_h2_chr.h2

# 8. --M defaults to the .l2.M_5_50 file next to --ref-ld.
$L --h2 trait1.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
   --out t_h2_mfile
diff -q t_h2.h2 t_h2_mfile.h2

# 9. A missing Z column is an error, not a silent wrong answer.
cut -f1,2,3,5 trait1.sumstats > no_z.sumstats
if $L --h2 no_z.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
      --M 6000 --out t_bad 2> tmp_err.txt; then
    echo "expected ldsc to reject a file with no Z column"
    exit 1
fi
grep -q "must have SNP, Z and N columns" tmp_err.txt

# 10. --rg needs at least two filesets.
if $L --rg trait1.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
      --M 6000 --out t_bad2 2> tmp_err2.txt; then
    echo "expected ldsc to reject a single-file --rg"
    exit 1
fi
grep -q "at least two summary statistic files" tmp_err2.txt

echo "All ldsc tests passed."
