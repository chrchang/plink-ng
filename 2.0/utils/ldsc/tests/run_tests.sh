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

# 9. Partitioned (stratified) LD Scores: one coefficient per annotation, and
#    the per-category output that follows from them.
$L --h2 part_trait.sumstats --ref-ld part_ref --w-ld part_w --out t_part
python3 check.py t_part.h2,t_part.results oracle.py h2 \
   --ref-ld part_ref.l2.ldscore --w-ld part_w.l2.ldscore \
   --sumstats part_trait.sumstats --M 3000,2000,1000
grep -q "Categories: 3" t_part.log
# The partitioned regression drops the two-step estimator for a chi^2 ceiling.
if grep -q "Using two-step estimator" t_part.log; then
    echo "the partitioned regression must not use the two-step estimator"
    exit 1
fi
grep -q "Removed .* variants with chi^2 >" t_part.log
python3 - << 'PYEOF'
import csv
row = next(csv.DictReader(open('t_part.h2'), delimiter='\t'))
est, se = float(row['h2']), float(row['h2_se'])
truth = 0.18 + 0.09 + 0.04
if abs(est - truth) > 4 * se:
    raise SystemExit('partitioned h2: %g is more than 4 SE (%g) from %g'
                     % (est, se, truth))
rows = list(csv.DictReader(open('t_part.results'), delimiter='\t'))
if len(rows) != 3:
    raise SystemExit('expected 3 category rows, got %d' % len(rows))
prop_sum = sum(float(r['Prop._h2']) for r in rows)
if abs(prop_sum - 1.0) > 1e-9:
    raise SystemExit('category h2 proportions sum to %g, not 1' % prop_sum)
print('partitioned h2 is within 4 SE, and its category proportions sum to 1')
PYEOF

# 10. --overlap-annot: annotations that overlap, corrected with the .annot
#     and .frq files.  The base annotation covers every variant, so its
#     proportion is 1 with no complement to test against.
$L --h2 ov_trait.sumstats --ref-ld ov_ref --w-ld ov_w --overlap-annot \
   --frqfile ov_ref --out t_ov
# The same run without the correction, to check the correction changes it.
$L --h2 ov_trait.sumstats --ref-ld ov_ref --w-ld ov_w --out t_ov_raw
M_OV=$(tr ' ' ',' < ov_ref.l2.M_5_50 | tr -d '\n')
python3 check.py t_ov.h2,t_ov.results oracle.py h2 \
   --ref-ld ov_ref.l2.ldscore --w-ld ov_w.l2.ldscore \
   --sumstats ov_trait.sumstats --M "$M_OV" --annot ov_ref.annot \
   --frqfile ov_ref.frq
python3 - << 'PYEOF'
import csv
rows = list(csv.DictReader(open('t_ov.results'), delimiter='\t'))
base = rows[0]
if abs(float(base['Prop._h2']) - 1.0) > 1e-9:
    raise SystemExit('the all-variant annotation must hold all the h2, got %s'
                     % base['Prop._h2'])
if base['Enrichment_p'] != 'NA':
    raise SystemExit('the all-variant annotation has no complement to test '
                     'against, so its enrichment p must be NA, got %s'
                     % base['Enrichment_p'])
# The correction has to actually change the answer.
raw = list(csv.DictReader(open('t_ov_raw.results'), delimiter='\t'))
changed = [abs(float(rows[i]['Prop._h2']) - float(raw[i]['Prop._h2']))
           for i in range(len(rows))]
if max(changed) < 1e-6:
    raise SystemExit('--overlap-annot left every proportion unchanged')
print('overlap correction: base holds all h2, its p is NA, and the '
      'corrected proportions differ from the uncorrected ones by up to %.3g'
      % max(changed))
PYEOF

# 11. --w-ld has to name exactly one LD Score column.
if $L --h2 part_trait.sumstats --ref-ld part_ref --w-ld part_ref \
      --out t_bad3 2> tmp_err3.txt; then
    echo "expected ldsc to reject multi-column --w-ld"
    exit 1
fi
grep -q "must name a single LD Score column" tmp_err3.txt

# 12. A missing Z column is an error, not a silent wrong answer.
cut -f1,2,3,5 trait1.sumstats > no_z.sumstats
if $L --h2 no_z.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
      --M 6000 --out t_bad 2> tmp_err.txt; then
    echo "expected ldsc to reject a file with no Z column"
    exit 1
fi
grep -q "must have SNP, Z and N columns" tmp_err.txt

# 13. --rg needs at least two filesets.
if $L --rg trait1.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
      --M 6000 --out t_bad2 2> tmp_err2.txt; then
    echo "expected ldsc to reject a single-file --rg"
    exit 1
fi
grep -q "at least two summary statistic files" tmp_err2.txt

echo "All ldsc tests passed."
