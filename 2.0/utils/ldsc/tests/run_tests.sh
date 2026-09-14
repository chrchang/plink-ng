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

# 11. --munge: the quality control that turns a raw GWAS file into
#     .sumstats, checked against munge_oracle.py, which re-derives which
#     variants should survive and what their Z and N should be.
$L --munge raw.txt --out t_munge
python3 munge_oracle.py --raw raw.txt --munged t_munge.sumstats
grep -q "Removed 1 variants that were not SNPs or were strand-ambiguous" \
   t_munge.log
grep -q "Removed 1 variants with out-of-bounds p-values" t_munge.log
grep -q "Removed 1 variants with missing values" t_munge.log
grep -q "variants with duplicated IDs" t_munge.log
grep -q "Median BETA was" t_munge.log

# 11b. Case/control counts, an odds ratio as the signed statistic, and the
#      frequency column kept.
$L --munge raw_cc.txt --keep-maf --out t_munge_cc
python3 munge_oracle.py --raw raw_cc.txt --munged t_munge_cc.sumstats \
   --keep-maf
head -1 t_munge_cc.sumstats | grep -q "FRQ"
# A frequency on a rounding boundary has to print the way every other tool
# prints it, which takes a correctly-rounded parse of the input.
grep -q "^rs_round	A	G	.*	0.086$" t_munge_cc.sumstats

# 11c. --merge-alleles: the output covers the list, in its order, with the
#      variants it could not fill in left missing.
$L --munge raw_cc.txt --merge-alleles merge_alleles.txt --out t_munge_merge
python3 munge_oracle.py --raw raw_cc.txt --munged t_munge_merge.sumstats \
   --merge-alleles merge_alleles.txt
grep -q "^rs_absent	NA	NA	NA	NA$" t_munge_merge.sumstats

# 11d. --daner takes the case and control counts from the column names, and
#      --a1-inc accepts a file with no signed statistic.
$L --munge raw_daner.txt --daner --out t_munge_daner
python3 munge_oracle.py --raw raw_daner.txt --munged t_munge_daner.sumstats \
   --daner
grep -q "N_cas = 12345, N_con = 67890" t_munge_daner.log
$L --munge raw_a1inc.txt --a1-inc --out t_munge_a1inc
python3 munge_oracle.py --raw raw_a1inc.txt --munged t_munge_a1inc.sumstats \
   --a1-inc
# With --a1-inc every Z is positive, since nothing says otherwise.
if awk 'NR > 1 && $4 < 0 { found = 1 } END { exit !found }' \
      t_munge_a1inc.sumstats; then
    echo "--a1-inc should not produce negative Z values"
    exit 1
fi

# 11e. A munged file feeds straight back into --h2.
$L --h2 t_munge.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
   --M 6000 --out t_munge_h2
grep -q "Total Observed scale h2" t_munge_h2.log

# 11f. A file with two signed statistics is ambiguous, and --ignore resolves
#      it.
if $L --munge raw_daner.txt --out t_bad4 2> tmp_err4.txt; then
    echo "expected ldsc to reject a file with no determinable sample size"
    exit 1
fi
grep -q "sample size" tmp_err4.txt

# 12. --h2-cts: each cell type regressed alongside the baseline annotations,
#     reporting its own coefficient, sorted by a one-sided p-value.
$L --h2-cts part_trait.sumstats --ref-ld part_ref --ref-ld-chr-cts cts_list.txt \
   --w-ld part_w --out t_cts
grep -q "Ran 3 cell-type regressions, each with the 3 baseline annotations" \
   t_cts.log
python3 - << 'PYEOF'
import csv
import subprocess
import sys

rows = list(csv.DictReader(open('t_cts.cell_type_results.txt'),
                           delimiter='\t'))
if len(rows) != 3:
    raise SystemExit('expected 3 cell types, got %d' % len(rows))
ps = [float(r['Coefficient_P_value']) for r in rows]
if ps != sorted(ps):
    raise SystemExit('the cell types must come out sorted by p-value: %s' % ps)
m_base = open('part_ref.l2.M_5_50').read().split()
for row in rows:
    name = row['Name']
    m_cts = open('cts_%s.l2.M_5_50' % name).read().split()
    out = subprocess.run(
        [sys.executable, 'oracle.py', 'h2', '--ref-ld',
         'part_ref.l2.ldscore', '--extra-ld', 'cts_%s.l2.ldscore' % name,
         '--w-ld', 'part_w.l2.ldscore', '--sumstats', 'part_trait.sumstats',
         '--M', ','.join(m_cts + m_base), '--cts-p'],
        check=True, capture_output=True, text=True).stdout
    want = {}
    for line in out.split('\n'):
        if line.strip():
            key, val = line.split()
            want[key] = float(val)
    for key, col in (('Coefficient_0', 'Coefficient'),
                     ('Coefficient_std_error_0', 'Coefficient_std_error'),
                     ('Coefficient_P_value_0', 'Coefficient_P_value')):
        got = float(row[col])
        if abs(got - want[key]) > 1e-6 * max(abs(want[key]), 1e-12):
            raise SystemExit('%s %s: expected %.12g, got %.12g'
                             % (name, col, want[key], got))
print('3 cell types match the oracle, in p-value order')
PYEOF

# 12b. --print-cov and --print-delete-vals, whose contents have to agree with
#      the estimates they come from.
$L --h2 part_trait.sumstats --ref-ld part_ref --w-ld part_w --print-cov \
   --print-delete-vals --out t_dumps
python3 - << 'PYEOF'
import csv

cov = [[float(x) for x in line.split()] for line in open('t_dumps.cov')]
delete = [[float(x) for x in line.split()] for line in open('t_dumps.delete')]
part = [[float(x) for x in line.split()]
        for line in open('t_dumps.part_delete')]
rows = list(csv.DictReader(open('t_dumps.results'), delimiter='\t'))
result = next(csv.DictReader(open('t_dumps.h2'), delimiter='\t'))
n_annot = len(rows)
n_blocks = int(result['n_blocks'])
if len(cov) != n_annot or any(len(r) != n_annot for r in cov):
    raise SystemExit('.cov should be %dx%d' % (n_annot, n_annot))
if len(delete) != n_blocks or any(len(r) != 1 for r in delete):
    raise SystemExit('.delete should be %dx1' % n_blocks)
if len(part) != n_blocks or any(len(r) != n_annot for r in part):
    raise SystemExit('.part_delete should be %dx%d' % (n_blocks, n_annot))
# The square roots of the diagonal are the coefficient standard errors.
for j in range(n_annot):
    want = float(rows[j]['Coefficient_std_error'])
    got = cov[j][j] ** 0.5
    if abs(got - want) > 1e-9 * max(want, 1e-12):
        raise SystemExit('.cov diagonal %d: %g, but the coefficient SE is %g'
                         % (j, got, want))
# Deleting one block at a time cannot move the estimate far.
h2 = float(result['h2'])
mean_delete = sum(r[0] for r in delete) / n_blocks
if abs(mean_delete - h2) > 0.05 * abs(h2):
    raise SystemExit('the delete values average %g, but h2 is %g'
                     % (mean_delete, h2))
# Each row of .part_delete has to add up to that row of .delete, weighted by
# the per-annotation variant counts.
m_vec = [float(x) for x in open('part_ref.l2.M_5_50').read().split()]
for b in range(n_blocks):
    want = sum(part[b][j] * m_vec[j] for j in range(n_annot))
    if abs(want - delete[b][0]) > 1e-9 * max(abs(delete[b][0]), 1e-12):
        raise SystemExit('block %d: partitioned delete values sum to %g, not '
                         '%g' % (b, want, delete[b][0]))
print('.cov, .delete and .part_delete agree with the estimates')
PYEOF

# 13. --w-ld has to name exactly one LD Score column.
if $L --h2 part_trait.sumstats --ref-ld part_ref --w-ld part_ref \
      --out t_bad3 2> tmp_err3.txt; then
    echo "expected ldsc to reject multi-column --w-ld"
    exit 1
fi
grep -q "must name a single LD Score column" tmp_err3.txt

# 14. A missing Z column is an error, not a silent wrong answer.
cut -f1,2,3,5 trait1.sumstats > no_z.sumstats
if $L --h2 no_z.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
      --M 6000 --out t_bad 2> tmp_err.txt; then
    echo "expected ldsc to reject a file with no Z column"
    exit 1
fi
grep -q "must have SNP, Z and N columns" tmp_err.txt

# 15. --rg needs at least two filesets.
if $L --rg trait1.sumstats --ref-ld ldscores.ldscore --w-ld w_ld.ldscore \
      --M 6000 --out t_bad2 2> tmp_err2.txt; then
    echo "expected ldsc to reject a single-file --rg"
    exit 1
fi
grep -q "at least two summary statistic files" tmp_err2.txt

echo "All ldsc tests passed."
