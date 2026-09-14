#!/bin/bash

# --ld-score-annot, against a brute-force computation from the exported
# genotypes.

set -exo pipefail

$1/plink2 $2 $3 --dummy 150 300 0.05 --seed 7 --make-pgen --out tmp_data
python3 make_fixture.py 1
cp tmp_data.pgen tmp_pos.pgen
cp tmp_data.psam tmp_pos.psam

RADIUS_KB=20
RADIUS_BP=20000

# 1. Partitioned LD Scores, one column per annotation, against the brute
#    force.
$1/plink2 $2 $3 --pfile tmp_pos --ld-score --ld-score-window-kb $RADIUS_KB \
   --ld-score-annot annot.txt --out plink2_part
head -1 plink2_part.ldscore | grep -q $'codingL2\tconservedL2'
$1/plink2 $2 $3 --pfile tmp_pos --export A --out plink2_geno
python3 brute_ldscore.py plink2_geno.raw tmp_pos.pvar annot.txt $RADIUS_BP \
   > brute.txt
python3 compare.py plink2_part.ldscore brute.txt

# 2. The variant counts behind them.
python3 check_m.py plink2_geno.raw annot.txt plink2_part.ldscore.M \
   plink2_part.ldscore.M_5_50

# 3. An all-ones annotation is the unpartitioned score, exactly.
$1/plink2 $2 $3 --pfile tmp_pos --ld-score --ld-score-window-kb $RADIUS_KB \
   --out plink2_plain
$1/plink2 $2 $3 --pfile tmp_pos --ld-score --ld-score-window-kb $RADIUS_KB \
   --ld-score-annot annot_ones.txt --out plink2_ones
diff -q <(tail -n +2 plink2_plain.ldscore) <(tail -n +2 plink2_ones.ldscore)
# And its .M is just the variant count.
test "$(cat plink2_plain.ldscore.M)" = "$(cat plink2_ones.ldscore.M)"
test "$(cat plink2_plain.ldscore.M)" = "$(grep -vc '^#' tmp_pos.pvar)"

# 4. The other window conventions carry over, and a narrower window is a
#    subset: its window counts can only be smaller.  (Its scores can be
#    either way round, since the unbiased estimator makes some contributions
#    negative.)
$1/plink2 $2 $3 --pfile tmp_pos --ld-score cols=chrom,pos,nobsi,l2 \
   --ld-score-window-kb 20 --ld-score-annot annot.txt --out plink2_wide
$1/plink2 $2 $3 --pfile tmp_pos --ld-score cols=chrom,pos,nobsi,l2 \
   --ld-score-window-kb 5 --ld-score-annot annot.txt --out plink2_narrow
test "$(head -1 plink2_narrow.ldscore | tr '\t' '\n' | grep -c 'L2$')" -eq 2
python3 - << 'PYEOF'
def read(path):
    out = {}
    with open(path) as f:
        f.readline()
        for line in f:
            fields = line.split()
            out[fields[2]] = int(fields[3])
    return out

wide = read('plink2_wide.ldscore')
narrow = read('plink2_narrow.ldscore')
strictly_smaller = 0
for vid, count in narrow.items():
    if count > wide[vid]:
        raise SystemExit('%s: the 5kb window counted %d variants, more than '
                         'the 20kb window\'s %d' % (vid, count, wide[vid]))
    if count < wide[vid]:
        strictly_smaller += 1
if not strictly_smaller:
    raise SystemExit('the two windows covered the same variants; the '
                     'comparison proves nothing')
print('%d windows shrank with the narrower radius, none grew'
      % strictly_smaller)
PYEOF

# 4b. The variant-count window convention also works with annotations.
$1/plink2 $2 $3 --pfile tmp_pos --ld-score --ld-score-window 5 \
   --ld-score-annot annot.txt --out plink2_win
test "$(head -1 plink2_win.ldscore | tr '\t' '\n' | grep -c 'L2$')" -eq 2

# 5. A variant with no annotation row is an error, not a partial answer.
if $1/plink2 $2 $3 --pfile tmp_pos --ld-score --ld-score-window-kb $RADIUS_KB \
      --ld-score-annot annot_short.txt --out plink2_bad 2> tmp_err.txt; then
    echo "expected --ld-score-annot to reject an incomplete annotation file"
    exit 1
fi
grep -q "is missing from annot_short.txt" tmp_err.txt

# 6. An annotation file with no variant ID column is an error.
awk 'NR == 1 { print "coding\tconserved"; next } { print $2 "\t" $3 }' \
   annot.txt > annot_noid.txt
if $1/plink2 $2 $3 --pfile tmp_pos --ld-score --ld-score-window-kb $RADIUS_KB \
      --ld-score-annot annot_noid.txt --out plink2_bad2 2> tmp_err2.txt; then
    echo "expected --ld-score-annot to reject a file with no SNP column"
    exit 1
fi
grep -q "no SNP (or ID) column" tmp_err2.txt

echo "TEST_LD_SCORE_ANNOT passed."
