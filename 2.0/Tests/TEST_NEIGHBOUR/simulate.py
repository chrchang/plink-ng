#!/usr/bin/env python3
"""Write a .ped/.map pair with a main cluster and a handful of outliers.

Each outlier is homozygous across its own block of variants, so each one ends
up alone in a different direction of PC space.  They have to be separated from
each other as well as from the cluster: the statistic is a *local* outlier
factor, and a tight group of identical samples reads as a small cluster, not
as outliers.  Deterministic, so the expected ranking does not move between
runs.
"""
import random
import sys

main_ct = 56
outlier_ct = 4
variant_ct = 2000
block_size = variant_ct // outlier_ct

random.seed(20260915)

with open('tmp_nb.map', 'w') as f:
    for vidx in range(variant_ct):
        f.write('1\tv%d\t0\t%d\n' % (vidx + 1, (vidx + 1) * 1000))

def common_genotype():
    # allele frequency 0.5, so the main cluster has no structure of its own
    return ' '.join(random.choice(['A', 'C']) for _ in range(2))

with open('tmp_nb.ped', 'w') as f:
    for sidx in range(main_ct + outlier_ct):
        is_outlier = sidx >= main_ct
        outlier_idx = sidx - main_ct
        iid = 'out%d' % (outlier_idx + 1) if is_outlier else 'samp%d' % (sidx + 1)
        row = ['fam%d' % (sidx + 1), iid, '0', '0', '1', '-9']
        block_start = outlier_idx * block_size
        block_end = block_start + block_size
        for vidx in range(variant_ct):
            if is_outlier and block_start <= vidx < block_end:
                row.append('C C')
            else:
                row.append(common_genotype())
        f.write('\t'.join(row) + '\n')

with open('tmp_nb_outliers.txt', 'w') as f:
    for oidx in range(outlier_ct):
        f.write('out%d\n' % (oidx + 1))

sys.stderr.write('simulate.py: %d samples, %d variants\n' % (main_ct + outlier_ct, variant_ct))
