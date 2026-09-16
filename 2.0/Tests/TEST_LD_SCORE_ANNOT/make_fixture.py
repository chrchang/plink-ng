#!/usr/bin/env python3
"""Give the dummy fileset real positions, and write an annotation file.

Deliberately numpy-free: the test runners do not have numpy.
"""
import random
import sys

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
rand = random.Random(seed)
lines = open('tmp_data.pvar').read().rstrip('\n').split('\n')
out = []
ids = []
pos = 0
for line in lines:
    if line.startswith('#'):
        out.append(line)
        continue
    fields = line.split('\t')
    # Spacings on the order of the window radius, so the windows vary in size.
    pos += rand.randint(500, 4000)
    fields[0] = '1'
    fields[1] = str(pos)
    ids.append(fields[2])
    out.append('\t'.join(fields))
open('tmp_pos.pvar', 'w').write('\n'.join(out) + '\n')

# One binary annotation and one continuous one, plus an all-ones column in a
# separate file: partitioning by that has to reproduce the unpartitioned
# score exactly.
with open('annot.txt', 'w') as f:
    f.write('SNP\tcoding\tconserved\n')
    for vid in ids:
        f.write('%s\t%d\t%.4f\n'
                % (vid, 1 if rand.random() < 0.3 else 0, rand.uniform(0, 1)))
with open('annot_ones.txt', 'w') as f:
    f.write('SNP\tall\n')
    for vid in ids:
        f.write('%s\t1\n' % vid)
# An annotation file missing one variant, which has to be an error rather
# than a silently partial answer.
with open('annot_short.txt', 'w') as f:
    f.write('SNP\tcoding\n')
    for vid in ids[:-1]:
        f.write('%s\t1\n' % vid)
print('%d variants' % len(ids))
