#!/usr/bin/env python3
"""Partitioned LD Scores computed the slow, obvious way, for comparison.

Reads the genotypes as exported by --export A, so it shares no code with
plink2's r^2 path, and applies the same unbiased estimator:
r^2 - (1 - r^2) / (n - 2), with n the pair's non-missing count.

Deliberately numpy-free: the test runners do not have numpy.
"""
import sys


def read_geno(path):
    with open(path) as f:
        header = f.readline().split()
        ids = [h.rsplit('_', 1)[0] for h in header[6:]]
        rows = [[None if x == 'NA' else float(x) for x in line.split()[6:]]
                for line in f if line.strip()]
    cols = [[row[j] for row in rows] for j in range(len(ids))]
    return ids, cols


def unbiased_r2(a, b):
    xs, ys = [], []
    for i in range(len(a)):
        if a[i] is None or b[i] is None:
            continue
        xs.append(a[i])
        ys.append(b[i])
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return None
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    r2 = min(sxy * sxy / (sxx * syy), 1.0)
    return r2 - (1.0 - r2) / (n - 2)


def main():
    raw_path, pvar_path, annot_path, radius = sys.argv[1:5]
    radius = int(radius)
    ids, cols = read_geno(raw_path)
    pos = {}
    for line in open(pvar_path):
        if line.startswith('#'):
            continue
        fields = line.split()
        pos[fields[2]] = int(fields[1])
    with open(annot_path) as f:
        annot_names = f.readline().split()[1:]
        annots = {}
        for line in f:
            fields = line.split()
            annots[fields[0]] = [float(x) for x in fields[1:]]

    print('ID\t' + '\t'.join(n + 'L2' for n in annot_names))
    for j, vid in enumerate(ids):
        present = [v for v in cols[j] if v is not None]
        if len(set(present)) < 2:
            # Monomorphic: no correlation with anything, so undefined.
            print('%s\t%s' % (vid, '\t'.join(['NA'] * len(annot_names))))
            continue
        # The variant's own weight, since r^2 with itself is 1 and the
        # correction vanishes.
        l2 = list(annots[vid])
        for k, other in enumerate(ids):
            if k == j or abs(pos[other] - pos[vid]) > radius:
                continue
            val = unbiased_r2(cols[j], cols[k])
            if val is None:
                continue
            for c in range(len(annot_names)):
                weight = annots[other][c]
                if weight:
                    l2[c] += weight * val
        print('%s\t%s' % (vid, '\t'.join('%.10g' % v for v in l2)))


if __name__ == '__main__':
    main()
