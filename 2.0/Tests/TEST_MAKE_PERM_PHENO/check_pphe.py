#!/usr/bin/env python3
"""Check a .pphe against the phenotype it was permuted from.

A permutation report cannot be compared value-by-value against PLINK 1.9, since
the two draw from different random streams.  What must hold is that every
column is a rearrangement of the same values, that missing stays missing, and
that the columns actually differ from each other.
"""
import sys
from collections import Counter


def main():
    pheno_path, pphe_path = sys.argv[1], sys.argv[2]
    base = {}
    with open(pheno_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            g = line.split()
            base[g[1]] = g[2]
    rows = [line.rstrip('\n').split('\t') for line in open(pphe_path)]
    hdr, body = rows[0], rows[1:]
    perm_ct = sum(1 for h in hdr if h.startswith('PERM'))
    assert perm_ct, 'no PERM columns'
    id_col = hdr.index('IID')
    first = id_col + 1

    def norm(v):
        try:
            return '%.6g' % float(v)
        except ValueError:
            return v

    orig = Counter(norm(v) for v in base.values() if v != 'NA')
    assert orig, 'no nonmissing phenotype values'
    for k in range(perm_ct):
        col = Counter(norm(r[first + k]) for r in body if base[r[id_col]] != 'NA')
        assert col == orig, 'column %d is not a rearrangement of the input' % (k + 1)
    for r in body:
        if base[r[id_col]] == 'NA':
            for k in range(perm_ct):
                assert r[first + k] == 'NA', 'sample %s gained a phenotype' % r[id_col]
    cols = [tuple(r[first + k] for r in body) for k in range(perm_ct)]
    assert len(set(cols)) == perm_ct, 'only %d of %d columns are distinct' % (len(set(cols)), perm_ct)
    print('%d rows x %d permutations verified' % (len(body), perm_ct))


if __name__ == '__main__':
    main()
