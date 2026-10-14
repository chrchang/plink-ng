#!/usr/bin/env python3
"""Recompute --neighbour's report from the .eigenvec it was derived from.

Usage: oracle.py <.eigenvec> <.nearest> <neighbor count>

Pure Python on purpose: the CI runners have no numpy.
"""
import math
import sys


def read_table(path):
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        assert header[0].startswith('#'), path
        header[0] = header[0][1:]
        rows = [line.rstrip('\n').split('\t') for line in f if line.strip()]
    return header, rows


def sample_key(header, row, id_cols):
    return tuple(row[header.index(col)] for col in id_cols)


def main():
    eigenvec_path, nearest_path, nn_ct = sys.argv[1], sys.argv[2], int(sys.argv[3])

    ev_header, ev_rows = read_table(eigenvec_path)
    pc_cols = [i for i, name in enumerate(ev_header) if name.startswith('PC')]
    assert pc_cols, 'no PC columns in ' + eigenvec_path
    coords = [[float(row[i]) for i in pc_cols] for row in ev_rows]
    sample_ct = len(coords)
    assert nn_ct < sample_ct

    # Same two passes as CalcNeighbour(): each sample's mean distance to its K
    # nearest neighbours, then the mean of those neighbours' own values.
    dist_self = []
    nn_idxs = []
    for i in range(sample_ct):
        dists = []
        for j in range(sample_ct):
            if j == i:
                continue
            d2 = sum((a - b) ** 2 for a, b in zip(coords[i], coords[j]))
            dists.append((d2, j))
        dists.sort()
        best = dists[:nn_ct]
        dist_self.append(sum(math.sqrt(d2) for d2, _ in best) / nn_ct)
        nn_idxs.append([j for _, j in best])

    expected = []
    for i in range(sample_ct):
        dist_nn = sum(dist_self[j] for j in nn_idxs[i]) / nn_ct
        stat = None if dist_nn == 0.0 else math.sqrt(dist_self[i] / dist_nn)
        expected.append((dist_self[i], dist_nn, stat))

    nb_header, nb_rows = read_table(nearest_path)
    assert len(nb_rows) == sample_ct, 'row count: %d vs %d' % (len(nb_rows), sample_ct)

    # cols= can drop FID or SID from the report, so compare on whichever ID
    # columns both files still carry.
    id_cols = [col for col in ('FID', 'IID', 'SID') if col in nb_header and col in ev_header]
    assert 'IID' in id_cols

    max_rel = 0.0
    for i, (row, (e_self, e_nn, e_stat)) in enumerate(zip(nb_rows, expected)):
        assert sample_key(nb_header, row, id_cols) == sample_key(ev_header, ev_rows[i], id_cols), \
            'sample order differs at row %d' % (i + 1)
        for col, want in (('DIST_SELF', e_self), ('DIST_NN', e_nn), ('STAT', e_stat)):
            if col not in nb_header:
                continue
            got_str = row[nb_header.index(col)]
            if want is None:
                assert got_str == 'NA', '%s row %d: expected NA, got %s' % (col, i + 1, got_str)
                continue
            got = float(got_str)
            rel = abs(got - want) / max(abs(want), 1e-300)
            max_rel = max(max_rel, rel)
            # The PC scores are re-read from the text .eigenvec, so this
            # recomputation starts from coordinates that have already been
            # rounded for printing; differences of that order are expected.
            assert rel < 1e-4, '%s row %d: %g vs %g (rel %g)' % (col, i + 1, got, want, rel)

    print('oracle.py: %d samples, K=%d, max relative difference %g' % (sample_ct, nn_ct, max_rel))


if __name__ == '__main__':
    main()
