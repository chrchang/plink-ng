#!/usr/bin/env python3
"""Compare plink2's .ldscore against the brute-force values."""
import sys


def read(path, id_col, val_start):
    out = {}
    order = []
    with open(path) as f:
        f.readline()
        for line in f:
            fields = line.split()
            out[fields[id_col]] = fields[val_start:]
            order.append(fields[id_col])
    return out, order


def main():
    got, got_order = read(sys.argv[1], 2, 3)
    want, want_order = read(sys.argv[2], 0, 1)
    if got_order != want_order:
        sys.exit('the two files list different variants, or in a different '
                 'order')
    checked = 0
    defined = 0
    worst = 0.0
    for vid, vals in want.items():
        for a, b in zip(got[vid], vals):
            if (a == 'NA') != (b == 'NA'):
                sys.exit('%s: got %s, expected %s' % (vid, a, b))
            if a == 'NA':
                checked += 1
                continue
            fa, fb = float(a), float(b)
            scale = max(abs(fa), abs(fb), 1e-300)
            rel = abs(fa - fb) / scale
            worst = max(worst, rel)
            # plink2 writes about six significant digits.
            if rel > 1e-5:
                sys.exit('%s: got %s, expected %s (relative %.2e)'
                         % (vid, a, b, rel))
            checked += 1
            defined += 1
    if not defined:
        sys.exit('every value was NA; the comparison proves nothing')
    print('%d LD Scores match the brute force (%d defined), worst relative '
          'difference %.2e' % (checked, defined, worst))


if __name__ == '__main__':
    main()
