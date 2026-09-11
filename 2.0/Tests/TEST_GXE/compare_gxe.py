#!/usr/bin/env python3
"""Compare plink2's .gxe against PLINK 1.9's .qassoc.gxe.

PLINK 1.9 regresses on the minor-allele count and plink2 on the ALT count, so
BETA can differ in sign; the magnitudes, and Z_GXE and P, cannot.  1.9 writes 4
significant digits, which sets the tolerance.
"""
import sys

TOL = 2e-3


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def main():
    p19 = {}
    with open(sys.argv[1]) as f:
        for line in f:
            g = line.split()
            if g[0] == 'CHR':
                continue
            p19[g[1]] = tuple(float(x) for x in (g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9]))
    hdr = None
    n = 0
    flipped = 0
    with open(sys.argv[2]) as f:
        for line in f:
            g = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                hdr = g
                continue
            d = dict(zip(hdr, g))
            vid = d['ID']
            assert vid in p19, 'variant %s missing from the PLINK 1.9 report' % vid
            n19_1, b1, s1, n19_2, b2, s2, z, p = p19[vid]
            n += 1
            assert int(d['OBS_CT1']) == int(n19_1), vid
            assert int(d['OBS_CT2']) == int(n19_2), vid
            g1, g2 = float(d['BETA1']), float(d['BETA2'])
            assert rel(abs(g1), abs(b1)) < TOL, '%s BETA1 %g vs %g' % (vid, g1, b1)
            assert rel(abs(g2), abs(b2)) < TOL, '%s BETA2 %g vs %g' % (vid, g2, b2)
            assert rel(float(d['SE1']), s1) < TOL, vid
            assert rel(float(d['SE2']), s2) < TOL, vid
            assert rel(abs(float(d['Z_GXE'])), abs(z)) < TOL, '%s Z %s vs %g' % (vid, d['Z_GXE'], z)
            assert rel(float(d['P']), p) < TOL, '%s P %s vs %g' % (vid, d['P'], p)
            if (g1 * b1) < 0:
                flipped += 1
    assert n >= 50, 'only %d variants compared' % n
    print('%d variants match PLINK 1.9 (%d with the expected BETA sign flip)' % (n, flipped))


if __name__ == '__main__':
    main()
