#!/usr/bin/env python3
"""Recompute --het's fhat1/fhat2/fhat3 straight from a VCF and an .afreq, and
compare against the report.  The multiallelic terms have no PLINK 1.9
equivalent to check against, so this is the independent implementation."""
import sys

vcf_path, afreq_path, het_path = sys.argv[1:4]

freqs = {}
with open(afreq_path) as f:
    hdr = f.readline().split()
    id_idx = hdr.index('ID')
    alt_idx = hdr.index('ALT_FREQS')
    for line in f:
        r = line.split()
        alt = [float(x) for x in r[alt_idx].split(',')]
        freqs[r[id_idx]] = [1.0 - sum(alt)] + alt

sums = {}
nobs = {}
samples = None
for line in open(vcf_path):
    if line.startswith('##'):
        continue
    r = line.rstrip('\n').split('\t')
    if line.startswith('#CHROM'):
        samples = r[9:]
        for s in samples:
            sums[s] = [0.0, 0.0, 0.0]
            nobs[s] = 0
        continue
    a = freqs[r[2]]
    ssq = sum(x * x for x in a)
    ehet = 1.0 - ssq
    if ehet < 2.0**-35:   # --het skips monomorphic variants
        continue
    f1_den = 2 * ehet
    f1_base = 4 * ssq
    for s, gt in zip(samples, r[9:]):
        if gt.startswith('.'):
            continue
        lo, hi = sorted(int(x) for x in gt.split('/'))
        if lo == hi:
            f1 = (f1_base + 4 - 8 * a[lo]) / f1_den
            f2 = 2.0
            f3 = 1.0 + (1.0 - a[lo]) / a[lo]
        else:
            f1 = (f1_base + 2 - 4 * (a[lo] + a[hi])) / f1_den
            f2 = 2.0 - 1.0 / ehet
            f3 = 0.0
        cur = sums[s]
        cur[0] += f1
        cur[1] += f2
        cur[2] += f3
        nobs[s] += 1

checked = 0
with open(het_path) as f:
    hdr = f.readline().split()
    assert hdr == ['#IID', 'OBS_CT', 'FHAT1', 'FHAT2', 'FHAT3'], hdr
    for line in f:
        r = line.split()
        s = r[0]
        n = int(r[1])
        if n != nobs[s]:
            sys.exit('OBS_CT on %s: %d vs %d' % (s, n, nobs[s]))
        for k in range(3):
            want = sums[s][k] / nobs[s] - 1.0
            got = float(r[2 + k])
            # the report prints 6 significant digits
            if abs(got - want) > 1e-5 * max(abs(want), 1.0):
                sys.exit('FHAT%d on %s: %g vs %g' % (k + 1, s, got, want))
        checked += 1
if not checked:
    sys.exit('no samples checked')
print('%d samples matched on multiallelic data' % checked)
