# Independent weighted-distance oracle, written from the documented semantics:
#   dist_ij = (W / (W - wmiss_i - wmiss_j + wboth_ij)) * sum_v w_v * d_v
# over the variants both samples have called, where d_v is 0, 1 or 2 allele
# differences, w_v is the --distance-wts weight, and the missingness weights
# mw_v are w_v times the Hardy-Weinberg expected contribution.
import sys
import numpy as np

prefix, wts_mode = sys.argv[1], sys.argv[2]
fam = [l.split() for l in open(prefix + '.fam')]
n = len(fam)
bim = [l.split() for l in open(prefix + '.bim')]
m = len(bim)
raw = open(prefix + '.bed', 'rb').read()[3:]
nb = (n + 3) // 4
lut = np.array([2.0, np.nan, 1.0, 0.0])
G = np.empty((m, n))
for j in range(m):
    b = np.frombuffer(raw[j * nb:(j + 1) * nb], dtype=np.uint8)
    codes = np.unpackbits(b, bitorder='little').reshape(-1, 2)
    G[j] = lut[(codes[:, 0] + 2 * codes[:, 1])[:n]]

called = ~np.isnan(G)
# Allele frequencies over the non-missing calls, as plink computes them here.
ac = np.nansum(np.where(called, G, 0.0), axis=1)
an = 2.0 * called.sum(axis=1)
with np.errstate(invalid='ignore', divide='ignore'):
    q = np.where(an > 0, ac / np.maximum(an, 1), 0.0)

if wts_mode.startswith('exp='):
    expo = float(wts_mode[4:])
    het = 2 * q * (1.0 - q)
    w = np.where(het != 0.0, np.power(np.maximum(het, 1e-300), -expo), 0.0)
else:
    w = np.full(m, -1.0)
    ids = {r[1]: j for j, r in enumerate(bim)}
    lines = open(wts_mode).read().split('\n')
    if len(sys.argv) > 3 and sys.argv[3] == 'noheader':
        start = 0
    else:
        start = 1
    for line in lines[start:]:
        t = line.split()
        if len(t) < 2:
            continue
        if t[0] in ids:
            w[ids[t[0]]] = float(t[1])

keep = w > 0.0
w = np.where(keep, w, 0.0)
mw = q * (1.0 - q) * (q * q - q + 1.0) * w
mw = np.where(keep, mw, 0.0)
W = mw.sum()

Gk = np.where(called, G, 0.0)[keep]
ck = called[keep]
wk = w[keep]
mwk = mw[keep]
wmiss = ((~ck) * mwk[:, None]).sum(axis=0)

out = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        both = ck[:, i] & ck[:, j]
        rawd = float((wk[both] * np.abs(Gk[both, i] - Gk[both, j])).sum())
        wboth = float(mwk[(~ck[:, i]) & (~ck[:, j])].sum())
        val = (W / (W - wmiss[i] - wmiss[j] + wboth)) * rawd
        out[i, j] = val
        out[j, i] = val
out.astype(np.float64).tofile(sys.argv[-1])
