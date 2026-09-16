# Independent weighted-distance oracle, written from the documented semantics:
#   dist_ij = (W / (W - wmiss_i - wmiss_j + wboth_ij)) * sum_v w_v * d_v
# over the variants both samples have called, where d_v is 0, 1 or 2 allele
# differences, w_v is the --distance-wts weight, and the missingness weights
# mw_v are w_v times the Hardy-Weinberg expected contribution.
#
# Standard library only: the test suite has no third-party dependencies.
import struct
import sys

prefix, wts_mode = sys.argv[1], sys.argv[2]
fam = [l.split() for l in open(prefix + '.fam') if l.strip()]
n = len(fam)
bim = [l.split() for l in open(prefix + '.bim') if l.strip()]
m = len(bim)
raw = open(prefix + '.bed', 'rb').read()[3:]
nb = (n + 3) // 4

# geno[j][i] is the ALT dosage, or None for a missing call
DOSAGE = (2.0, None, 1.0, 0.0)
geno = []
for j in range(m):
    row = [None] * n
    base = j * nb
    for i in range(n):
        code = (raw[base + (i >> 2)] >> (2 * (i & 3))) & 3
        row[i] = DOSAGE[code]
    geno.append(row)

# Allele frequencies over the non-missing calls, as plink computes them here.
freq = [0.0] * m
for j in range(m):
    called = [g for g in geno[j] if g is not None]
    freq[j] = (sum(called) / (2.0 * len(called))) if called else 0.0

if wts_mode.startswith('exp='):
    expo = float(wts_mode[4:])
    wts = []
    for q in freq:
        het = 2 * q * (1.0 - q)
        wts.append(max(het, 1e-300) ** (-expo) if het != 0.0 else 0.0)
else:
    wts = [-1.0] * m
    ids = {r[1]: j for j, r in enumerate(bim)}
    lines = open(wts_mode).read().split('\n')
    start = 0 if (len(sys.argv) > 3 and sys.argv[3] == 'noheader') else 1
    for line in lines[start:]:
        t = line.split()
        if len(t) >= 2 and t[0] in ids:
            wts[ids[t[0]]] = float(t[1])

keep = [j for j in range(m) if wts[j] > 0.0]
mw = {}
for j in keep:
    q = freq[j]
    mw[j] = q * (1.0 - q) * (q * q - q + 1.0) * wts[j]
total_mw = sum(mw.values())

wmiss = [0.0] * n
for j in keep:
    row = geno[j]
    for i in range(n):
        if row[i] is None:
            wmiss[i] += mw[j]

out = [[0.0] * n for _ in range(n)]
for i in range(n):
    for k in range(i + 1, n):
        rawd = 0.0
        wboth = 0.0
        for j in keep:
            gi = geno[j][i]
            gk = geno[j][k]
            if gi is None:
                if gk is None:
                    wboth += mw[j]
                continue
            if gk is None:
                continue
            rawd += wts[j] * abs(gi - gk)
        val = (total_mw / (total_mw - wmiss[i] - wmiss[k] + wboth)) * rawd
        out[i][k] = val
        out[k][i] = val

with open(sys.argv[-1], 'wb') as f:
    for row in out:
        f.write(struct.pack('<%dd' % n, *row))
