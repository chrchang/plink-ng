#!/usr/bin/env python3
"""Independent least-squares oracle for --epistasis.

Reads a PLINK 1 fileset and an optional covariate file, fits

  phenotype ~ 1 + A + B + A*B + covariates

over the samples with both genotypes called, and prints the interaction
term's coefficient, standard error and t-statistic for every pair.

Standard library only: the test suite has no third-party dependencies, so the
normal equations are formed directly and solved by Gaussian elimination with
partial pivoting.  The matrix is (4 + covariate count) square, so this costs
nothing next to the pass over the samples.
"""
import math
import sys


def solve(a, b):
    """Solves a x = b in place, or returns None if a is singular."""
    k = len(b)
    for col in range(k):
        pivot = max(range(col, k), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-12:
            return None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            b[col], b[pivot] = b[pivot], b[col]
        inv = 1.0 / a[col][col]
        for row in range(col + 1, k):
            factor = a[row][col] * inv
            if factor == 0.0:
                continue
            for j in range(col, k):
                a[row][j] -= factor * a[col][j]
            b[row] -= factor * b[col]
    x = [0.0] * k
    for col in range(k - 1, -1, -1):
        acc = b[col]
        for j in range(col + 1, k):
            acc -= a[col][j] * x[j]
        x[col] = acc / a[col][col]
    return x


def inverse_diag(a, idx):
    """Element idx of the inverse's diagonal, from a fresh copy of a."""
    k = len(a)
    unit = [1.0 if i == idx else 0.0 for i in range(k)]
    col = solve([row[:] for row in a], unit)
    return None if col is None else col[idx]


prefix = sys.argv[1]
covar_fname = sys.argv[2] if len(sys.argv) > 2 else None

fam = [l.split() for l in open(prefix + '.fam') if l.strip()]
sample_ids = [(f[0], f[1]) for f in fam]
pheno = [float(f[5]) for f in fam]
n = len(fam)
bim = [l.split() for l in open(prefix + '.bim') if l.strip()]
ids = [b[1] for b in bim]
m = len(bim)

raw = open(prefix + '.bed', 'rb').read()
assert raw[:3] == b'\x6c\x1b\x01'
stride = (n + 3) // 4
# plink1 code 0 = hom A1, 1 = missing, 2 = het, 3 = hom A2; plink2 counts A1
DOSAGE = (2.0, None, 1.0, 0.0)
geno = []
for v in range(m):
    base = 3 + v * stride
    geno.append([DOSAGE[(raw[base + (i >> 2)] >> (2 * (i & 3))) & 3] for i in range(n)])

keep = [i for i in range(n) if pheno[i] == pheno[i]]
covars = [[] for _ in range(n)]
if covar_fname:
    rows = {}
    for line in open(covar_fname).read().split('\n')[1:]:
        f = line.split()
        if f:
            rows[(f[0], f[1])] = [float('nan') if x in ('NA', 'nan', '-9') else float(x)
                                  for x in f[2:]]
    covars = [rows.get(sid, [float('nan')]) for sid in sample_ids]
    keep = [i for i in keep if all(x == x for x in covars[i])]

covar_ct = len(covars[keep[0]]) if (covar_fname and keep) else 0
param_ct = 4 + covar_ct

print('ID1\tID2\tBETA_INT\tSE\tT_STAT')
out = []
for vi in range(m):
    for vj in range(vi + 1, m):
        gi = geno[vi]
        gj = geno[vj]
        rows = [i for i in keep if gi[i] is not None and gj[i] is not None]
        if len(rows) <= param_ct:
            continue
        xtx = [[0.0] * param_ct for _ in range(param_ct)]
        xty = [0.0] * param_ct
        yty = 0.0
        for i in rows:
            a = gi[i]
            b = gj[i]
            x = [1.0, a, b, a * b] + covars[i]
            y = pheno[i]
            yty += y * y
            for p in range(param_ct):
                xp = x[p]
                if xp == 0.0:
                    continue
                xty[p] += xp * y
                row = xtx[p]
                for q in range(p, param_ct):
                    row[q] += xp * x[q]
        for p in range(param_ct):
            for q in range(p):
                xtx[p][q] = xtx[q][p]
        beta = solve([row[:] for row in xtx], xty[:])
        if beta is None:
            continue
        rss = yty - sum(beta[p] * xty[p] for p in range(param_ct))
        df = len(rows) - param_ct
        if df <= 0 or rss <= 0.0:
            continue
        inv_diag = inverse_diag(xtx, 3)
        if inv_diag is None or inv_diag <= 0.0:
            continue
        se = math.sqrt(inv_diag * rss / df)
        out.append('%s\t%s\t%.10g\t%.10g\t%.10g' % (ids[vi], ids[vj], beta[3], se, beta[3] / se))
print('\n'.join(out))
