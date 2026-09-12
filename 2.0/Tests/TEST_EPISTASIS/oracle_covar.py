#!/usr/bin/env python3
"""Independent logistic oracle for --epistasis-boost's covariate adjustment.

Reads a PLINK 1 fileset and a covariate file, and for every pair of variants
fits

  logit P(case) ~ 1 + covariates + genotype-of-A + genotype-of-B      (main)
  logit P(case) ~ 1 + covariates + one term per occupied (A, B) cell  (full)

printing the likelihood ratio statistic and its degrees of freedom.

Standard library only: the fits are plain Newton-Raphson, and the normal
equations of each step are solved by Gaussian elimination with partial
pivoting.  The design is at most 9 + covariate count columns wide, so this
costs nothing next to the pass over the samples.

The coding here is deliberately not the one plink2 uses: the likelihood ratio
does not depend on which cell is taken as the reference, so a different choice
is a check on plink2's rather than a copy of it.
"""
import math
import sys


def solve(a, b):
    """Solves a x = b in place, or returns None if a is singular."""
    k = len(b)
    for col in range(k):
        pivot = max(range(col, k), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-11:
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


def loglik(design, y, beta):
    total = 0.0
    for row, yi in zip(design, y):
        eta = sum(bj * xj for bj, xj in zip(beta, row))
        if eta > 0.0:
            total += (yi - 1.0) * eta - math.log1p(math.exp(-eta))
        else:
            total += yi * eta - math.log1p(math.exp(eta))
    return total


def logistic_fit(design, y):
    """Newton-Raphson.  Returns the log-likelihood, or None if it fails."""
    k = len(design[0])
    beta = [0.0] * k
    prev = None
    for _ in range(50):
        xtwx = [[0.0] * k for _ in range(k)]
        grad = [0.0] * k
        for row, yi in zip(design, y):
            eta = sum(bj * xj for bj, xj in zip(beta, row))
            if eta > 30.0:
                p = 1.0 - 1e-13
            elif eta < -30.0:
                p = 1e-13
            else:
                p = 1.0 / (1.0 + math.exp(-eta))
            w = p * (1.0 - p)
            resid = yi - p
            for a in range(k):
                xa = row[a]
                if xa == 0.0:
                    continue
                grad[a] += xa * resid
                for b in range(a, k):
                    xtwx[a][b] += xa * w * row[b]
        for a in range(k):
            for b in range(a):
                xtwx[a][b] = xtwx[b][a]
        step = solve([r[:] for r in xtwx], grad[:])
        if step is None:
            return None
        beta = [bj + sj for bj, sj in zip(beta, step)]
        cur = loglik(design, y, beta)
        if prev is not None and abs(cur - prev) < 1e-11 * (abs(cur) + 0.1):
            return cur
        prev = cur
    return None


def main():
    prefix = sys.argv[1]
    covar_fname = sys.argv[2]
    fam = [l.split() for l in open(prefix + '.fam') if l.strip()]
    sample_ids = [(f[0], f[1]) for f in fam]
    n = len(fam)
    # plink1 .fam case/control: 1 = control, 2 = case, anything else missing.
    pheno = [1.0 if f[5] == '2' else (0.0 if f[5] == '1' else None) for f in fam]
    bim = [l.split() for l in open(prefix + '.bim') if l.strip()]
    ids = [b[1] for b in bim]
    m = len(bim)

    raw = open(prefix + '.bed', 'rb').read()
    assert raw[:3] == b'\x6c\x1b\x01'
    stride = (n + 3) // 4
    # plink1 code 0 = hom A1, 1 = missing, 2 = het, 3 = hom A2
    CODE = (2, None, 1, 0)
    geno = []
    for v in range(m):
        base = 3 + v * stride
        geno.append([CODE[(raw[base + (i >> 2)] >> (2 * (i & 3))) & 3] for i in range(n)])

    rows = {}
    for line in open(covar_fname).read().split('\n')[1:]:
        f = line.split()
        if f:
            rows[(f[0], f[1])] = [float('nan') if x in ('NA', 'nan', '-9') else float(x)
                                  for x in f[2:]]
    covars = [rows.get(sid, [float('nan')]) for sid in sample_ids]
    keep = [i for i in range(n)
            if pheno[i] is not None and all(x == x for x in covars[i])]
    covar_ct = len(covars[keep[0]]) if keep else 0
    # A covariate that is constant over the kept samples is collinear with the
    # intercept, and plink2 drops it.
    live = [c for c in range(covar_ct)
            if any(covars[i][c] != covars[keep[0]][c] for i in keep)]

    print('ID1\tID2\tSTAT\tDF')
    out = []
    for vi in range(m):
        for vj in range(vi + 1, m):
            gi = geno[vi]
            gj = geno[vj]
            members = [i for i in keep if gi[i] is not None and gj[i] is not None]
            cells = {}
            for i in members:
                cells.setdefault((gi[i], gj[i]), []).append(i)
            row_levels = sorted({g for g, _ in cells})
            col_levels = sorted({g for _, g in cells})
            if len(row_levels) < 2 or len(col_levels) < 2:
                continue
            cell_keys = sorted(cells)
            main_ct = 1 + len(live) + (len(row_levels) - 1) + (len(col_levels) - 1)
            full_ct = 1 + len(live) + len(cell_keys) - 1
            df = full_ct - main_ct
            if df < 1 or len(members) <= full_ct:
                continue
            # The reference level is the last one here, and the reference cell
            # the last occupied one, where plink2 takes the first of each.
            y = []
            main = []
            full = []
            for i in members:
                base = [1.0] + [covars[i][c] for c in live]
                a_dummies = [1.0 if gi[i] == g else 0.0 for g in row_levels[:-1]]
                b_dummies = [1.0 if gj[i] == g else 0.0 for g in col_levels[:-1]]
                cell_dummies = [1.0 if (gi[i], gj[i]) == k else 0.0
                                for k in cell_keys[:-1]]
                main.append(base + a_dummies + b_dummies)
                full.append(base + cell_dummies)
                y.append(pheno[i])
            ll_main = logistic_fit(main, y)
            ll_full = logistic_fit(full, y)
            if ll_main is None or ll_full is None:
                continue
            stat = 2.0 * (ll_full - ll_main)
            out.append('%s\t%s\t%.10g\t%u' % (ids[vi], ids[vj], max(stat, 0.0), df))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
