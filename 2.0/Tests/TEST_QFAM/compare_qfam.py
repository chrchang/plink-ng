#!/usr/bin/env python3
"""Compare a PLINK 2.0 .qfam report against PLINK 1.9's.

PLINK 1.9 counts its own A1 (the minor allele) and PLINK 2.0 counts ALT, so
BETA and STAT flip sign wherever the two disagree; that is aligned here rather
than papered over with an absolute value, so a genuine sign error would still
show up.  EMP1 comes from a different random stream in each program, so it is
only checked for agreement within Monte Carlo error.
"""
import math
import sys


def main():
    p19_path, p19_perm_path, p2_path, tol, emp1_tol = sys.argv[1:6]
    tol = float(tol)
    emp1_tol = float(emp1_tol)
    p19 = {}
    with open(p19_path) as f:
        f.readline()
        for line in f:
            g = line.split()
            p19[g[1]] = (g[3], g[5], g[6], g[7], g[8])
    p19_emp1 = {}
    with open(p19_perm_path) as f:
        f.readline()
        for line in f:
            g = line.split()
            p19_emp1[g[1]] = g[2]
    rows = []
    with open(p2_path) as f:
        header = f.readline().split()
        idx = {name.lstrip('#'): i for i, name in enumerate(header)}
        for line in f:
            g = line.split()
            rows.append(tuple(g[idx[name]] for name in
                              ('ID', 'A1', 'NIND', 'BETA', 'STAT', 'RAW_P', 'EMP1')))
    if len(rows) != len(p19):
        sys.exit('%d PLINK 2.0 rows, %d PLINK 1.9 rows' % (len(rows), len(p19)))
    flipped = 0
    checked = 0
    emp1_diffs = []
    for vid, a1, nind, beta, stat, rawp, emp1 in rows:
        a1_19, nind_19, beta_19, stat_19, rawp_19 = p19[vid]
        if nind != nind_19:
            sys.exit('%s: NIND %s, PLINK 1.9 says %s' % (vid, nind, nind_19))
        if (beta == 'NA') != (beta_19 == 'NA'):
            sys.exit('%s: BETA %s, PLINK 1.9 says %s' % (vid, beta, beta_19))
        if beta == 'NA':
            continue
        sign = 1.0 if a1 == a1_19 else -1.0
        if sign < 0:
            flipped += 1
        for name, got, want in (('BETA', sign * float(beta), float(beta_19)),
                                ('STAT', sign * float(stat), float(stat_19)),
                                ('RAW_P', float(rawp), float(rawp_19))):
            if abs(got - want) > tol * max(abs(want), 1e-9):
                sys.exit('%s: %s %g, PLINK 1.9 says %g' % (vid, name, got, want))
        checked += 1
        if emp1 != 'NA' and p19_emp1[vid] != 'NA':
            emp1_diffs.append(float(emp1) - float(p19_emp1[vid]))
    if checked * 2 < len(rows):
        sys.exit('only %d of %d variants were comparable' % (checked, len(rows)))
    if not flipped:
        sys.exit('no variant had opposite A1 assignments; the sign alignment is '
                 'untested')
    mean = sum(emp1_diffs) / len(emp1_diffs)
    sd = math.sqrt(sum((d - mean) ** 2 for d in emp1_diffs) / (len(emp1_diffs) - 1))
    if abs(mean) > emp1_tol:
        sys.exit('EMP1 differs systematically: mean difference %g' % mean)
    if sd > 4 * emp1_tol:
        sys.exit('EMP1 scatter %g is larger than Monte Carlo error explains' % sd)
    print('%d variants verified (%d with opposite A1), EMP1 mean diff %.4f, '
          'sd %.4f' % (checked, flipped, mean, sd))


if __name__ == '__main__':
    main()
