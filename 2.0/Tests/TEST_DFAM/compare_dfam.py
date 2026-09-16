#!/usr/bin/env python3
"""Compare a .dfam report against ref_dfam.py's output."""
import sys


def load_ref(path):
    out = {}
    with open(path) as f:
        f.readline()
        for line in f:
            g = line.split()
            out[g[0]] = (int(g[1]), float(g[2]), g[3])
    return out


def load_plink2(path):
    rows = []
    with open(path) as f:
        header = f.readline().split()
        idx = {name.lstrip('#'): i for i, name in enumerate(header)}
        for line in f:
            g = line.split()
            rows.append((g[idx['ID']], int(g[idx['OBS_CT']]),
                         float(g[idx['EXP_CT']]), g[idx['CHISQ']]))
    return rows


def load_plink19(path):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            g = line.split()
            rows.append((g[1], int(g[4]), float(g[5]), g[6]))
    return rows


def main():
    which, report_path, ref_path, tol = sys.argv[1:5]
    tol = float(tol)
    ref = load_ref(ref_path)
    rows = load_plink2(report_path) if which == 'plink2' else load_plink19(report_path)
    if len(rows) != len(ref):
        sys.exit('%d report rows, %d reference rows' % (len(rows), len(ref)))
    informative = 0
    for vid, obs, exp, chisq in rows:
        want_obs, want_exp, want_chisq = ref[vid]
        if obs != want_obs:
            sys.exit('%s: OBS_CT %d, expected %d' % (vid, obs, want_obs))
        if abs(exp - want_exp) > tol * max(1.0, abs(want_exp)):
            sys.exit('%s: EXP_CT %g, expected %g' % (vid, exp, want_exp))
        if (chisq == 'NA') != (want_chisq == 'NA'):
            sys.exit('%s: CHISQ %s, expected %s' % (vid, chisq, want_chisq))
        if chisq == 'NA':
            continue
        informative += 1
        if abs(float(chisq) - float(want_chisq)) > tol * max(1e-12, abs(float(want_chisq))):
            sys.exit('%s: CHISQ %s, expected %s' % (vid, chisq, want_chisq))
    if informative * 2 < len(rows):
        sys.exit('only %d of %d variants produced a statistic' %
                 (informative, len(rows)))
    print('%d variants verified (%d with a statistic)' % (len(rows), informative))


if __name__ == '__main__':
    main()
