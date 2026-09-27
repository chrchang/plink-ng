#!/usr/bin/env python3
"""Independent check of --flip-scan dprime's D' values.

Usage: dprime_oracle.py <.raw from --export A> <.flipscan.verbose>

For every line of the verbose report, recomputes case-only and control-only
signed D' from the ALT-allele counts in the .raw file and compares them with
D_PRIME_A and D_PRIME_U.  The haplotype frequencies are fitted by a direct
search over how the double heterozygotes split between the two phasings (grid,
then golden section), rather than by the cubic solve plink2 uses.  Standard
library only.
"""

import math
import sys


def em_dprime(g0, g1):
    # 3x3 genotype table over the samples nonmissing at both variants.
    tab = [[0] * 3 for _ in range(3)]
    for a, b in zip(g0, g1):
        if a is not None and b is not None:
            tab[a][b] += 1
    n = sum(sum(row) for row in tab)
    if n < 2:
        return None
    # Haplotypes known from everything but the double heterozygotes.
    # Index: 0 = (0,0), 1 = (0,1), 2 = (1,0), 3 = (1,1), counting ALT copies.
    known = [0.0] * 4
    for a in range(3):
        for b in range(3):
            c = tab[a][b]
            if not c or (a == 1 and b == 1):
                continue
            # With at most one heterozygous variant, both haplotypes are
            # determined: the homozygous variant's allele sits on both.
            alleles0 = {0: (0, 0), 1: (0, 1), 2: (1, 1)}[a]
            alleles1 = {0: (0, 0), 1: (0, 1), 2: (1, 1)}[b]
            for u, v in zip(alleles0, alleles1):
                known[2 * u + v] += c
    hh = tab[1][1]
    tot = 2.0 * n
    # ALT frequencies at the two variants.
    p0 = sum(a * tab[a][b] for a in range(3) for b in range(3)) / tot
    p1 = sum(b * tab[a][b] for a in range(3) for b in range(3)) / tot
    if min(p0, 1 - p0, p1, 1 - p1) < 1e-12:
        return None

    def freqs(share):
        # share = fraction of the double heterozygotes resolved as
        # (0,0)+(1,1) rather than (0,1)+(1,0).
        return [(known[0] + hh * share) / tot, (known[1] + hh * (1 - share)) / tot,
                (known[2] + hh * (1 - share)) / tot, (known[3] + hh * share) / tot]

    def lnlike(share):
        # Observed-data log-likelihood (up to a constant).  Every
        # maximum-likelihood haplotype frequency vector is of the freqs()
        # form, so maximizing over share in [0, 1] finds it.
        f = freqs(share)
        total = 0.0
        for k in range(4):
            if known[k]:
                if f[k] <= 0:
                    return -math.inf
                total += known[k] * math.log(f[k])
        if hh:
            cross = f[0] * f[3] + f[1] * f[2]
            if cross <= 0:
                return -math.inf
            total += hh * math.log(cross)
        return total

    # Grid (the likelihood can have two local maxima), then golden-section
    # refinement around every local maximum on it; the best one wins.
    grid_ct = 400
    shares = [i / grid_ct for i in range(grid_ct + 1)]
    vals = [lnlike(x) for x in shares]
    ratio = (math.sqrt(5) - 1) / 2
    best_share = None
    best_val = -math.inf
    for i in range(grid_ct + 1):
        if (i and vals[i - 1] > vals[i]) or (i < grid_ct and vals[i + 1] > vals[i]):
            continue
        lo = shares[max(i - 1, 0)]
        hi = shares[min(i + 1, grid_ct)]
        for _ in range(80):
            m1 = hi - ratio * (hi - lo)
            m2 = lo + ratio * (hi - lo)
            if lnlike(m1) < lnlike(m2):
                lo = m1
            else:
                hi = m2
        for cand in (shares[i], lo, hi, 0.5 * (lo + hi)):
            val = lnlike(cand)
            if val > best_val:
                best_share = cand
                best_val = val
    f = freqs(best_share)
    d = f[3] - p0 * p1
    if d >= 0:
        dmax = min(p0 * (1 - p1), (1 - p0) * p1)
    else:
        dmax = min(p0 * p1, (1 - p0) * (1 - p1))
    if dmax <= 0:
        return None
    return d / dmax


def main():
    raw_fname, verbose_fname = sys.argv[1], sys.argv[2]
    with open(raw_fname) as f:
        header = f.readline().split()
        rows = [line.split() for line in f]
    col_of = {}
    for i, name in enumerate(header[6:], start=6):
        col_of[name.rsplit("_", 1)[0]] = i
    cases = [r for r in rows if r[5] == "2"]
    ctrls = [r for r in rows if r[5] == "1"]

    def geno(group, vid):
        c = col_of[vid]
        return [None if r[c] == "NA" else int(r[c]) for r in group]

    with open(verbose_fname) as f:
        vheader = f.readline().lstrip("#").rstrip("\n").split("\t")
        ci = vheader.index("ID_INDEX")
        cp = vheader.index("ID_PAIR")
        ca = vheader.index("D_PRIME_A")
        cu = vheader.index("D_PRIME_U")
        checked = 0
        worst = 0.0
        for line in f:
            fields = line.rstrip("\n").split("\t")
            for group, col in ((cases, ca), (ctrls, cu)):
                expected = em_dprime(geno(group, fields[ci]), geno(group, fields[cp]))
                got = float(fields[col])
                if expected is None:
                    print("oracle has no D' for %s/%s but plink2 reported %s" % (fields[ci], fields[cp], fields[col]))
                    sys.exit(1)
                err = abs(expected - got)
                worst = max(worst, err)
                if err > 2e-5:
                    print("D' mismatch for %s/%s: plink2 %s, oracle %.8g" % (fields[ci], fields[cp], fields[col], expected))
                    sys.exit(1)
                checked += 1
    if checked == 0:
        print("no verbose lines to check")
        sys.exit(1)
    print("checked %d D' values, max abs difference %.3g" % (checked, worst))


main()
