#!/usr/bin/env python3
"""Independent reimplementation of plink2's phenotype/covariate transforms.

Usage: check_transform.py <original file> <transformed file> <spec...>
Each spec is <column>=<op>[+<op>...], with ops applied left to right:
  vs: variance-standardize (mean 0, sample variance 1)
  qn: quantile-normalize (Phi^-1((rank - 0.5) / n), ties get the mean rank)
  id: unchanged
Every column not named in a spec must be unchanged.  NA stays NA.
"""

import statistics
import sys


def load(fname):
    with open(fname) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f]
    return header, rows


def vs(vals):
    present = [v for v in vals if v is not None]
    mean = statistics.fmean(present)
    sd = statistics.stdev(present)
    return [None if v is None else (v - mean) / sd for v in vals]


def qn(vals):
    present = sorted(v for v in vals if v is not None)
    n = len(present)
    nd = statistics.NormalDist()
    rank = {}
    i = 0
    while i < n:
        j = i
        while j < n and present[j] == present[i]:
            j += 1
        rank[present[i]] = (i + 1 + j) / 2.0
        i = j
    return [None if v is None else nd.inv_cdf((rank[v] - 0.5) / n) for v in vals]


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    header, orig = load(sys.argv[1])
    new_header, new = load(sys.argv[2])
    if header != new_header or len(orig) != len(new):
        sys.exit("check_transform.py: shape mismatch between %s and %s" % (sys.argv[1], sys.argv[2]))
    specs = dict(spec.split("=") for spec in sys.argv[3:])
    for col_idx, name in enumerate(header):
        old_col = [row[col_idx] for row in orig]
        new_col = [row[col_idx] for row in new]
        ops = specs.get(name, "id").split("+")
        if ops == ["id"]:
            if old_col != new_col:
                sys.exit("check_transform.py: column %s changed: %s -> %s" % (name, old_col, new_col))
            continue
        vals = [None if x == "NA" else float(x) for x in old_col]
        for op in ops:
            vals = {"vs": vs, "qn": qn}[op](vals)
        for i, (expected, got) in enumerate(zip(vals, new_col)):
            if expected is None:
                if got != "NA":
                    sys.exit("check_transform.py: %s row %d is %s, expected NA" % (name, i, got))
            elif abs(float(got) - expected) > 1e-5:
                sys.exit("check_transform.py: %s row %d is %s, expected %.6g" % (name, i, got, expected))


main()
