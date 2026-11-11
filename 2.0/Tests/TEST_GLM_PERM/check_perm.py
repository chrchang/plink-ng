"""Checks for --glm permutation test output.

  check_perm.py report <dump> <mperm report>
      EMP1_CT and EMP2_CT in a 'perm-count' max(T) report must match the
      counts recomputed from the --mperm-save-all dump.
  check_perm.py ties <dump> <mperm report>
      Like 'report', but every near tie must count as exactly 1/2.  Use this
      only for a phenotype where near ties are exact ties, e.g. one with two
      values, where a permutation that reproduces the original genotype x
      phenotype table gives the same statistic up to rounding.
  check_perm.py adaptive <dump> <aperm report> <allele test ct> <max perms>
      Replays adaptive permutation (default --aperm parameters other than the
      maximum) on the permutations in a max(T) dump from a run with the same
      --seed, and checks each allele's EMP1_CT and PERM_CT.
  check_perm.py identity <mperm report> <perm ct>
      For identity permutations (--permute-within with one sample per
      category), every permutation reproduces the original statistics, so
      each valid allele must have EMP1_CT = N/2, and EMP2_CT must be N/2 for
      the most extreme allele and N for the others.

Dumps contain |ln p| values (higher is more extreme) with 6 significant
digits, so values within a relative 1e-6 of the original statistic are "near
ties" that may count as more extreme, tied, or less extreme; the reported
count must fall in the resulting range.
"""
import sys
from statistics import NormalDist


def read_dump(fname):
    rows = [line.split() for line in open(fname)]
    to_f = lambda vals: [None if v == "NA" else float(v) for v in vals]
    return to_f(rows[0][1:]), [to_f(r[1:]) for r in rows[1:]]


def read_report(fname):
    lines = [line.rstrip("\n").split("\t") for line in open(fname)]
    hdr = lines[0]
    return [dict(zip(hdr, line)) for line in lines[1:]]


def near(v, o):
    return abs(v - o) <= 1e-6 * abs(o) + 1e-12


def count_range(vals, o, strict=False):
    """(low, high) x2 counts of values at least as extreme as o, ties = 1/2."""
    lo = hi = 0
    for v in vals:
        if v is None:
            continue
        if near(v, o):
            if strict:
                lo += 1
                hi += 1
            else:
                hi += 2
        elif v > o:
            lo += 2
            hi += 2
    return lo, hi


def check_report(dump_fname, report_fname, strict=False):
    orig, perms = read_dump(dump_fname)
    report = read_report(report_fname)
    bests = [max((v for v in p if v is not None), default=None) for p in perms]
    errs = []
    if len(report) != len(orig):
        errs.append("%d report rows, %d dump columns" % (len(report), len(orig)))
    for k, (row, o) in enumerate(zip(report, orig)):
        if o is None:
            if row["EMP1_CT"] != "NA":
                errs.append("%s (column %d): expected NA" % (row["ID"], k + 1))
            continue
        e1 = 2 * float(row["EMP1_CT"])
        lo, hi = count_range([p[k] for p in perms], o, strict)
        if not lo <= e1 <= hi:
            errs.append("%s: EMP1_CT %s, dump gives %g-%g" % (row["ID"], row["EMP1_CT"], lo / 2, hi / 2))
        e2 = 2 * float(row["EMP2_CT"])
        lo, hi = count_range(bests, o, strict)
        if not lo <= e2 <= hi:
            errs.append("%s: EMP2_CT %s, dump gives %g-%g" % (row["ID"], row["EMP2_CT"], lo / 2, hi / 2))
    return errs


def check_adaptive(dump_fname, report_fname, allele_test_ct, max_perms):
    # plink2 defaults: --aperm 5 (i.e. min 6), alpha 0, beta 1e-4, initial
    # interval 1, slope 0.001.
    amin, alpha, beta, intercept, slope = 6, 0.0, 1e-4, 1.0, 0.001 * (1 + 2.0 ** -44)
    zt = NormalDist().inv_cdf(1 - beta / (2.0 * allele_test_ct))
    orig, perms = read_dump(dump_fname)
    perms = perms[:max_perms]
    checks = set()
    c = amin
    while c <= max_perms:
        checks.add(c)
        c += int(intercept + c * slope)
    report = read_report(report_fname)
    errs = []
    for k, (row, o) in enumerate(zip(report, orig)):
        if o is None:
            if row["EMP1_CT"] != "NA":
                errs.append("%s (column %d): expected NA" % (row["ID"], k + 1))
            continue
        # Track the low and high ends of the count range; a stop decision
        # that differs between them can't be checked.
        lo = hi = 0
        stop = None
        for pidx, p in enumerate(perms, 1):
            v = p[k]
            if v is not None:
                if near(v, o):
                    hi += 1
                elif v > o:
                    lo += 2
                    hi += 2
            if pidx in checks:
                decisions = set()
                for ctx2 in (lo, hi):
                    emp1 = (ctx2 + 2) / (2.0 * (pidx + 1))
                    radius = zt * (emp1 * (1 - emp1) / pidx) ** 0.5
                    decisions.add((emp1 - radius > alpha) or (emp1 + radius < alpha))
                if len(decisions) == 2:
                    stop = "ambiguous"
                    break
                if decisions.pop():
                    stop = pidx
                    break
        if stop == "ambiguous":
            continue
        if stop is None:
            stop = max_perms
        e1 = 2 * float(row["EMP1_CT"])
        if int(row["PERM_CT"]) != stop or not lo <= e1 <= hi:
            errs.append("%s: EMP1_CT %s, PERM_CT %s; replay gives %g-%g at %d" % (row["ID"], row["EMP1_CT"], row["PERM_CT"], lo / 2, hi / 2, stop))
    return errs


def check_identity(report_fname, n):
    report = [r for r in read_report(report_fname) if r["EMP1_CT"] != "NA"]
    errs = []
    if not report:
        return ["no valid alleles"]
    half = n / 2
    top_ct = 0
    for r in report:
        if float(r["EMP1_CT"]) != half:
            errs.append("%s: EMP1_CT %s, expected %g" % (r["ID"], r["EMP1_CT"], half))
        e2 = float(r["EMP2_CT"])
        if e2 == half:
            top_ct += 1
        elif e2 != n:
            errs.append("%s: EMP2_CT %s, expected %g or %g" % (r["ID"], r["EMP2_CT"], half, n))
    if top_ct == 0:
        errs.append("no allele has EMP2_CT %g" % half)
    return errs


def main():
    mode = sys.argv[1]
    if mode == "report":
        errs = check_report(sys.argv[2], sys.argv[3])
    elif mode == "ties":
        errs = check_report(sys.argv[2], sys.argv[3], True)
    elif mode == "adaptive":
        errs = check_adaptive(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    elif mode == "identity":
        errs = check_identity(sys.argv[2], int(sys.argv[3]))
    else:
        sys.exit("unknown mode " + mode)
    for e in errs[:10]:
        print(e)
    sys.exit(1 if errs else 0)


main()
