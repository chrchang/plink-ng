#!/usr/bin/env python3

# Usage:
#   check_glm_out.py ci <.glm.* file> <ci size>
#     Checks the L<pct>/U<pct> columns against BETA +/- z * SE (linear) or
#     exp(log(OR) +/- z * LOG(OR)_SE) (logistic/Firth), z = Phi^-1((1+ci)/2).
#   check_glm_out.py corr <psam> <vcf> <variant ID> <covariate>
#     Prints the genotype-covariate correlation and the corresponding
#     two-predictor VIF, 1 / (1 - r^2).
#   check_glm_out.py counts <counts file> <EMP file> <perm ct or 'col'>
#     Checks EMP = (CT + 1) / (PERM_CT + 1) for every EMP1/EMP2 column.
#     PERM_CT is read from the counts file's PERM_CT column if the third
#     argument is 'col'.
#   check_glm_out.py dump <dump.all> <dump.best> <counts .mperm file>
#     Checks that each .dump.best row is the maximum of the corresponding
#     .dump.all row, and recomputes EMP1_CT/EMP2_CT.  (Row 0 holds the
#     original statistics.)

import math
import statistics
import sys


def read_table(fname):
    with open(fname) as f:
        header = f.readline().rstrip("\n").lstrip("#").split("\t")
        return header, [line.rstrip("\n").split("\t") for line in f]


def close(a, b, tol=2e-5):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def check_ci(fname, ci_size):
    header, rows = read_table(fname)
    pct = round(100 * ci_size)
    lcol = header.index(f"L{pct}")
    ucol = header.index(f"U{pct}")
    z = statistics.NormalDist().inv_cdf((1 + ci_size) / 2)
    is_linear = "BETA" in header
    for row in rows:
        if is_linear:
            est = float(row[header.index("BETA")])
            se = float(row[header.index("SE")])
            lo, hi = est - z * se, est + z * se
        else:
            log_or = math.log(float(row[header.index("OR")]))
            se = float(row[header.index("LOG(OR)_SE")])
            lo, hi = math.exp(log_or - z * se), math.exp(log_or + z * se)
        if not (close(float(row[lcol]), lo) and close(float(row[ucol]), hi)):
            sys.exit(f"{fname}: bad CI for {row[header.index('ID')]}: {row[lcol]} {row[ucol]} (expected {lo:g} {hi:g})")


def corr(psam, vcf, vid, covar):
    header, rows = read_table(psam)
    cvals = [float(row[header.index(covar)]) for row in rows]
    with open(vcf) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2 and fields[2] == vid:
                gvals = [{"0/0": 0, "0/1": 1, "1/1": 2}[x] for x in fields[9:]]
    r = statistics.correlation(gvals, cvals)
    print(f"{r:.6f} {1 / (1 - r * r):.6f}")


def counts(ct_fname, emp_fname, perm_ct_arg):
    ct_header, ct_rows = read_table(ct_fname)
    emp_header, emp_rows = read_table(emp_fname)
    assert len(ct_rows) == len(emp_rows)
    for col in ("EMP1", "EMP2"):
        if col not in emp_header:
            continue
        for ct_row, emp_row in zip(ct_rows, emp_rows):
            if perm_ct_arg == "col":
                perm_ct = int(ct_row[ct_header.index("PERM_CT")])
            else:
                perm_ct = int(perm_ct_arg)
            ct = float(ct_row[ct_header.index(col + "_CT")])
            if not close(float(emp_row[emp_header.index(col)]), (ct + 1) / (perm_ct + 1)):
                sys.exit(f"{emp_fname}: {col} mismatch for {emp_row[emp_header.index('ID')]}")


def dump(all_fname, best_fname, ct_fname):
    with open(all_fname) as f:
        all_rows = [[float(x) for x in line.split()[1:]] for line in f]
    with open(best_fname) as f:
        best_rows = [line.split() for line in f]
    assert len(all_rows) == len(best_rows)
    for idx, (all_row, best_row) in enumerate(zip(all_rows, best_rows)):
        assert int(best_row[0]) == idx
        if not close(max(all_row), float(best_row[1]), 1e-6):
            sys.exit(f"{best_fname}: row {idx} is not the maximum")
    orig = all_rows[0]
    perms = all_rows[1:]
    ct_header, ct_rows = read_table(ct_fname)
    assert len(ct_rows) == len(orig)
    for vidx, row in enumerate(ct_rows):
        emp1_ct = sum(p[vidx] >= orig[vidx] for p in perms)
        emp2_ct = sum(max(p) >= orig[vidx] for p in perms)
        if (int(row[ct_header.index("EMP1_CT")]), int(row[ct_header.index("EMP2_CT")])) != (emp1_ct, emp2_ct):
            sys.exit(f"{ct_fname}: counts mismatch for {row[ct_header.index('ID')]}")


def main():
    mode = sys.argv[1]
    if mode == "ci":
        check_ci(sys.argv[2], float(sys.argv[3]))
    elif mode == "corr":
        corr(*sys.argv[2:6])
    elif mode == "counts":
        counts(*sys.argv[2:5])
    elif mode == "dump":
        dump(*sys.argv[2:5])
    else:
        sys.exit("bad mode")


main()
