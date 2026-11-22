# Usage: compare_rescaled.py <auto glm file> <reference glm file> [TEST=mult ...]
# The reference run used covariates that were already on a safe scale, with
# covariate "X" there equal to (original X) * mult, possibly shifted.  So
# reported effects for the original covariate (and its interaction terms) must
# equal the reference's times mult (on the log scale for OR), and everything
# else (genotype terms, P-values, joint tests) must match directly.
import math
import sys

def load(fn):
    header = None
    rows = {}
    for line in open(fn):
        t = line.rstrip("\n").split("\t")
        if header is None:
            header = t
            continue
        rec = dict(zip(header, t))
        key = (rec["#CHROM"], rec["ID"], rec["A1"], rec["TEST"])
        assert key not in rows, key
        rows[key] = rec
    return header, rows

ha, a = load(sys.argv[1])
hb, b = load(sys.argv[2])
mults = {}
for arg in sys.argv[3:]:
    k, v = arg.split("=")
    mults[k] = float(v)
assert ha == hb, (ha, hb)
assert sorted(a) == sorted(b), (sorted(a), sorted(b))

# Values are printed with 6 significant digits, so after taking logs of
# OR-scale values, an absolute tolerance is needed.
def close(x, y, tol, abs_tol=0.0):
    return abs(x - y) <= tol * max(1e-3, abs(y)) + abs_tol

tol = 2e-5
for key in sorted(a):
    ra = a[key]
    rb = b[key]
    test = key[3]
    mult = 1.0
    for k, v in mults.items():
        if test == k or test.endswith("x" + k):
            mult = v
    for col in ha:
        va = ra[col]
        vb = rb[col]
        if col in ("BETA", "SE", "LOG(OR)_SE", "OR", "L95", "U95"):
            if va == "NA" or vb == "NA":
                assert va == vb, (key, col, va, vb)
                continue
            x = float(va)
            y = float(vb)
            abs_tol = 0.0
            if col == "OR" or ((col in ("L95", "U95")) and ("OR" in ha)):
                x = math.log(x)
                y = math.log(y)
                abs_tol = 1e-5
            assert close(x, y * mult, tol, abs_tol), (key, col, va, vb, mult)
        elif col in ("P", "T_STAT", "Z_STAT", "LOG10_P"):
            if va == "NA" or vb == "NA":
                assert va == vb, (key, col, va, vb)
                continue
            x = float(va)
            y = float(vb)
            if col == "P":
                assert close(math.log(x), math.log(y), tol) or close(x, y, tol), (key, col, va, vb)
            else:
                assert close(x, y, tol), (key, col, va, vb)
        else:
            assert va == vb, (key, col, va, vb)
print("ok", len(a))
