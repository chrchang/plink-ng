# Usage: compare.py [--loose] <glm file A> <glm file B (dense reference)>
# File B's run coded the indicator covariates (BATCH_*, CENTER_*, FLAG_*) as
# -1 = in category, 0 = not in category.  --glm only groups two-valued
# columns whose larger values mark disjoint sample sets, so these take the
# generic dense path.  File A's run used categorical covariates (TEST names
# "X=level", mapped to "X_level") or 0/1 or 1/2 indicators, which are
# recognized as groups.  After internal rescaling, the 0/-1 coding is the
# exact negation of the 0/1 coding (and of the 1/2 coding when centering is
# in effect), so for indicator terms, including interactions with them, BETA
# and Z_STAT flip sign and the CI bounds swap; SE and P are unchanged.  All
# other rows must match.  Values must agree to 6-significant-digit output
# precision; ERRCODE exactly.  --loose (for single-precision runs, where the
# fast path accumulates in double precision but the dense reference doesn't)
# instead allows coefficient/CI differences up to 1e-4 SEs, and relative
# SE/Z/P differences up to 1e-4.
import math
import sys

INDICATOR_PREFIXES = ("BATCH_", "CENTER_", "FLAG_")

def load(fn):
    header = None
    rows = {}
    for line in open(fn):
        t = line.rstrip("\n").split("\t")
        if header is None:
            header = t
            continue
        rec = dict(zip(header, t))
        test = rec["TEST"].replace("=", "_")
        key = (rec["#CHROM"], rec["ID"], rec["A1"], test)
        assert key not in rows, key
        rows[key] = rec
    return header, rows

args = sys.argv[1:]
loose = (args[0] == "--loose")
if loose:
    args = args[1:]
ha, a = load(args[0])
hb, b = load(args[1])
assert ha == hb, (ha, hb)
assert sorted(a) == sorted(b), sorted(set(a) ^ set(b))[:10]
is_or = "OR" in ha
cols = [c for c in ha if c in ("A1_FREQ", "OBS_CT", "OR", "BETA", "LOG(OR)_SE", "SE", "L95", "U95", "Z_STAT", "P")]

def close(x, y, scale=None):
    if loose:
        if scale is not None:
            return abs(x - y) <= 1e-4 * scale
        return abs(x - y) <= 1e-4 * max(abs(y), 1e-300) + 1e-300
    return abs(x - y) <= 2e-5 * max(abs(y), 1e-300) + 1e-300

nonna = 0
for key in a:
    ra = a[key]
    rb = b[key]
    assert ra["ERRCODE"] == rb["ERRCODE"], (key, ra["ERRCODE"], rb["ERRCODE"])
    flip = any(p in key[3] for p in INDICATOR_PREFIXES)
    se_col = "SE" if "SE" in rb else "LOG(OR)_SE"
    se = float(rb[se_col]) if rb[se_col] != "NA" else None
    for col in cols:
        src = col
        if flip and col in ("L95", "U95"):
            src = "U95" if col == "L95" else "L95"
        x = ra[col]
        y = rb[src]
        if x == "NA" or y == "NA":
            assert x == y, (key, col, x, y)
            continue
        x = float(x)
        y = float(y)
        if is_or and col in ("OR", "L95", "U95"):
            lx = math.log(x)
            ly = -math.log(y) if flip else math.log(y)
            if loose:
                assert abs(lx - ly) <= 1e-4 * se + 2e-5 * max(1.0, abs(ly)), (key, col, x, y)
            else:
                assert abs(lx - ly) <= 2e-5 * max(1.0, abs(ly)), (key, col, x, y)
        else:
            if flip and col in ("BETA", "L95", "U95", "Z_STAT"):
                y = -y
            if loose and col in ("BETA", "L95", "U95"):
                assert close(x, y, se), (key, col, x, y)
            elif loose and col == "P":
                assert close(math.log(x), math.log(y)), (key, col, x, y)
            else:
                assert close(x, y), (key, col, x, y)
    nonna += (ra["P"] != "NA")
assert nonna > len(a) // 2, (nonna, len(a))
print("ok", len(a), "rows,", nonna, "non-NA")
