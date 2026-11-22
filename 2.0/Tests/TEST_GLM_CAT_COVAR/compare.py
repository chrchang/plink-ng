# Usage: compare.py <glm file from categorical covariates> <glm file from
#   pre-expanded 0/1 covariates>
# TEST names "X=level" in the first file correspond to "X_level" in the second.
# Every row must match: ERRCODE exactly, numeric columns to 6-significant-digit
# output precision.
import sys

def load(fn, rename):
    header = None
    rows = {}
    for line in open(fn):
        t = line.rstrip("\n").split("\t")
        if header is None:
            header = t
            continue
        rec = dict(zip(header, t))
        test = rec["TEST"].replace("=", "_") if rename else rec["TEST"]
        key = (rec["#CHROM"], rec["ID"], rec["A1"], test)
        assert key not in rows, key
        rows[key] = rec
    return header, rows

ha, a = load(sys.argv[1], True)
hb, b = load(sys.argv[2], False)
assert ha == hb, (ha, hb)
assert sorted(a) == sorted(b), sorted(set(a) ^ set(b))[:10]
numeric = [c for c in ha if c in ("A1_FREQ", "OBS_CT", "OR", "BETA", "LOG(OR)_SE", "SE", "L95", "U95", "Z_STAT", "P")]
nonna = 0
for key in a:
    ra = a[key]
    rb = b[key]
    assert ra["ERRCODE"] == rb["ERRCODE"], (key, ra["ERRCODE"], rb["ERRCODE"])
    for col in numeric:
        x = ra[col]
        y = rb[col]
        if x == "NA" or y == "NA":
            assert x == y, (key, col, x, y)
            continue
        x = float(x)
        y = float(y)
        assert abs(x - y) <= 2e-5 * max(abs(y), 1e-300) + 1e-300, (key, col, x, y)
    nonna += (ra["P"] != "NA")
assert nonna > len(a) // 2, (nonna, len(a))
print("ok", len(a), "rows,", nonna, "non-NA")
