# Usage: compare_glm.py <glm file A> <glm file B> [old=new ...]
# Compares (ID, TEST, effect, SE, P) rows after renaming TEST substrings in A,
# with rows sorted within each file.  Local covariates are stored in single
# precision, so a small tolerance is allowed.
import math
import sys

def load(fn, ren):
    header = None
    rows = []
    for line in open(fn):
        t = line.rstrip("\n").split("\t")
        if header is None:
            header = {k: i for i, k in enumerate(t)}
            continue
        test = t[header["TEST"]]
        for k, v in ren.items():
            test = test.replace(k, v)
        cols = [c for c in ("BETA", "OR", "SE", "LOG(OR)_SE", "P") if c in header]
        rows.append((t[header["ID"]], test, [t[header[c]] for c in cols]))
    return sorted(rows)

a = load(sys.argv[1], dict(x.split("=") for x in sys.argv[3:]))
b = load(sys.argv[2], {})
assert len(a) == len(b), (len(a), len(b))
for (ida, ta, va), (idb, tb, vb) in zip(a, b):
    assert (ida, ta) == (idb, tb), (ida, ta, idb, tb)
    for x, y in zip(va, vb):
        assert x != "NA" and y != "NA", (ida, ta, va, vb)
        x, y = float(x), float(y)
        assert abs(x - y) <= 1e-4 * max(1.0, abs(y)) or abs(math.log(x) - math.log(y)) < 1e-3, (ida, ta, va, vb)
print("ok", len(a))
