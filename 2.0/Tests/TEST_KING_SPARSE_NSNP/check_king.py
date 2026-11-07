# Brute-force check of --make-king-table counts.  Pure Python (no numpy), so
# it runs on CI images without third-party packages: each sample's genotype
# classes are stored as bitmasks over variants, and pairwise counts are
# popcounts of their intersections.
import sys

vcf, kin0 = sys.argv[1], sys.argv[2]
req = set(open(sys.argv[3]).read().split()) if len(sys.argv) > 3 else None

ids = None
cols = None
for line in open(vcf):
    if line.startswith("##"):
        continue
    t = line.rstrip("\n").split("\t")
    if line.startswith("#"):
        ids = t[9:]
        cols = [[] for _ in ids]
        continue
    for i, x in enumerate(t[9:]):
        cols[i].append(x)


def mask(calls, wanted):
    m = 0
    for v, x in enumerate(calls):
        if x in wanted:
            m |= 1 << v
    return m


nm = [mask(c, ("0/0", "0/1", "1/1")) for c in cols]
het = [mask(c, ("0/1",)) for c in cols]
hom0 = [mask(c, ("0/0",)) for c in cols]
hom2 = [mask(c, ("1/1",)) for c in cols]
hom = [a | b for a, b in zip(hom0, hom2)]

try:
    (0).bit_count()

    def popcount(x):
        return x.bit_count()
except AttributeError:  # Python < 3.10

    def popcount(x):
        return bin(x).count("1")


idx = {s: i for i, s in enumerate(ids)}
header = None
seen = set()
bad = 0
for line in open(kin0):
    t = line.split()
    if header is None:
        header = [h.lstrip("#") for h in t]
        continue
    r = dict(zip(header, t))
    i, j = idx[r["IID1"]], idx[r["IID2"]]
    seen.add(frozenset((i, j)))
    want = {
        "NSNP": popcount(nm[i] & nm[j]),
        "HETHET": popcount(het[i] & het[j]),
        "IBS0": popcount((hom0[i] & hom2[j]) | (hom2[i] & hom0[j])),
        "HET1_HOM2": popcount(het[i] & hom[j]),
        "HET2_HOM1": popcount(het[j] & hom[i]),
    }
    for k, w in want.items():
        if int(r[k]) != w:
            bad += 1
            if bad < 4:
                print("mismatch", r["IID1"], r["IID2"], k, r[k], w)
n = len(ids)
expected = set()
for a in range(n):
    for b in range(a):
        if req is not None and ((ids[a] in req) == (ids[b] in req)):
            continue
        expected.add(frozenset((a, b)))
assert seen == expected, (len(seen), len(expected))
assert not bad, bad
print("ok", len(seen))
