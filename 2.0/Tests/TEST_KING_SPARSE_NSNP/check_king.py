import sys
import numpy as np
vcf, kin0 = sys.argv[1], sys.argv[2]
req = set(open(sys.argv[3]).read().split()) if len(sys.argv) > 3 else None
code = {"0/0": 0, "0/1": 1, "1/1": 2, "./.": -1}
rows = []
for line in open(vcf):
    if line.startswith("##"):
        continue
    t = line.rstrip("\n").split("\t")
    if line.startswith("#"):
        ids = t[9:]
        continue
    rows.append([code[x] for x in t[9:]])
g = np.array(rows)
nm = (g >= 0).astype(np.int64)
het = (g == 1).astype(np.int64)
hom0 = (g == 0).astype(np.int64)
hom2 = (g == 2).astype(np.int64)
nsnp = nm.T @ nm
hethet = het.T @ het
ibs0 = hom0.T @ hom2 + hom2.T @ hom0
het_hom = het.T @ (hom0 + hom2)   # [i, j]: i het, j hom
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
    want = {"NSNP": nsnp[i, j], "HETHET": hethet[i, j], "IBS0": ibs0[i, j], "HET1_HOM2": het_hom[i, j], "HET2_HOM1": het_hom[j, i]}
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
