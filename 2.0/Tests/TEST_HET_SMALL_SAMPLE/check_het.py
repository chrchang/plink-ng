# Expected E(HOM) for --het small-sample on fully-called data: the sum over
# founder-polymorphic variants of 1 - 2pq * n/(n-1), where p is the founder
# ALT frequency and n is the number of founder alleles.  Founders are the
# samples whose IDs start with 'f'.
import sys

vcf_fn, het_fn = sys.argv[1:3]
rows = []
for line in open(vcf_fn):
    if line.startswith("##"):
        continue
    t = line.rstrip("\n").split("\t")
    if line.startswith("#"):
        ids = t[9:]
        continue
    rows.append([{"0/0": 0, "0/1": 1, "1/1": 2}[x] for x in t[9:]])
founders = [i for i, s in enumerate(ids) if s.startswith("f")]
n = 2 * len(founders)
expected = 0.0
for r in rows:
    p = sum(r[i] for i in founders) / n
    if 0.0 < p < 1.0:
        expected += 1 - 2 * p * (1 - p) * n / (n - 1)
header = None
checked = 0
for line in open(het_fn):
    t = line.split()
    if header is None:
        header = t
        continue
    got = float(t[header.index("E(HOM)")])
    assert abs(got - expected) < 0.0015, (t[0], got, expected)
    checked += 1
assert checked == len(ids), checked
print("ok", round(expected, 3))
