import sys
import numpy as np
from scipy import stats

vcf_fn, psam_fn, glm_fn = sys.argv[1:4]
genos = {}
for line in open(vcf_fn):
    if line.startswith("#"):
        continue
    t = line.split()
    genos[t[2]] = np.array([{"0/0": 0, "0/1": 1, "1/1": 2}[x] for x in t[9:]], dtype=float)
ph = np.genfromtxt(psam_fn, comments=None, skip_header=1, usecols=(2, 3, 4))
y, covs = ph[:, 0], ph[:, 1:]

def rss(x):
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    r = y - x @ beta
    return r @ r

header = None
checked = 0
for line in open(glm_fn):
    t = line.rstrip("\n").split("\t")
    if header is None:
        header = t
        continue
    row = dict(zip(header, t))
    if row["TEST"] != "GENO_2DF":
        continue
    g = genos[row["ID"]]
    n = len(y)
    one = np.ones(n)
    full = np.column_stack([one, g, (g == 1).astype(float), covs])
    reduced = np.column_stack([one, covs])
    rss1, rss0 = rss(full), rss(reduced)
    df2 = n - full.shape[1]
    fstat = ((rss0 - rss1) / 2) / (rss1 / df2)
    p_expected = stats.f.sf(fstat, 2, df2)
    p = float(row["P"])
    assert abs(np.log(p) - np.log(p_expected)) < 1e-4, (row["ID"], p, p_expected)
    checked += 1
assert checked == 3, checked
print("ok")
