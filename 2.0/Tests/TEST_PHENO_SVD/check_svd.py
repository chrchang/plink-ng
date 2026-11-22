# Compares --pheno-svd output against ref_svd.txt, numpy's SVD of the same
# phenotype matrix as captured by make_ref.py.  Pure Python, so it needs no
# third-party packages.  The comparison mirrors the original numpy check:
# np.allclose (rtol=1e-5, atol=1e-5) on each of the first two left/right
# singular vectors, after aligning their sign, plus a 1e-4 relative tolerance
# on the singular values.
import hashlib
import sys

psam_fn, pheno_fn, wts_fn, ref_fn = sys.argv[1:5]

ref = {}
ref_sha256 = None
for line in open(ref_fn):
    if line.startswith("#"):
        continue
    t = line.rstrip("\n").split("\t")
    if t[0] == "input_sha256":
        ref_sha256 = t[1]
        continue
    ref[(t[0], t[1], int(t[2]))] = [float(x) for x in t[3:]]

psam_sha256 = hashlib.sha256(open(psam_fn, "rb").read()).hexdigest()
assert psam_sha256 == ref_sha256, (
    "%s differs from the input ref_svd.txt was generated from; regenerate it "
    "as described in make_ref.py" % psam_fn)

pheno_cols = [[], []]
for line in open(pheno_fn):
    if line.startswith("#"):
        continue
    t = line.split()
    pheno_cols[0].append(float(t[1]))
    pheno_cols[1].append(float(t[2]))
wts_lines = [line.split() for line in open(wts_fn) if not line.startswith("#")]
svs = [float(t[1]) for t in wts_lines]
vt = [[float(x) for x in t[2:]] for t in wts_lines]


def allclose(a, b, rtol=1e-5, atol=1e-5):
    return len(a) == len(b) and all(abs(x - y) <= atol + rtol * abs(y) for x, y in zip(a, b))


def matches(name):
    for k in range(2):
        u = ref[(name, "u", k)]
        vh = ref[(name, "vh", k)]
        (s,) = ref[(name, "s", k)]
        sign = 1.0 if sum(x * y for x, y in zip(u, pheno_cols[k])) >= 0 else -1.0
        if not allclose(pheno_cols[k], [sign * x for x in u]):
            return False
        if not allclose(vt[k], [sign * x for x in vh]):
            return False
        if abs(svs[k] - s) > 1e-4 * s:
            return False
    return True


assert matches("raw") or matches("centered"), "SVD mismatch"
print("ok")
