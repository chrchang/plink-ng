import sys
import numpy as np

psam_fn, pheno_fn, wts_fn = sys.argv[1:4]
psam = np.genfromtxt(psam_fn, comments=None, skip_header=1, usecols=(2, 3, 4))
# Compare against both the raw and the column-centered matrix, since only the
# layout of the result is under test here.
pheno = np.genfromtxt(pheno_fn, comments=None, skip_header=1, usecols=(1, 2))
wts_lines = [line.split() for line in open(wts_fn) if not line.startswith("#")]
svs = np.array([float(t[1]) for t in wts_lines])
vt = np.array([[float(x) for x in t[2:]] for t in wts_lines])

def matches(mat):
    u, s, vh = np.linalg.svd(mat, full_matrices=False)
    for k in range(2):
        sign = np.sign(u[:, k] @ pheno[:, k])
        if not np.allclose(pheno[:, k], sign * u[:, k], atol=1e-5):
            return False
        if not np.allclose(vt[k], sign * vh[k], atol=1e-5):
            return False
        if abs(svs[k] - s[k]) > 1e-4 * s[k]:
            return False
    return True

assert matches(psam) or matches(psam - psam.mean(axis=0)), "SVD mismatch"
print("ok")
