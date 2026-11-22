# Compares --glm genotypic joint-test (GENO_2DF) p-values against
# ref_joint.txt, the numpy/scipy F-test p-values captured by make_ref.py.  Pure
# Python, so it needs no third-party packages.  As in the original numpy/scipy
# check, each p-value must match to within 1e-4 on the log scale.
import hashlib
import math
import sys

vcf_fn, psam_fn, glm_fn, ref_fn = sys.argv[1:5]

ref = {}
ref_sha256 = None
for line in open(ref_fn):
    if line.startswith("#"):
        continue
    t = line.rstrip("\n").split("\t")
    if t[0] == "input_sha256":
        ref_sha256 = t[1]
    else:
        ref[t[0]] = float(t[1])

input_hash = hashlib.sha256()
for fn in (vcf_fn, psam_fn):
    input_hash.update(open(fn, "rb").read())
assert input_hash.hexdigest() == ref_sha256, (
    "%s/%s differ from the inputs ref_joint.txt was generated from; regenerate "
    "it as described in make_ref.py" % (vcf_fn, psam_fn))

header = None
checked = set()
for line in open(glm_fn):
    t = line.rstrip("\n").split("\t")
    if header is None:
        header = t
        continue
    row = dict(zip(header, t))
    if row["TEST"] != "GENO_2DF":
        continue
    vid = row["ID"]
    p = float(row["P"])
    p_expected = ref[vid]
    assert abs(math.log(p) - math.log(p_expected)) < 1e-4, (vid, p, p_expected)
    checked.add(vid)
assert checked == set(ref), (sorted(checked), sorted(ref))
print("ok")
