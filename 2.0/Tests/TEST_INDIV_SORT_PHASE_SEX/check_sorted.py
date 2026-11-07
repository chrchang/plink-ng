import sys

orig_fn, psam_fn, sorted_fn, shm_fn = sys.argv[1:5]

def load(fn):
    calls = {}
    for line in open(fn):
        if line.startswith("##"):
            continue
        t = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            ids = t[9:]
            continue
        calls[t[2]] = dict(zip(ids, t[9:]))
    return calls

sex = {}
for line in open(psam_fn):
    if not line.startswith("#"):
        iid, s = line.split()
        sex[iid] = s
orig = load(orig_fn)
for fn, set_hh in ((sorted_fn, False), (shm_fn, True)):
    got = load(fn)
    assert got.keys() == orig.keys()
    for vid in orig:
        for iid, gt in orig[vid].items():
            expected = gt
            if vid.startswith("x") and sex[iid] == "1" and len(gt) == 3:
                # A male diploid call on chrX: hets are invalid, homs are
                # written as haploid.
                if gt[0] != gt[2]:
                    expected = "." if set_hh else gt
                else:
                    expected = gt[0] if set_hh else gt
            g = got[vid][iid]
            if not set_hh and vid.startswith("x") and sex[iid] == "1" and len(gt) == 3 and gt[0] == gt[2]:
                # Without --set-invalid-haploid-missing, male homs may be
                # written as haploid too.
                ok = g in (gt, gt[0])
            elif expected == ".":
                ok = g in (".", "./.")
            else:
                ok = (g == expected)
            assert ok, (fn, vid, iid, gt, g)
print("ok")
