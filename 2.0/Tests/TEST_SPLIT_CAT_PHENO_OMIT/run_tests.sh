#!/bin/bash

# --split-cat-pheno omit-most/omit-last when the omitted category has the
# highest index among the observed categories.
#
# The per-category pointer array was sized by the largest *kept* category
# index, but was also written at the omitted category's index.  That
# out-of-bounds write usually lands in allocation slack, so this mainly
# guards the output (and sanitizer builds).

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
n = 200
ids = ["S%d" % i for i in range(n)]
with open("tmp_data.psam", "w") as f:
    f.write("#IID\tSEX\tCAT\n")
    for i, s in enumerate(ids):
        f.write("%s\t%d\tc%03d\n" % (s, 1 + i % 2, 2 if i in (50, 51) else i % 100))
with open("tmp_data.vcf", "w") as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    f.write("1\t1\tv1\tA\tG\t.\t.\t.\tGT\t" + "\t".join(["0/1"] * n) + "\n")
# Keep categories c001 (2 samples), c002 (4 samples), c099 (2 samples).
with open("tmp_keep.txt", "w") as f:
    for i, s in enumerate(ids):
        if i % 100 in (1, 2, 99) or i in (50, 51):
            f.write(s + "\n")
PY
$plink2 --vcf tmp_data.vcf --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

# omit-last drops c099; omit-most drops c002.
$plink2 --pfile tmp_data --keep tmp_keep.txt --split-cat-pheno omit-last --make-just-psam --out tmp_last > /dev/null
[[ "$(head -n 1 tmp_last.psam)" == "$(printf '#IID\tSEX\tCAT=c001\tCAT=c002')" ]]
$plink2 --pfile tmp_data --keep tmp_keep.txt --split-cat-pheno omit-most --make-just-psam --out tmp_most > /dev/null
[[ "$(head -n 1 tmp_most.psam)" == "$(printf '#IID\tSEX\tCAT=c001\tCAT=c099')" ]]
python3 - <<'PY'
cats = {}
for line in open("tmp_data.psam"):
    if not line.startswith("#"):
        iid, _, cat = line.split()
        cats[iid] = cat
for fn, kept in (("tmp_last.psam", ("c001", "c002")), ("tmp_most.psam", ("c001", "c099"))):
    rows = [line.split() for line in open(fn)]
    assert len(rows) == 9, (fn, len(rows))
    for r in rows[1:]:
        want = ["2" if cats[r[0]] == c else "1" for c in kept]
        assert r[2:] == want, (fn, r, want)
print("ok")
PY
