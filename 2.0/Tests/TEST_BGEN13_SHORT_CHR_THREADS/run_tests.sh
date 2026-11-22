#!/bin/bash

# BGEN v1.3 import when every selected variant fits in the first scan block
# (here, --chr picks a chromosome with one or two variants).
#
# In that case the first scan job is only launched after the read loop, and
# the old code gave a variant range to just the next worker thread; any
# further threads ran on uninitialized ranges, typically failing with a
# ZSTD_decompress() error.  The imported dosages must match a
# single-threaded import.
#
# The scan uses (workspace size / 4 GiB) threads, so several scan threads,
# and thus reproducing the old failure, needs a large --memory reservation.
# That part is skipped if the reservation can't be made.

set -exo pipefail

plink2="$1/plink2 $2"

python3 - <<'PY'
import random
random.seed(13)
n = 3000
ids = ["s%d" % i for i in range(n)]
with open("tmp_data.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c in (1, 2, 3):
        f.write("##contig=<ID=%d,length=10000000>\n" % c)
    f.write('##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    # chr1: 50 variants; chr2: 1 variant; chr3: 7 variants.
    for c, nv in ((1, 50), (2, 1), (3, 7)):
        for v in range(nv):
            f.write("%d\t%d\tc%dv%d\tA\tG\t.\t.\t.\tDS\t%s\n" % (c, 1000 * (v + 1), c, v, "\t".join("%.3f" % (random.random() * 2) for _ in ids)))
PY
$plink2 --vcf tmp_data.vcf dosage=DS --make-pgen --out tmp_data > /dev/null
# plink2's .bgen export stores ALT first.
$plink2 --pfile tmp_data --export bgen-1.3 --out tmp_data > /dev/null
import_args="--bgen tmp_data.bgen ref-last --sample tmp_data.sample"

for c in 2 3; do
    $plink2 $import_args --threads 1 --chr $c --make-pgen --out tmp_ref_chr$c > /dev/null
    $plink2 --pfile tmp_ref_chr$c --export A --out tmp_ref_chr$c > /dev/null
    # Sanity check: the reference matches the source dosages (BGEN stores
    # them with 16-bit precision).
    $plink2 --pfile tmp_data --chr $c --export A --out tmp_src_chr$c > /dev/null
    python3 -c "
def load(fn):
    return [line.split()[6:] for line in open(fn)][1:]
a, b = load('tmp_ref_chr$c.raw'), load('tmp_src_chr$c.raw')
assert len(a) == len(b) == 3000
for ra, rb in zip(a, b):
    for x, y in zip(ra, rb):
        assert abs(float(x) - float(y)) < 2e-4, (x, y)
"

    # 1. Default workspace size.
    for t in 2 4 8; do
        $plink2 $import_args --threads $t --chr $c --make-pgen --out tmp_t${t}_chr$c > /dev/null
        $plink2 --pfile tmp_t${t}_chr$c --export A --out tmp_t${t}_chr$c > /dev/null
        cmp tmp_ref_chr$c.raw tmp_t${t}_chr$c.raw || exit 1
    done

    # 2. Large workspace, so that up to 7 scan threads are used.  32-bit
    #    builds reject --memory 30000 outright; since they can't reserve even
    #    one 4 GiB scan-thread's worth of workspace, there is nothing for this
    #    part to check there.
    for t in 4 8 16; do
        if ! $plink2 $import_args --threads $t --memory 30000 --chr $c --make-pgen --out tmp_big_t${t}_chr$c > tmp_big_t${t}_chr$c.out 2>&1; then
            if grep -q 'argument too large for 32-bit version' tmp_big_t${t}_chr$c.out; then
                echo "32-bit build; skipping large-workspace check."
                break
            fi
            cat tmp_big_t${t}_chr$c.out
            exit 1
        fi
        if grep -q 'after larger attempt(s) failed' tmp_big_t${t}_chr$c.log; then
            echo "--memory 30000 reservation unavailable; skipping large-workspace check."
            break
        fi
        $plink2 --pfile tmp_big_t${t}_chr$c --export A --out tmp_big_t${t}_chr$c > /dev/null
        cmp tmp_ref_chr$c.raw tmp_big_t${t}_chr$c.raw || exit 1
    done
done
