#!/bin/bash

# KING-robust matrix and sparse-GRM flags with no other test coverage:
# --make-king (all shape/format modifiers), --king-cutoff (from a .king.bin
# fileset and computed on the fly), --king-table-require, and
# --make-grm-sparse.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Founders A-D; E and F are full sibs (children of A and B); G is a noisy
# duplicate of C.
python3 - <<'PY'
import random
random.seed(7)
n_var = 400
def founder():
    return [[random.random() < 0.3 for _ in range(2)] for _ in range(n_var)]
g = {k: founder() for k in "ABCD"}
def child(p, q):
    return [[random.choice(p[i]), random.choice(q[i])] for i in range(n_var)]
g["E"] = child(g["A"], g["B"])
g["F"] = child(g["A"], g["B"])
g["G"] = [list(x) if random.random() > 0.02 else [random.random() < 0.3, random.random() < 0.3] for x in g["C"]]
ids = list("ABCDEFG")
with open("tmp_data.vcf", "w") as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=10000000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for i in range(n_var):
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t" % (1000 * (i + 1), i) + "\t".join("%d/%d" % tuple(g[s][i]) for s in ids) + "\n")
PY
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null
$plink2 --pfile tmp_data --make-king-table --out tmp_table > /dev/null

# Sanity-check the relationships themselves: parent-child and full-sib pairs
# near 0.25, the duplicate near 0.5, and unrelated pairs near 0.
awk -F '\t' 'NR > 1 {
    pair = $1 $2
    rel = (pair ~ /^(EA|EB|FA|FB|FE)$/)? 0.25 : ((pair == "GC")? 0.5 : 0)
    if ($6 < rel - 0.1 || $6 > rel + 0.1) { print "unexpected kinship " $0; exit 1 }
}' tmp_table.kin0

# --make-king text shapes.
$plink2 --pfile tmp_data --make-king --out tmp_tri > /dev/null
test "$(tail -n +2 tmp_tri.king.id | tr '\n' ' ')" = "A B C D E F G "
python3 check_king.py king tmp_table.kin0 tmp_tri.king.id triangle tmp_tri.king
$plink2 --pfile tmp_data --make-king triangle --out tmp_tri2 > /dev/null
cmp tmp_tri.king tmp_tri2.king
$plink2 --pfile tmp_data --make-king square --out tmp_sq > /dev/null
python3 check_king.py king tmp_table.kin0 tmp_sq.king.id square tmp_sq.king
$plink2 --pfile tmp_data --make-king square0 --out tmp_sq0 > /dev/null
python3 check_king.py king tmp_table.kin0 tmp_sq0.king.id square0 tmp_sq0.king

# Binary shapes.  The default binary shape is square.
$plink2 --pfile tmp_data --make-king bin --out tmp_bsq > /dev/null
python3 check_king.py king tmp_table.kin0 tmp_bsq.king.id bin8sq tmp_bsq.king.bin
$plink2 --pfile tmp_data --make-king bin triangle --out tmp_btri > /dev/null
python3 check_king.py king tmp_table.kin0 tmp_btri.king.id bin8tri tmp_btri.king.bin
$plink2 --pfile tmp_data --make-king bin4 triangle --out tmp_b4tri > /dev/null
python3 check_king.py king tmp_table.kin0 tmp_b4tri.king.id bin4tri tmp_b4tri.king.bin
# 'zs' compresses the text matrix.
$plink2 --pfile tmp_data --make-king zs --out tmp_zs > /dev/null
test -s tmp_zs.king.zst
test ! -e tmp_zs.king

# --king-cutoff: kinship > 0.177 links E-A, E-B, F-A, F-B, F-E, and G-C.
# Every input route must make the same choice.
for route in "--king-cutoff 0.177" "--king-cutoff tmp_btri 0.177" "--king-cutoff tmp_b4tri 0.177" "--king-cutoff-table tmp_table.kin0 0.177"; do
    $plink2 --pfile tmp_data $route --out tmp_cut > /dev/null
    test "$(tail -n +2 tmp_cut.king.cutoff.in.id | tr '\n' ' ')" = "A B C D "
    test "$(tail -n +2 tmp_cut.king.cutoff.out.id | tr '\n' ' ')" = "E F G "
done
# The exclusions apply to the rest of the run.
$plink2 --pfile tmp_data --king-cutoff tmp_btri 0.177 --make-just-psam --out tmp_cut_psam > /dev/null
test "$(tail -n +2 tmp_cut_psam.psam | cut -f 1 | tr '\n' ' ')" = "A B C D "
# A looser cutoff only separates the duplicates.
$plink2 --pfile tmp_data --king-cutoff tmp_btri 0.35 --out tmp_cut2 > /dev/null
test "$(tail -n +2 tmp_cut2.king.cutoff.out.id | wc -l)" -eq 1
grep -qx -e C -e G tmp_cut2.king.cutoff.out.id
# A square .bin is rejected.
if $plink2 --pfile tmp_data --king-cutoff tmp_bsq 0.177 --out tmp_cut_sq > /dev/null 2>&1; then
    exit 1
fi

# --king-table-require: only pairs involving E or C.
printf 'E\nC\n' > tmp_require.txt
$plink2 --pfile tmp_data --make-king-table --king-table-require tmp_require.txt --out tmp_req > /dev/null
awk -F '\t' 'NR > 1 && ($1 == "E" || $2 == "E" || $1 == "C" || $2 == "C")' tmp_table.kin0 > tmp_req_expected.txt
tail -n +2 tmp_req.kin0 | diff - tmp_req_expected.txt
test "$(wc -l < tmp_req_expected.txt)" -eq 11
# --king-table-require-xor wants exactly one of them, which drops E-C.
$plink2 --pfile tmp_data --make-king-table --king-table-require-xor tmp_require.txt --out tmp_xor > /dev/null
test "$(tail -n +2 tmp_xor.kin0 | wc -l)" -eq 10
if grep -q "$(printf '^E\tC\t')" tmp_xor.kin0; then
    exit 1
fi

# --make-grm-sparse: exactly the --make-grm-list entries at or above the
# threshold, 0-based.  Seven samples are too few to estimate allele
# frequencies, so read them in, and drop monomorphic variants.
$plink2 --pfile tmp_data --freq --out tmp_freq > /dev/null
grm_flags="--read-freq tmp_freq.afreq --maf 0.01"
$plink2 --pfile tmp_data $grm_flags --make-grm-list --out tmp_list > /dev/null
$plink2 --pfile tmp_data $grm_flags --make-grm-list cov --out tmp_list_cov > /dev/null
for threshold in 0.2 -0.1 2; do
    $plink2 --pfile tmp_data $grm_flags --make-grm-sparse $threshold --out tmp_sp > /dev/null
    python3 check_king.py sparse tmp_list.grm tmp_sp.grm.sp $threshold
    $plink2 --pfile tmp_data $grm_flags --make-grm-sparse $threshold cov --out tmp_sp_cov > /dev/null
    python3 check_king.py sparse tmp_list_cov.grm tmp_sp_cov.grm.sp $threshold
done
test "$(wc -l < tmp_sp.grm.sp)" -eq 0
# The 0.2 run keeps the diagonal plus the G-C duplicate pair.
$plink2 --pfile tmp_data $grm_flags --make-grm-sparse 0.2 --out tmp_sp2 > /dev/null
test "$(awk -F '\t' '$1 != $2' tmp_sp2.grm.sp | cut -f 1,2)" = "$(printf '6\t2')"
# ID file variants.
test "$(head -n 1 tmp_sp2.grm.id)" = "$(printf '0\tA')"
$plink2 --pfile tmp_data $grm_flags --make-grm-sparse 0.2 iid-only --out tmp_sp_iid > /dev/null
test "$(tr '\n' ' ' < tmp_sp_iid.grm.id)" = "A B C D E F G "
$plink2 --pfile tmp_data $grm_flags --make-grm-sparse 0.2 id-header --out tmp_sp_hdr > /dev/null
test "$(head -n 2 tmp_sp_hdr.grm.id | tr '\n' ' ')" = "#IID A "
cmp tmp_sp2.grm.sp tmp_sp_hdr.grm.sp
$plink2 --pfile tmp_data $grm_flags --make-grm-sparse 0.2 zs --out tmp_sp_zs > /dev/null
test -s tmp_sp_zs.grm.sp.zst
