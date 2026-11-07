#!/bin/bash

# --mendel on chrY with non-trio samples present, and per-sample
# observation counts.
#
# 1. The female mask (used on chrY, where females' calls are set to
#    missing) was built by masking a collapsed (trio-sample) bitarray with
#    the uncollapsed male bitarray, so trio males could be treated as
#    female because of the samples at the same raw positions.
# 2. The per-trio missing-call accumulator was flushed every 255 affected
#    variants but never reset, so .imendel's OBS_CT was wrong whenever one
#    thread processed more than 255 such variants.  The result must not
#    depend on the thread count.

set -exo pipefail

plink2="$1/plink2 $2"

# 1. Three unrelated females precede father P / mother M / son K.  P is
#    REF and K is ALT on y1: one Mendel error.
printf '#IID\tPAT\tMAT\tSEX\nU1\t0\t0\t2\nU2\t0\t0\t2\nU3\t0\t0\t2\nP\t0\t0\t1\nM\t0\t0\t2\nK\tP\tM\t1\n' > tmp_y.psam
cat > tmp_y.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=Y,length=50000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	U1	U2	U3	P	M	K
Y	5000000	y1	A	G	.	.	.	GT	.	.	.	0	.	1
VCF
$plink2 --vcf tmp_y.vcf --psam tmp_y.psam --make-pgen --out tmp_y > /dev/null
$plink2 --pfile tmp_y --mendel --out tmp_y > /dev/null
[[ "$(grep -vc '^#' tmp_y.mendel)" -eq 1 ]] || exit 1
grep -v '^#' tmp_y.mendel | grep -q $'\ty1\t'

# 2. 3000 variants with no errors.  K is missing on v % 4 == 0, P on
#    v % 3 == 1, M on v % 11 == 5.  A variant counts for a trio member
#    unless the child, or both parents, or that member is missing.
printf '#IID\tPAT\tMAT\tSEX\nP\t0\t0\t1\nM\t0\t0\t2\nK\tP\tM\t1\n' > tmp_m.psam
python3 - <<'PY'
with open("tmp_m.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\tM\tK\n")
    obs = {"P": 0, "M": 0, "K": 0}
    for v in range(3000):
        p_miss, m_miss, k_miss = (v % 3 == 1), (v % 11 == 5), (v % 4 == 0)
        gts = ["./." if x else "0/0" for x in (p_miss, m_miss, k_miss)]
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 * (v + 1), v, "\t".join(gts)))
        if not (k_miss or (p_miss and m_miss)):
            obs["K"] += 1
            obs["P"] += not p_miss
            obs["M"] += not m_miss
with open("tmp_expected_obs.txt", "w") as f:
    f.write("P:0:%d M:0:%d K:0:%d \n" % (obs["P"], obs["M"], obs["K"]))
PY
$plink2 --vcf tmp_m.vcf --psam tmp_m.psam --make-pgen --out tmp_m > /dev/null
for t in 1 2 8; do
    $plink2 --pfile tmp_m --threads $t --mendel summaries-only cols=+nobsi --out tmp_m_t$t > /dev/null
    [[ "$(tail -n +2 tmp_m_t$t.imendel | tr '\n\t' ' :')" == "$(cat tmp_expected_obs.txt)" ]] || exit 1
done
