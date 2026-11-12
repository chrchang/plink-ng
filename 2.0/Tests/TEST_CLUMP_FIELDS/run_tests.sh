#!/bin/bash

# --clump options with no other test coverage: --clump-a1-field,
# --clump-test-field, --clump-test, --clump-log10 (all modes),
# --clump-log10-p1/-p2, --clump-bins, and --clump-unphased.
#
# Also a regression test: --clump-bins never enforced increasing order (the
# previous value wasn't tracked), so unsorted boundaries were accepted and
# produced wrong bin counts.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Phased data: v1-v4 in strong LD (r^2 > 0.8), v5 independent, v6-v8 in
# moderate LD (r^2 ~0.6).
python3 - <<'PY'
import random
random.seed(11)
n = 200
def hap_block(k, noise):
    base = [random.random() < 0.4 for _ in range(2 * n)]
    return [[b if random.random() > noise else (random.random() < 0.4) for b in base] for _ in range(k)]
haps = hap_block(4, 0.05) + [[random.random() < 0.3 for _ in range(2 * n)]] + hap_block(3, 0.1)
with open("tmp_data.vcf", "w") as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=1000000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("s%d" % i for i in range(n)) + "\n")
    for v, h in enumerate(haps):
        f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t" % (10000 * (v + 1), v + 1) + "\t".join("%d|%d" % (h[2*i], h[2*i+1]) for i in range(n)) + "\n")
PY
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

# Two tests per variant.
printf '#CHROM\tPOS\tID\tA1\tTEST\tP\n' > tmp_assoc.txt
for v in "1 ADD 1e-8" "1 DOM 0.5" "2 ADD 1e-5" "2 DOM 1e-9" "3 ADD 0.003" "3 DOM 0.9" "4 ADD 0.03" "4 DOM 0.04" "5 ADD 1e-6" "5 DOM 0.2" "6 ADD 0.0005" "6 DOM 1e-7" "7 ADD 0.008" "7 DOM 0.6" "8 ADD 0.2" "8 DOM 0.002"; do
    set -- $v
    printf "1\t$((10000 * $1))\tv$1\tG\t$2\t$3\n" >> tmp_assoc.txt
done

# Default: only ADD rows.  v1 is the index; v2 (1e-5), v3 (0.003), v4 (0.03)
# join it.  v6 (5e-4) is above the 1e-4 index threshold.
$plink2 --pfile tmp_data --clump tmp_assoc.txt --out tmp_add > /dev/null
diff tmp_add.clumps - <<'EOF'
#CHROM	POS	ID	P	TOTAL	NONSIG	S0.05	S0.01	S0.001	S0.0001	SP2
1	10000	v1	1e-08	3	0	1	1	0	1	v2,v3
1	50000	v5	1e-06	0	0	0	0	0	0	.
EOF

# --clump-test DOM.
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-test DOM --out tmp_dom > /dev/null
diff tmp_dom.clumps - <<'EOF'
#CHROM	POS	ID	P	TOTAL	NONSIG	S0.05	S0.01	S0.001	S0.0001	SP2
1	20000	v2	1e-09	3	2	1	0	0	0	.
1	60000	v6	1e-07	2	1	0	1	0	0	v8
EOF
# Both tests pooled: index v2/DOM picks up every other v1-v4 row, including
# its own ADD row.  (SP2 omits rows for the index variant itself, so it's
# left out of this comparison.)
pooled='#CHROM	POS	ID	P	TOTAL	NONSIG	S0.05	S0.01	S0.001	S0.0001
1	20000	v2	1e-09	7	2	2	1	0	2
1	60000	v6	1e-07	5	2	0	2	1	0
1	50000	v5	1e-06	1	1	0	0	0	0'
$plink2 --pfile tmp_data --clump cols=-sp2 tmp_assoc.txt --clump-test ADD DOM --out tmp_both > /dev/null
test "$(cat tmp_both.clumps)" = "$pooled"
# An empty --clump-test-field disables test filtering, which also pools.
$plink2 --pfile tmp_data --clump cols=-sp2 tmp_assoc.txt --clump-test-field --out tmp_notest > /dev/null
test "$(cat tmp_notest.clumps)" = "$pooled"

# Renamed columns.  (MARKER/PVAL are not among the default names.)
awk 'BEGIN {FS = OFS = "\t"} NR == 1 {print "#CHROM\tPOS\tMARKER\tEFFA\tMODEL\tPVAL"; next} {print}' tmp_assoc.txt > tmp_ren.txt
ren_fields="--clump-id-field MARKER --clump-p-field PVAL"
$plink2 --pfile tmp_data --clump tmp_ren.txt $ren_fields --clump-a1-field EFFA --clump-test-field MODEL --out tmp_ren > /dev/null
cmp tmp_add.clumps tmp_ren.clumps
$plink2 --pfile tmp_data --clump tmp_ren.txt $ren_fields --clump-test-field MODEL --clump-test DOM --out tmp_ren_dom > /dev/null
cmp tmp_dom.clumps tmp_ren_dom.clumps
# Without --clump-test-field, MODEL isn't a test column, so all rows pool.
$plink2 --pfile tmp_data --clump cols=-sp2 tmp_ren.txt $ren_fields --out tmp_ren_pool > /dev/null
test "$(cat tmp_ren_pool.clumps)" = "$pooled"
# Field-name search order: earlier names win, missing names fall through.
$plink2 --pfile tmp_data --clump tmp_ren.txt --clump-id-field NOPE1 MARKER --clump-p-field PVAL --clump-test-field NOPE2 MODEL --out tmp_ren_order > /dev/null
cmp tmp_add.clumps tmp_ren_order.clumps

# --clump-a1-field, observed via --clump-force-a1 (biallelic A1s are
# otherwise not tracked).  v1's A1 is changed to the REF allele.
awk 'BEGIN {FS = OFS = "\t"} NR > 1 && $3 == "v1" {$4 = "A"} {print}' tmp_ren.txt > tmp_ren_a1.txt
$plink2 --pfile tmp_data --clump cols=+a1 tmp_ren_a1.txt $ren_fields --clump-test-field MODEL --clump-a1-field EFFA --clump-force-a1 --out tmp_a1 > /dev/null
test "$(tail -n +2 tmp_a1.clumps | cut -f 3,4 | tr '\n' ' ')" = "$(printf 'v1\tA v5\tG ')"
test "$(awk -F '\t' '$3 == "v1" {print $NF}' tmp_a1.clumps)" = "v2(G),v3(G)"
if $plink2 --pfile tmp_data --clump tmp_ren_a1.txt $ren_fields --clump-test-field MODEL --clump-force-a1 --out tmp_a1_missing > /dev/null 2>&1; then
    exit 1
fi

# --clump-log10.
awk 'BEGIN {FS = OFS = "\t"} NR == 1 {print "#CHROM\tPOS\tID\tA1\tTEST\tLOG10_P"; next} {$6 = -log($6) / log(10); print}' tmp_assoc.txt > tmp_log10.txt
$plink2 --pfile tmp_data --clump tmp_log10.txt --clump-log10 input-only --out tmp_log_in > /dev/null
cmp tmp_add.clumps tmp_log_in.clumps
log_out='#CHROM	POS	ID	NEG_LOG10_P	TOTAL	NONSIG	S0.05	S0.01	S0.001	S0.0001	SP2
1	10000	v1	8	3	0	1	1	0	1	v2,v3
1	50000	v5	6	0	0	0	0	0	0	.'
$plink2 --pfile tmp_data --clump tmp_log10.txt --clump-log10 --out tmp_log_both > /dev/null
test "$(cat tmp_log_both.clumps)" = "$log_out"
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-log10 output-only --out tmp_log_out > /dev/null
test "$(cat tmp_log_out.clumps)" = "$log_out"
# Read as raw p-values, the -log10 values are out of range.
if $plink2 --pfile tmp_data --clump tmp_log10.txt --clump-p-field LOG10_P --out tmp_log_misread > /dev/null 2>&1; then
    exit 1
fi
grep -q 'p-value > 1' tmp_log_misread.log

# --clump-log10-p1/-p2 are equivalent to --clump-p1/-p2.
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-log10-p1 4 --clump-log10-p2 2 --out tmp_lp > /dev/null
cmp tmp_add.clumps tmp_lp.clumps
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-log10-p1 7 --out tmp_lp7 > /dev/null
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-p1 1e-7 --out tmp_p7 > /dev/null
cmp tmp_lp7.clumps tmp_p7.clumps
test "$(tail -n +2 tmp_lp7.clumps | cut -f 3)" = "v1"
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-log10-p2 1.5 --out tmp_lp2 > /dev/null
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-p2 0.0316227766 --out tmp_p2 > /dev/null
cmp tmp_lp2.clumps tmp_p2.clumps
# v4 (0.03) is now within --clump-p2, so SP2 gains it.
test "$(awk -F '\t' '$3 == "v1" {print $NF}' tmp_lp2.clumps)" = "v2,v3,v4"
for flags in "--clump-p1 1e-7 --clump-log10-p1 7" "--clump-p2 0.01 --clump-log10-p2 2"; do
    if $plink2 --pfile tmp_data --clump tmp_assoc.txt $flags --out tmp_bad > /dev/null 2>&1; then
        exit 1
    fi
done

# --clump-bins: custom boundaries; commas and spaces both separate values.
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-bins 0.001 0.01 --out tmp_bins > /dev/null
diff tmp_bins.clumps - <<'EOF'
#CHROM	POS	ID	P	TOTAL	NONSIG	S0.01	S0.001	SP2
1	10000	v1	1e-08	3	1	1	1	v2,v3
1	50000	v5	1e-06	0	0	0	0	.
EOF
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-bins 0.001,0.01 --out tmp_bins2 > /dev/null
cmp tmp_bins.clumps tmp_bins2.clumps
$plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-bins 0.5 --out tmp_bins3 > /dev/null
test "$(head -n 2 tmp_bins3.clumps | tail -n 1 | cut -f 5-7)" = "$(printf '3\t0\t3')"
# Unsorted or repeated boundaries, and boundaries >= 1, are rejected.
for bins in "0.01 0.001" "0.001,0.001" "0.01 1"; do
    if $plink2 --pfile tmp_data --clump tmp_assoc.txt --clump-bins $bins --out tmp_bad_bins > /dev/null 2>&1; then
        exit 1
    fi
done
if $plink2 --pfile tmp_data --clump cols=-bins tmp_assoc.txt --clump-bins 0.01 --out tmp_bad_bins > /dev/null 2>&1; then
    exit 1
fi

# --clump-unphased.  w1/w2 haplotypes are 0-0, 1-1, 0-1, 1-0 in equal
# numbers (phased r^2 = 0), but the genotype dosages are perfectly
# anticorrelated (unphased r^2 = 1).
{
    printf '##fileformat=VCFv4.2\n##contig=<ID=2,length=1000000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    for i in $(seq 1 60); do printf "\tt$i"; done
    printf '\n2\t1000\tw1\tA\tG\t.\t.\t.\tGT'
    for i in $(seq 1 60); do
        if [ $i -le 30 ]; then printf '\t0|1'; elif [ $i -le 45 ]; then printf '\t0|0'; else printf '\t1|1'; fi
    done
    printf '\n2\t2000\tw2\tA\tG\t.\t.\t.\tGT'
    for i in $(seq 1 60); do
        if [ $i -le 30 ]; then printf '\t0|1'; elif [ $i -le 45 ]; then printf '\t1|1'; else printf '\t0|0'; fi
    done
    printf '\n'
} > tmp_phase.vcf
$plink2 --vcf tmp_phase.vcf --make-pgen --out tmp_phase > /dev/null
printf '#CHROM\tPOS\tID\tP\n2\t1000\tw1\t1e-6\n2\t2000\tw2\t1e-3\n' > tmp_phase_assoc.txt
$plink2 --pfile tmp_phase --clump tmp_phase_assoc.txt --out tmp_phased > /dev/null
test "$(tail -n +2 tmp_phased.clumps | cut -f 3,5,11)" = "$(printf 'w1\t0\t.')"
$plink2 --pfile tmp_phase --clump tmp_phase_assoc.txt --clump-unphased --out tmp_unphased > /dev/null
test "$(tail -n +2 tmp_unphased.clumps | cut -f 3,5,11)" = "$(printf 'w1\t1\tw2')"
