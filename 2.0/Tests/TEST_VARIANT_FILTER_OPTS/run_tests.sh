#!/bin/bash

# Variant filter flags with no other test coverage: --extract-intersect,
# --force-intersect, --bed-border-bp, --bed-border-kb, --mach-r2-filter,
# --minimac3-r2-filter, and --maf-succ.
#
# Regression check: --mach-r2-filter with no arguments (documented defaults
# 0.1 and 2.0) set only the minimum, and since a zero maximum is what marks
# the filter as inactive, it silently did nothing.
# --bed-border-bp/-kb without a BED-file filter used to be silently ignored;
# it's now an error.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# v1..v8 are 1 kb apart on chr1; w1/w2 are on chr2.  S1 carries the only ALT
# allele at v1..v8, and w2 is monomorphic.
{
    printf '##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##contig=<ID=2,length=100000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n'
    for i in 1 2 3 4 5 6 7 8; do
        printf '1\t%d000\tv%d\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\t0/0\t0/0\n' $i $i
    done
    printf '2\t1000\tw1\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/0\t0/0\n'
    printf '2\t5000\tw2\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0\n'
} > tmp_data.vcf
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

# $1: dataset prefix; $2: expected space-separated variant IDs; remaining
# args: plink2 flags.
check_variants() {
    local prefix=$1
    local expected=$2
    shift 2
    rm -f tmp_v.snplist
    $plink2 --pfile $prefix "$@" --write-snplist --out tmp_v > /dev/null
    test "$(tr '\n' ' ' < tmp_v.snplist)" = "$expected "
}

# --- --extract-intersect / --force-intersect ---
printf 'v1\nv2\nv3\nv4\nw1\n' > tmp_a.txt
printf 'v3\nv4\nv5\nw1\n' > tmp_b.txt
printf 'v4\nw1\nw2\n' > tmp_c.txt
check_variants tmp_data "v1 v2 v3 v4 v5 w1" --extract tmp_a.txt tmp_b.txt
check_variants tmp_data "v3 v4 w1" --extract-intersect tmp_a.txt tmp_b.txt
check_variants tmp_data "v4 w1" --extract-intersect tmp_a.txt tmp_b.txt tmp_c.txt
check_variants tmp_data "v1 v2 v3 v4 w1" --extract-intersect tmp_a.txt
# Multiple inclusion filters need --force-intersect, which then intersects.
for flags in "--extract tmp_a.txt --extract-intersect tmp_b.txt" "--extract tmp_a.txt --snps v2-v6"; do
    if $plink2 --pfile tmp_data $flags --write-snplist --out tmp_bad > /dev/null 2>&1; then
        exit 1
    fi
done
check_variants tmp_data "v3 v4 w1" --extract tmp_a.txt --extract-intersect tmp_b.txt --force-intersect
check_variants tmp_data "v2 v3 v4" --extract tmp_a.txt --snps v2-v6 --force-intersect
check_variants tmp_data "v3 v4" --extract-intersect tmp_a.txt tmp_b.txt --from v2 --to v6 --force-intersect

# --- BED intervals, --bed-border-bp/-kb ---
# bed0 is 0-based half-open, bed1 1-based closed.
# tmp_r0: chr1 2501..5500 (bed0) or 2500..5500 (bed1), chr2 1..1000.
# tmp_r1: chr1 3001..6000 (bed0) or 3000..6000 (bed1), chr2 901..1000.
printf '1\t2500\t5500\n2\t0\t1000\n' > tmp_r0.bed
printf '1\t3000\t6000\n2\t900\t1000\n' > tmp_r1.bed
check_variants tmp_data "v3 v4 v5 w1" --extract bed0 tmp_r0.bed
check_variants tmp_data "v4 v5 w1" --extract-intersect bed0 tmp_r0.bed tmp_r1.bed
check_variants tmp_data "v3 v4 v5 w1" --extract-intersect bed1 tmp_r0.bed tmp_r1.bed
# Borders stretch both ends: bed0 [2000, 6000) includes 6000 but not 2000.
check_variants tmp_data "v3 v4 v5 v6 w1" --extract bed0 tmp_r0.bed --bed-border-bp 500
check_variants tmp_data "v3 v4 v5 w1" --extract bed0 tmp_r0.bed --bed-border-bp 499
check_variants tmp_data "v2 v3 v4 v5 v6 w1" --extract bed1 tmp_r0.bed --bed-border-kb 1
check_variants tmp_data "v1 v7 v8 w2" --exclude bed1 tmp_r0.bed --bed-border-bp 1000
check_variants tmp_data "v3 v4 v5 v6 w1" --extract-intersect bed0 tmp_r0.bed tmp_r1.bed --bed-border-bp 1000
# A border with no BED filter (including an ID-list --extract) is an error.
for flags in "--bed-border-bp 500" "--bed-border-kb 0" "--extract tmp_a.txt --bed-border-bp 500"; do
    if $plink2 --pfile tmp_data $flags --write-snplist --out tmp_bad > /dev/null 2>&1; then
        exit 1
    fi
    grep -q 'must be used with a BED-file' tmp_bad.log
done

# --- --maf-succ ---
# S1's lone ALT call gives v1..v8 a MAF of 1/8; w2 has MAF 0, which becomes
# (0+1)/(8+2) = 0.1 with a pseudocount of 1.
check_variants tmp_data "v1 v2 v3 v4 v5 v6 v7 v8 w1" --maf 0.05
check_variants tmp_data "v1 v2 v3 v4 v5 v6 v7 v8 w1 w2" --maf 0.05 --maf-succ
check_variants tmp_data "v1 v2 v3 v4 v5 v6 v7 v8 w1 w2" --maf 0.05 --af-pseudocount 1
# v1..v8: (1+1)/(8+2) = 0.2 passes --max-maf 0.2 only with the pseudocount.
check_variants tmp_data "w2" --max-maf 0.15 --maf-succ
if $plink2 --pfile tmp_data --maf-succ --af-pseudocount 1 --maf 0.05 --write-snplist --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# --- --mach-r2-filter / --minimac3-r2-filter ---
# Phased dosages with increasing uncertainty, plus a monomorphic variant
# (r2 = nan, never excluded).
python3 - <<'EOF'
import random
random.seed(4)
n = 40
with open("tmp_r2.vcf", "w") as f:
    f.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haploid dosages">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join("S%d" % i for i in range(n)) + "\n")
    for v in range(10):
        p = random.uniform(0.1, 0.5)
        noise = v / 9.0
        fields = []
        for s in range(n):
            h = []
            for _ in range(2):
                a = 1 if random.random() < p else 0
                d = a * (1 - noise) + noise * p + random.uniform(-0.05, 0.05) * noise
                h.append(min(1, max(0, d)))
            fields.append("%d|%d:%.3f,%.3f" % (round(h[0]), round(h[1]), h[0], h[1]))
        f.write("1\t%d\tr%d\tA\tG\t.\t.\t.\tGT:HDS\t" % (1000 * (v + 1), v) + "\t".join(fields) + "\n")
    f.write("1\t20000\tmono\tA\tG\t.\t.\t.\tGT:HDS\t" + "\t".join(["0|0:0,0"] * n) + "\n")
EOF
$plink2 --vcf tmp_r2.vcf dosage=HDS --make-pgen --out tmp_r2 > /dev/null
$plink2 --pfile tmp_r2 --freq cols=+machr2 --out tmp_mach > /dev/null
$plink2 --pfile tmp_r2 --freq cols=+minimac3r2 --out tmp_mm3 > /dev/null
# $1: .afreq file; $2/$3: min/max.  Prints the IDs a filter should keep.
r2_expected() {
    awk -F '\t' -v lo="$2" -v hi="$3" 'NR > 1 {if ($6 == "nan" || ($6 + 0 >= lo && $6 + 0 <= hi)) {printf "%s%s", sep, $2; sep = " "}}' "$1"
}
# The fixture spans both sides of each threshold.
test "$(r2_expected tmp_mach.afreq 0.1 2.0)" != "$(r2_expected tmp_mach.afreq 0 1000)"
test "$(r2_expected tmp_mach.afreq 0.5 0.9)" != "$(r2_expected tmp_mach.afreq 0.5 2.0)"
check_variants tmp_r2 "$(r2_expected tmp_mach.afreq 0.1 2.0)" --mach-r2-filter
check_variants tmp_r2 "$(r2_expected tmp_mach.afreq 0.5 2.0)" --mach-r2-filter 0.5
check_variants tmp_r2 "$(r2_expected tmp_mach.afreq 0.5 0.9)" --mach-r2-filter 0.5 0.9
check_variants tmp_r2 "$(r2_expected tmp_mm3.afreq 0.5 1000)" --minimac3-r2-filter 0.5
check_variants tmp_r2 "$(r2_expected tmp_mm3.afreq 0.3 0.8)" --minimac3-r2-filter 0.3 0.8
for flags in "--mach-r2-filter --minimac3-r2-filter 0.5" "--mach-r2-filter 0.9 0.5" "--mach-r2-filter -1"; do
    if $plink2 --pfile tmp_r2 $flags --write-snplist --out tmp_bad > /dev/null 2>&1; then
        exit 1
    fi
done
