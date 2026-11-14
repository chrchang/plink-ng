#!/bin/bash

# Import and output flags with no other test coverage:
# --dosage-erase-threshold, --lax-bgen-import, --lax-chrx-import,
# --keep-autoconv, --rename-chrs, --merge-par, --alt1-allele, --no-id-header,
# --make-just-fam, --var-id-multi, --var-id-multi-nonsnp, and --zst-level.
#
# Regression checks for --rename-chrs:
# - Renaming onto a chromosome code that another (not renamed) chromosome
#   already had left two chromosome blocks sharing one code, and the sorted
#   .pvar write then hit an assertion (or, in release builds, wrote
#   inconsistent output).  This is now an error.
# - When the input had split chromosomes, the rename was silently ignored.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# --- --lax-chrx-import, --keep-autoconv, --no-id-header, --make-just-fam ---
cat > tmp_x.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##contig=<ID=X,length=200000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
1	100	v1	A	G	.	.	.	GT	0/1	0/0	1/1
X	3000000	x1	G	A	.	.	.	GT	1	0/1	0/0
VCF
# chrX without sex information is an error unless --lax-chrx-import.
if $plink2 --vcf tmp_x.vcf --make-pgen --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
$plink2 --vcf tmp_x.vcf --lax-chrx-import --make-pgen --out tmp_x > /dev/null
test "$(grep -v '^#' tmp_x.pvar | cut -f 1 | tr '\n' ' ')" = "1 X "

# The autoconverted fileset is normally deleted at the end of the run.
rm -f tmp_ka.*
$plink2 --vcf tmp_x.vcf --lax-chrx-import --freq --out tmp_ka > /dev/null
test ! -e tmp_ka.pgen
$plink2 --vcf tmp_x.vcf --lax-chrx-import --keep-autoconv --freq --out tmp_ka > /dev/null
diff -q tmp_ka.pgen tmp_x.pgen
# (The autoconverted .pvar keeps the QUAL/FILTER columns.)
diff -q <(grep -v '^##' tmp_ka.pvar | cut -f 1-5) <(grep -v '^##' tmp_x.pvar)
diff -q tmp_ka.psam tmp_x.psam
rm -f tmp_kz.*
$plink2 --vcf tmp_x.vcf --lax-chrx-import --keep-autoconv vzs --freq --out tmp_kz > /dev/null
test ! -e tmp_kz.pvar
diff -q <(zstd -dc tmp_kz.pvar.zst | grep -v '^##' | cut -f 1-5) <(grep -v '^##' tmp_x.pvar)
if $plink2 --vcf tmp_x.vcf --lax-chrx-import --keep-autoconv --make-pgen --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# .id output: '#IID' header by default; --no-id-header drops it and forces
# FID/IID columns, unless 'iid-only'.
$plink2 --pfile tmp_x --write-samples --out tmp_ws > /dev/null
diff tmp_ws.id <(printf '#IID\nA\nB\nC\n')
$plink2 --pfile tmp_x --write-samples --no-id-header --out tmp_ws2 > /dev/null
diff tmp_ws2.id <(printf '0\tA\n0\tB\n0\tC\n')
$plink2 --pfile tmp_x --write-samples --no-id-header iid-only --out tmp_ws3 > /dev/null
diff tmp_ws3.id <(printf 'A\nB\nC\n')

# --make-just-fam writes just the .fam that --make-bed would.
$plink2 --pfile tmp_x --make-just-fam --out tmp_mjf > /dev/null
$plink2 --pfile tmp_x --make-bed --out tmp_mb > /dev/null
diff -q tmp_mjf.fam tmp_mb.fam
test ! -e tmp_mjf.bed
test ! -e tmp_mjf.bim

# --- --zst-level ---
$plink2 --dummy 20 3000 --seed 3 --make-pgen --out tmp_dummy > /dev/null
$plink2 --pfile tmp_dummy --zst-level 1 --make-pgen vzs --out tmp_z1 > /dev/null
$plink2 --pfile tmp_dummy --zst-level 19 --make-pgen vzs --out tmp_z19 > /dev/null
diff -q <(zstd -dc tmp_z1.pvar.zst) <(zstd -dc tmp_z19.pvar.zst)
test "$(wc -c < tmp_z19.pvar.zst)" -lt "$(wc -c < tmp_z1.pvar.zst)"
for lvl in 0 23; do
    if $plink2 --pfile tmp_dummy --zst-level $lvl --make-pgen vzs --out tmp_bad > /dev/null 2>&1; then
        exit 1
    fi
done

# --- --dosage-erase-threshold ---
# Dosages 1.05 and 0.08 are within 0.1 of a hardcall; 1.85 and 0.5 aren't.
cat > tmp_ds.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C	D
1	100	v1	A	G	.	.	.	GT:DS	0/1:1.05	0/0:0.08	1/1:1.85	0/1:0.5
VCF
# $1: output prefix; remaining args: import flags.  Prints the re-exported
# sample fields.
dosage_export() {
    local out=$1
    shift
    $plink2 --vcf tmp_ds.vcf dosage=DS "$@" --make-pgen --out $out > /dev/null
    $plink2 --pfile $out --export vcf vcf-dosage=DS --out $out > /dev/null
    grep -v '^#' $out.vcf | cut -f 10- | tr '\t' ' '
}
test "$(dosage_export tmp_de0)" = "0/1:1.05 0/0:0.08 ./.:1.85 ./.:0.5"
test "$(dosage_export tmp_de1 --dosage-erase-threshold 0.1)" = "0/1 0/0 ./.:1.85 ./.:0.5"
test "$(dosage_export tmp_de2 --hard-call-threshold 0.2 --dosage-erase-threshold 0.1)" = "0/1 0/0 1/1:1.85 ./.:0.5"
# It can't exceed the hardcall threshold (default 0.1).
if $plink2 --vcf tmp_ds.vcf dosage=DS --dosage-erase-threshold 0.2 --make-pgen --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi

# --- --lax-bgen-import ---
# Overstate the variant count in a .bgen header (offset 8).
$plink2 --pfile tmp_dummy --export bgen-1.2 --out tmp_bgen > /dev/null
python3 -c "
import struct
b = bytearray(open('tmp_bgen.bgen', 'rb').read())
m = struct.unpack_from('<I', b, 8)[0]
struct.pack_into('<I', b, 8, m + 5)
open('tmp_bgen_over.bgen', 'wb').write(b)
"
$plink2 --bgen tmp_bgen.bgen ref-first --sample tmp_bgen.sample --make-pgen --out tmp_b0 > /dev/null
if $plink2 --bgen tmp_bgen_over.bgen ref-first --sample tmp_bgen.sample --make-pgen --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
$plink2 --bgen tmp_bgen_over.bgen ref-first --sample tmp_bgen.sample --lax-bgen-import --make-pgen --out tmp_b1 > /dev/null
diff -q tmp_b0.pgen tmp_b1.pgen
diff -q tmp_b0.pvar tmp_b1.pvar

# --- --alt1-allele, --var-id-multi[-nonsnp], --rename-chrs ---
cat > tmp_m.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=10000>
##contig=<ID=scafB,length=10000>
##contig=<ID=scafA,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
1	100	s1	A	G	.	.	.	GT	0/1	0/0	1/1
1	200	m1	C	T,G	.	.	.	GT	1/2	2/2	0/1
1	300	i1	AT	A,ATT	.	.	.	GT	0/1	0/2	0/0
1	400	mx	A	C,AG	.	.	.	GT	0/1	0/2	1/1
scafB	50	b1	G	A	.	.	.	GT	0/1	0/0	0/0
scafA	70	a1	T	C	.	.	.	GT	0/0	0/1	0/0
VCF
$plink2 --vcf tmp_m.vcf --make-pgen --out tmp_m > /dev/null
# $1: dataset prefix.  Prints ID/REF/ALT and genotypes from a VCF re-export.
vcf_rows() {
    $plink2 --pfile $1 --export vcf --out $1_v > /dev/null
    grep -v '^#' $1_v.vcf | cut -f 3-5,10- | tr '\t' ' '
}
# --alt1-allele moves the named allele to ALT1, and genotypes follow.
printf 'm1\tG\ni1\tATT\nmx\tAG\n' > tmp_a1.txt
$plink2 --pfile tmp_m --alt1-allele tmp_a1.txt 2 1 --make-pgen --out tmp_ma > /dev/null
diff <(vcf_rows tmp_ma) - <<'EOF'
s1 A G 0/1 0/0 1/1
m1 C G,T 1/2 1/1 0/2
i1 AT ATT,A 0/2 0/1 0/0
mx A AG,C 0/2 0/1 2/2
b1 G A 0/1 0/0 0/0
a1 T C 0/0 0/1 0/0
EOF
# Moving a known REF allele to ALT1 requires 'force'.
printf 's1\tA\n' > tmp_a1ref.txt
if $plink2 --pfile tmp_m --alt1-allele tmp_a1ref.txt 2 1 --make-pgen --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
$plink2 --pfile tmp_m --alt1-allele force tmp_a1ref.txt 2 1 --make-pgen --out tmp_maf > /dev/null
test "$(vcf_rows tmp_maf | head -n 1)" = "s1 G A 0/1 1/1 0/0"

# Multiallelic variants get --var-id-multi, or --var-id-multi-nonsnp when any
# allele isn't a single base.
# $1: expected IDs; remaining args: extra flags.
check_ids() {
    local expected=$1
    shift
    $plink2 --pfile tmp_m --set-all-var-ids '@:#:$r:$a' "$@" --make-pgen --out tmp_vi > /dev/null
    test "$(grep -v '^#' tmp_vi.pvar | cut -f 3 | tr '\n' ' ')" = "$expected "
}
check_ids "1:100:A:G 1:200:C:T 1:300:AT:A 1:400:A:C scafB:50:G:A scafA:70:T:C"
check_ids "1:100:A:G 1:200:multi 1:300:multi 1:400:multi scafB:50:G:A scafA:70:T:C" --var-id-multi @:#:multi
check_ids "1:100:A:G 1:200:multi 1:300:nsmulti 1:400:nsmulti scafB:50:G:A scafA:70:T:C" --var-id-multi @:#:multi --var-id-multi-nonsnp @:#:nsmulti
check_ids "1:100:A:G 1:200:C:T 1:300:nsmulti 1:400:nsmulti scafB:50:G:A scafA:70:T:C" --var-id-multi-nonsnp @:#:nsmulti

# --rename-chrs requires --sort-vars, which then sorts by the new names.
# scafB precedes scafA in the input; after renaming, aa_scaf (was scafA)
# sorts first.
printf 'scafA\taa_scaf\nscafB\tzz_scaf\n' > tmp_rn.txt
if $plink2 --pfile tmp_m --rename-chrs tmp_rn.txt --make-pgen --out tmp_bad > /dev/null 2>&1; then
    exit 1
fi
$plink2 --pfile tmp_m --rename-chrs tmp_rn.txt --sort-vars --make-pgen --out tmp_rn > /dev/null
test "$(grep -v '^#' tmp_rn.pvar | cut -f 1,3 | tr '\t\n' ': ')" = "1:s1 1:m1 1:i1 1:mx aa_scaf:a1 zz_scaf:b1 "
# Genotypes move with their variants.
diff <(vcf_rows tmp_rn | sort) <(vcf_rows tmp_m | sort)
# Swapping two names is fine.
printf 'scafA\tscafB\nscafB\tscafA\n' > tmp_rn_swap.txt
$plink2 --pfile tmp_m --rename-chrs tmp_rn_swap.txt --sort-vars --make-pgen --out tmp_rn_swap > /dev/null
test "$(grep -v '^#' tmp_rn_swap.pvar | cut -f 1,3 | tail -n 2 | tr '\t\n' ': ')" = "scafA:b1 scafB:a1 "
# Merging two chromosomes into one name is an error, not a crash.
printf 'scafA\tscafB\n' > tmp_rn_merge.txt
printf 'scafA\t1\n' > tmp_rn_merge2.txt
for f in tmp_rn_merge.txt tmp_rn_merge2.txt; do
    if $plink2 --pfile tmp_m --rename-chrs $f --sort-vars --make-pgen --out tmp_rn_bad > tmp_rn_bad.out 2>&1; then
        exit 1
    fi
    grep -q 'would merge two chromosomes' tmp_rn_bad.log
done
# Input with split chromosomes (only allowed with --sort-vars).
cat > tmp_split.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=scafA,length=10000>
##contig=<ID=1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
scafA	70	a1	T	C	.	.	.	GT	0/0	0/1	0/0
1	100	s1	A	G	.	.	.	GT	0/1	0/0	1/1
scafA	90	a2	T	C	.	.	.	GT	0/1	0/1	1/1
1	50	s0	A	G	.	.	.	GT	1/1	0/0	0/1
VCF
printf 'scafA\tnewname\n' > tmp_rn_split.txt
$plink2 --vcf tmp_split.vcf --rename-chrs tmp_rn_split.txt --sort-vars --make-pgen --out tmp_rn_split > /dev/null
test "$(grep -v '^#' tmp_rn_split.pvar | cut -f 1,3 | tr '\t\n' ': ')" = "1:s0 1:s1 newname:a1 newname:a2 "
diff <(vcf_rows tmp_rn_split | sort) <(grep -v '^#' tmp_split.vcf | cut -f 3-5,10- | tr '\t' ' ' | sort)

# --- --merge-par ---
# --split-par b38 moves the first two and the last variant to PAR1/PAR2;
# --merge-par moves them back, and a second split recovers the first one
# exactly.
cat > tmp_par.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=X,length=156040895>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M1	F1	M2
X	60001	p1a	A	G	.	.	.	GT	0/1	0/0	1/1
X	2000000	p1b	C	T	.	.	.	GT	0/0	0/1	0/1
X	3000000	x1	G	A	.	.	.	GT	1	0/1	0
X	100000000	x2	T	C	.	.	.	GT	0	1/1	1
X	155800000	p2a	A	C	.	.	.	GT	0/1	0/1	0/0
VCF
printf '#IID\tSEX\nM1\t1\nF1\t2\nM2\t1\n' > tmp_par.psam
$plink2 --vcf tmp_par.vcf --psam tmp_par.psam --split-par b38 --make-pgen --out tmp_ps > /dev/null
test "$(grep -v '^#' tmp_ps.pvar | cut -f 1 | tr '\n' ' ')" = "PAR1 PAR1 X X PAR2 "
$plink2 --pfile tmp_ps --merge-par --make-pgen --out tmp_pm > /dev/null
test "$(grep -v '^#' tmp_pm.pvar | cut -f 1-3 | tr '\t\n' ': ')" = "X:60001:p1a X:2000000:p1b X:3000000:x1 X:100000000:x2 X:155800000:p2a "
$plink2 --pfile tmp_pm --split-par b38 --make-pgen --out tmp_ps2 > /dev/null
diff -q tmp_ps.pgen tmp_ps2.pgen
diff -q <(grep -v '^##' tmp_ps.pvar) <(grep -v '^##' tmp_ps2.pvar)
