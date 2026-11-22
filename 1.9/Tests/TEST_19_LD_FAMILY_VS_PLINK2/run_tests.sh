#!/bin/bash

# PLINK 1.9 relationship-matrix, LD, clumping and family (Mendel error)
# commands and options, checked against the equivalent plink2 flags.
#
# Usage: ./run_tests.sh {plink 1.9 dir} {plink2 build dir} {up to 2 plink2 args}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  {plink2 build dir} defaults to
# ../../../2.0/build_dynamic.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}
plink2="${2:-../../../2.0/build_dynamic}/plink2 ${3:-} ${4:-}"
compare_cols="../compare_cols.awk"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
spaces() {
    tr '\t' ' ' < "$1"
}

# 203 samples, 150 variants in LD blocks (see make_ld.awk): 100 on chr1 and 50
# on chr2, 2 kb and 0.05 cM apart.
awk -v n=203 -v m=150 -v seed=17 -f make_ld.awk > tmp_ld.ped
awk 'BEGIN { for (j = 1; j <= 150; ++j) print (j <= 100)? 1 : 2, "v" j, j * 0.05, 2000 * j }' > tmp_ld.map
plink --file tmp_ld --keep-allele-order --make-bed --out tmp_ld > /dev/null 2>&1

# --- Relationship matrix ---

# Prints every value of a matrix file, one per line.
values() {
    awk '{ for (i = 1; i <= NF; ++i) print $i }' "$1"
}
rel_ok=0
for shape in "" square; do
    plink --bfile tmp_ld --make-rel $shape --out tmp_rel19 > /dev/null 2>&1
    $plink2 --bfile tmp_ld --make-rel $shape --out tmp_rel2 > /dev/null 2>&1
    awk -v c1=1 -v c2=1 -v skip1=0 -v skip2=0 -f $compare_cols <(values tmp_rel19.rel) <(values tmp_rel2.rel) || rel_ok=1
done
cmp -s tmp_rel19.rel.id <(grep -v '^#' tmp_rel2.rel.id) || rel_ok=1
check "--make-rel (triangle, square) matches plink2" $rel_ok
# The binary matrices may differ in the last bits.
plink --bfile tmp_ld --make-rel bin --out tmp_relb19 > /dev/null 2>&1
$plink2 --bfile tmp_ld --make-rel bin --out tmp_relb2 > /dev/null 2>&1
binvalues() {
    od -A n -v -t f8 "$1" | tr -s ' ' '\n' | grep -v '^$'
}
awk -v c1=1 -v c2=1 -v skip1=0 -v skip2=0 -v tol=1e-9 -v abstol=1e-12 -f $compare_cols <(binvalues tmp_relb19.rel.bin) <(binvalues tmp_relb2.rel.bin) && \
    test "$(binvalues tmp_relb19.rel.bin | wc -l)" -eq $((203 * 203))
check "--make-rel bin matches plink2" $?

# --- LD ---

plink --bfile tmp_ld --indep 20 5 2 --out tmp_in19 > /dev/null 2>&1
$plink2 --bfile tmp_ld --indep 20 5 2 --out tmp_in2 > /dev/null 2>&1
cmp -s tmp_in19.prune.in tmp_in2.prune.in && cmp -s tmp_in19.prune.out tmp_in2.prune.out && test -s tmp_in19.prune.out
check "--indep matches plink2" $?

# --ld-snp/--ld-snps: PLINK 1.9 also reports each listed variant with itself,
# and lists every listed variant's whole window, while plink2 reports each
# pair once; so unordered pairs other than self-pairs are compared.
ld_pairs() {
    awk 'NR > 1 && $3 != $6 {
        a = $3
        b = $6
        if (substr(a, 2) + 0 > substr(b, 2) + 0) {
            t = a
            a = b
            b = t
        }
        print a, b, $7
    }' "$1" | LC_ALL=C sort -u
}
# $1: label, $2: PLINK 1.9 flags, $3: plink2 flags.
same_ld() {
    rm -f tmp_r19.ld tmp_r2.vcor
    plink --bfile tmp_ld $2 --out tmp_r19 > /dev/null 2>&1
    $plink2 --bfile tmp_ld $3 --out tmp_r2 > /dev/null 2>&1
    test -s tmp_r19.ld && awk -v c1=1,2,3 -v c2=1,2,3 -v skip1=0 -v skip2=0 -f $compare_cols <(ld_pairs tmp_r19.ld) <(ld_pairs tmp_r2.vcor)
    check "$1" $?
}
same_ld "--ld-snp matches plink2" "--r2 --ld-snp v15 --ld-window-r2 0" "--r2-unphased --ld-snp v15 --ld-window 10 --ld-window-r2 0"
same_ld "--ld-snps matches plink2" "--r2 --ld-snps v10-v20,v105 --ld-window-r2 0" "--r2-unphased --ld-snps v10-v20,v105 --ld-window 10 --ld-window-r2 0"
same_ld "--ld-window-cm matches plink2" "--r2 --ld-window-cm 0.3 --ld-window 99999 --ld-window-kb 100000 --ld-window-r2 0" "--r2-unphased --ld-window-cm 0.3 --ld-window 99999 --ld-window-kb 100000 --ld-window-r2 0"

# --- Clumping ---

# Synthetic p-values: tiny for block positions 2 and 3 (bit 0 and bit 1
# variants, nearly independent), moderate for position 4 (in LD with both),
# smaller for every 7th variant.  PALT is an alternative p-value column;
# tmp_res_x.txt has scrambled SNP IDs, with the real ones in MYID.
awk 'BEGIN { s = 23; print "CHR SNP BP P PALT MYID" } {
    s = (s * 69069 + 1) % 4294967296; u = int(s / 65536) % 10000 / 10000
    k = NR % 10
    p = (k == 2 || k == 3)? u / 1e5 : (k == 4)? u / 1e3 : (NR % 7 == 0)? u / 50 : u
    s = (s * 69069 + 1) % 4294967296; q = int(s / 65536) % 10000 / 1e6
    print $1, $2, $4, p, q, $2
}' tmp_ld.bim > tmp_res.txt
awk 'NR == 1 { print; next } { $2 = "x" $2; print }' tmp_res.txt > tmp_res_x.txt
# Index variant, P, TOTAL, NSIG, S05, S01, S001, S0001, SP2 from each report.
# (PLINK 1.9 appends "(1)" to SP2 entries and writes NONE for an empty one.)
clumps19() {
    awk 'NR > 1 && NF > 5 { sp = $12; gsub(/\([0-9]+\)/, "", sp); if (sp == "NONE") sp = "."; print $3, $5, $6, $7, $8, $9, $10, $11, sp }' "$1"
}
clumps2() {
    awk '!/^#/ { print $3, $4, $5, $6, $7, $8, $9, $10, $11 }' "$1"
}
plink --bfile tmp_ld --clump tmp_res.txt --out tmp_cdef > /dev/null 2>&1
# $1: label, $2: results file, $3: PLINK 1.9 flags, $4: plink2 flags, $5:
# nonempty if the flags must change PLINK 1.9's report from the defaults.
same_clump() {
    rm -f tmp_c19.clumped tmp_c2.clumps
    plink --bfile tmp_ld --clump "$2" $3 --out tmp_c19 > /dev/null 2>&1
    $plink2 --bfile tmp_ld --clump "$2" $4 --out tmp_c2 > /dev/null 2>&1
    test -s tmp_c19.clumped && awk -v c1=1,2,3,4,5,6,7,8,9 -v c2=1,2,3,4,5,6,7,8,9 -v skip1=0 -v skip2=0 -f $compare_cols <(clumps19 tmp_c19.clumped) <(clumps2 tmp_c2.clumps) && \
        { [[ -z "${5:-}" ]] || ! cmp -s tmp_c19.clumped tmp_cdef.clumped; }
    check "$1" $?
}
same_clump "--clump defaults match plink2" tmp_res.txt "" ""
same_clump "--clump-kb matches plink2" tmp_res.txt "--clump-kb 5" "--clump-kb 5" changes
same_clump "--clump-r2 matches plink2" tmp_res.txt "--clump-r2 0.9" "--clump-r2 0.9" changes
same_clump "--clump-p1/--clump-p2 match plink2" tmp_res.txt "--clump-p1 1e-3 --clump-p2 0.1" "--clump-p1 1e-3 --clump-p2 0.1" changes
same_clump "--clump-field matches plink2 --clump-p-field" tmp_res.txt "--clump-field PALT --clump-p1 1e-3" "--clump-p-field PALT --clump-p1 1e-3" changes
same_clump "--clump-snp-field matches plink2 --clump-id-field" tmp_res_x.txt "--clump-snp-field MYID" "--clump-id-field MYID"
plink --bfile tmp_ld --clump tmp_res_x.txt --out tmp_cx > /dev/null 2>&1
test ! -e tmp_cx.clumped
check "--clump without --clump-snp-field finds no scrambled IDs" $?
plink --bfile tmp_ld --clump tmp_res.txt --clump-r2 0.2 --out tmp_cnov > /dev/null 2>&1
same_clump "--clump-allow-overlap matches plink2" tmp_res.txt "--clump-r2 0.2 --clump-allow-overlap" "--clump-r2 0.2 --clump-allow-overlap"
! cmp -s tmp_c19.clumped tmp_cnov.clumped
check "--clump-allow-overlap changes the clumps" $?

# --clump-range/--clump-range-border: PLINK 1.9 writes the ranges and clump
# bounds to .clumped.ranges, and plink2 adds them as columns.  GB and GC
# overlap.
printf '1 1000 21000 GA\n1 50000 75000 GB\n1 70000 90000 GC\n2 210000 240000 GD\n' > tmp_ranges.txt
ranges19() {
    awk 'NR > 1 { split($5, a, ":"); split(a[2], b, "\\.\\."); r = $7; gsub(/[][]/, "", r); if (r == "") r = "."; print $2, $3, b[1], b[2], r }' "$1" | LC_ALL=C sort
}
ranges2() {
    awk '!/^#/ { print $3, $4, $6, $7, $14 }' "$1" | LC_ALL=C sort
}
range_ok=0
for border in 0 3; do
    plink --bfile tmp_ld --clump tmp_res.txt --clump-range tmp_ranges.txt --clump-range-border $border --out tmp_cr19_$border > /dev/null 2>&1
    $plink2 --bfile tmp_ld --clump cols=+bounds tmp_res.txt --clump-range tmp_ranges.txt --clump-range-border $border --out tmp_cr2 > /dev/null 2>&1
    awk -v c1=1,2,3,4,5 -v c2=1,2,3,4,5 -v skip1=0 -v skip2=0 -f $compare_cols <(ranges19 tmp_cr19_$border.clumped.ranges) <(ranges2 tmp_cr2.clumps) || range_ok=1
done
! cmp -s tmp_cr19_0.clumped.ranges tmp_cr19_3.clumped.ranges && grep -q 'GB,GC' tmp_cr19_3.clumped.ranges || range_ok=1
check "--clump-range/--clump-range-border match plink2" $range_ok

# --- Mendel errors and founders ---

# 61 trios over 120 variants (see ../TEST_19_TDT_EPI_CLUSTER/make_trios.awk),
# with every 5th mother removed (leaving duos).
awk -v n_trio=61 -v n_var=120 -v seed=29 -f ../TEST_19_TDT_EPI_CLUSTER/make_trios.awk > tmp_trio.ped
awk 'BEGIN { for (j = 1; j <= 120; ++j) print 1, "v" j, 0, 1000 * j }' > tmp_trio.map
awk '$2 ~ /^mo/ && substr($2, 3) % 5 == 0 { next } { print }' tmp_trio.ped > tmp_duo.ped
cp tmp_trio.map tmp_duo.map
plink --file tmp_duo --keep-allele-order --make-bed --out tmp_duo > /dev/null 2>&1

# The .mendel ERROR column formats differ, so only the first 5 columns are
# compared, with each error code folded into its mirror-image class (1/2,
# 3/6, 4/7, 5/8), since the two versions may orient a variant's alleles
# differently.  The other reports are compared without headers, sorted.
mendel_rows() {
    tail -n +2 "$1" | tr -s ' \t' ' ' | sed 's/^ //; s/ $//' | \
        awk -v err="$2" '{ if (err) { c = $5; $5 = (c == 2)? 1 : (c == 6)? 3 : (c == 7)? 4 : (c == 8)? 5 : c; NF = 5 } print }' | LC_ALL=C sort
}
same_mendel() {
    local ok=0
    for ext in mendel imendel fmendel lmendel; do
        cmp -s <(mendel_rows tmp_me19.$ext $([[ $ext == mendel ]] && echo 1 || echo 0)) <(mendel_rows tmp_me2.$ext $([[ $ext == mendel ]] && echo 1 || echo 0)) || ok=1
    done
    test -s tmp_me19.mendel || ok=1
    check "$1" $ok
}
plink --bfile tmp_duo --mendel --out tmp_me19 > /dev/null 2>&1
$plink2 --bfile tmp_duo --mendel --out tmp_me2 > /dev/null 2>&1
same_mendel "--mendel matches plink2"
cp tmp_me19.lmendel tmp_nodu.lmendel
plink --bfile tmp_duo --mendel --mendel-duos --out tmp_me19 > /dev/null 2>&1
$plink2 --bfile tmp_duo --mendel --mendel-duos --out tmp_me2 > /dev/null 2>&1
same_mendel "--mendel-duos matches plink2"
! cmp -s tmp_nodu.lmendel tmp_me19.lmendel
check "--mendel-duos changes the error counts" $?

# --me, --me-exclude-one and --set-me-missing.  One father's and one mother's
# genotypes are partly flipped to homozygous, so that their families have
# many errors.
awk 'BEGIN { s = 3 } {
    if ($2 == "fa7" || $2 == "mo12") {
        for (i = 7; i <= NF; i += 2) {
            s = (s * 69069 + 1) % 4294967296
            if (int(s / 65536) % 3 == 0) {
                $i = ($i == "A")? "C" : "A"
                $(i + 1) = $i
            }
        }
    }
    print
}' tmp_trio.ped > tmp_me.ped
cp tmp_trio.map tmp_me.map
plink --file tmp_me --keep-allele-order --make-bed --out tmp_me > /dev/null 2>&1
# $1: label, $2: input fileset, $3: flags for both.  The output must match
# plink2's and differ from the input.
same_family_edit() {
    plink --bfile "$2" --keep-allele-order $3 --make-bed --out tmp_fe19 > /dev/null 2>&1
    $plink2 --bfile "$2" $3 --make-bed --out tmp_fe2 > /dev/null 2>&1
    test -s tmp_fe19.bed && cmp -s tmp_fe19.bed tmp_fe2.bed && cmp -s tmp_fe19.bim tmp_fe2.bim && cmp -s <(spaces tmp_fe19.fam) <(spaces tmp_fe2.fam) && \
        ! { cmp -s tmp_fe19.bed "$2.bed" && cmp -s tmp_fe19.fam "$2.fam"; }
    check "$1" $?
}
same_family_edit "--me matches plink2" tmp_me "--me 0.05 0.015"
test "$(wc -l < tmp_fe19.fam)" -lt 183 && test "$(wc -l < tmp_fe19.bim)" -lt 120
check "--me removes both samples and variants" $?
same_family_edit "--me --me-exclude-one matches plink2" tmp_me "--me 0.05 0.1 --me-exclude-one"
same_family_edit "--set-me-missing matches plink2" tmp_me "--me 1 1 --set-me-missing"

# --make-founders, with some fathers removed: 'first' applies it before
# --remove, and 'require-2-missing' only clears both parents.  All four
# combinations must differ, and match plink2.
awk '$2 ~ /^fa/ && substr($2, 3) % 4 == 1 { print $1, $2 }' tmp_duo.fam > tmp_rm.txt
mf_ok=0
: > tmp_mf_sums.txt
for mode in "" require-2-missing first "require-2-missing first"; do
    plink --bfile tmp_duo --keep-allele-order --remove tmp_rm.txt --make-founders $mode --make-just-fam --out tmp_mf19 > /dev/null 2>&1
    $plink2 --bfile tmp_duo --remove tmp_rm.txt --make-founders $mode --make-bed --out tmp_mf2 > /dev/null 2>&1
    test -s tmp_mf19.fam && cmp -s <(spaces tmp_mf19.fam) <(spaces tmp_mf2.fam) || mf_ok=1
    cksum < tmp_mf19.fam >> tmp_mf_sums.txt
done
test "$(sort -u tmp_mf_sums.txt | wc -l)" -eq 4 || mf_ok=1
check "--make-founders (all modifier combinations) matches plink2" $mf_ok

exit $failed
