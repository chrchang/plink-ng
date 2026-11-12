#!/bin/bash

# PLINK 1.9 flags with no plink2 or PLINK 1.07 reference implementation,
# checked against independent recalculations: --indep-pairphase, --lasso and
# --lasso-select-covars, --set-p/--set-r2/--set-test-lambda (with --set-max),
# --je-cellmin and --distance-wts.
#
# Usage: ./run_tests.sh {plink 1.9 dir}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}

# Exits nonzero unless two numeric files agree entry by entry after the first
# is multiplied by $3 (default 1), to a relative 1e-5 (6-significant-digit
# output).
num_match() {
    paste "$1" "$2" | awk -v f="${3:-1}" '{
        n = NF / 2
        if (n != int(n)) exit 1
        for (k = 1; k <= n; ++k) {
            a = $k * f
            b = $(k + n)
            d = a - b
            if (d < 0) d = -d
            s = (a < 0)? -a : a
            if (d > 1e-5 * s + 1e-9) exit 1
        }
        ++rows
    } END { exit !rows }'
}

# Exits nonzero unless a --lasso report matches lasso.awk output: the same
# nonzero coefficients, to a relative 1e-5.
lasso_match() {
    awk 'NR == FNR { e[$1] = $2; ++ne; next }
    FNR > 1 {
        k = ($1 == "COV")? "COV:" $2 : $2
        if ($4 == 0 || $4 == "NA") next
        if (!(k in e)) exit 1
        d = $4 - e[k]
        if (d < 0) d = -d
        s = ($4 < 0)? -$4 : $4
        if (d > 1e-5 * s + 1e-9) exit 1
        ++no
    } END { exit !(no == ne && ne > 0) }' "$1" "$2"
}

# Exits zero iff the command fails and its output contains the given text.
fails_with() {
    local msg="$1"
    shift
    "$@" > tmp_err.out 2>&1 && return 1
    grep -q -- "$msg" tmp_err.out
}

# 400 unrelated samples: six 15-variant LD blocks on chromosome 1 and 30
# independent variants on chromosome 2, 1% missing calls, and a case/control
# phenotype driven by three of the variants; and an all-homozygous copy.
awk -v n=400 -v n_block=6 -v block_len=15 -v n_indep=30 -v seed=5 -v miss=0.01 -v map=tmp_bl.map -f make_blocks.awk > tmp_bl.ped
awk '{ for (i = 7; i < NF; i += 2) $(i + 1) = $i; print }' tmp_bl.ped > tmp_hz.ped
cp tmp_bl.map tmp_hz.map
for f in bl hz; do
    plink --file tmp_$f --make-bed --out tmp_$f > /dev/null 2>&1
done

# --- --indep-pairphase ---

# Same windowing and greedy pruning as --indep-pairwise, but on r^2 from
# EM-phased haplotype frequencies, which "--r2 dprime" also reports.
# pair_prune.awk reproduces the pruning from an r^2 table; it's validated
# here on --indep-pairwise too (whose results are compared with plink2's
# elsewhere).
plink --bfile tmp_bl --freq counts --out tmp_bl > /dev/null 2>&1
r2_all="--ld-window 99999 --ld-window-kb 99999 --ld-window-r2 0"
plink --bfile tmp_bl --r2 $r2_all --out tmp_blg > /dev/null 2>&1
plink --bfile tmp_bl --r2 dprime $r2_all --out tmp_blh > /dev/null 2>&1
ip_ok=0
pw_ok=0
any_diff=1
for spec in "10 3 0.3" "20 5 0.5" "7 1 0.2" "50 10 0.4" "10 kb 1 0.3" "6 kb 2 0.2" "25 kb 1 0.25"; do
    set -- $spec
    kb=0
    if [[ "$2" == kb ]]; then
        kb=1
        set -- "$1" "$3" "$4"
    fi
    win=$1
    [[ $kb -eq 1 ]] && win="$1 kb"
    rm -f tmp_pp.prune.in tmp_pw.prune.in
    plink --bfile tmp_bl --indep-pairphase $win $2 $3 --out tmp_pp > /dev/null 2>&1
    awk -v kb=$kb -v win=$1 -v step=$2 -v thr=$3 -f pair_prune.awk tmp_bl.frq.counts tmp_blh.ld tmp_bl.bim > tmp_pp_exp.txt || ip_ok=1
    test -s tmp_pp.prune.in && cmp -s tmp_pp.prune.in tmp_pp_exp.txt || ip_ok=1
    plink --bfile tmp_bl --indep-pairwise $win $2 $3 --out tmp_pw > /dev/null 2>&1
    awk -v kb=$kb -v win=$1 -v step=$2 -v thr=$3 -f pair_prune.awk tmp_bl.frq.counts tmp_blg.ld tmp_bl.bim > tmp_pw_exp.txt || pw_ok=1
    test -s tmp_pw.prune.in && cmp -s tmp_pw.prune.in tmp_pw_exp.txt || pw_ok=1
    # Something must be pruned.
    [[ $(wc -l < tmp_pp.prune.in) -lt $(wc -l < tmp_bl.bim) ]] || ip_ok=1
    cmp -s tmp_pp.prune.in tmp_pw.prune.in || any_diff=0
done
check "--indep-pairwise matches the pruning recalculation (validating it)" $pw_ok
check "--indep-pairphase matches the pruning recalculation on haplotype r^2" $ip_ok
check "--indep-pairphase and --indep-pairwise differ on unphased data" $any_diff

# With no heterozygous calls, haplotypes are known, so the two agree.
hz_ok=0
for spec in "10 3 0.3" "7 1 0.2" "8 kb 1 0.25"; do
    plink --bfile tmp_hz --indep-pairwise $spec --out tmp_hzw > /dev/null 2>&1
    plink --bfile tmp_hz --indep-pairphase $spec --out tmp_hzp > /dev/null 2>&1
    test -s tmp_hzp.prune.in && cmp -s tmp_hzw.prune.in tmp_hzp.prune.in || hz_ok=1
done
check "--indep-pairphase equals --indep-pairwise on all-homozygous data" $hz_ok

# --- --lasso, --lasso-select-covars ---

# 200 samples with a quantitative phenotype (and a case/control copy) over
# four 10-variant blocks and 20 independent variants; the last variant is
# made monomorphic.  Three covariates: noise, the phenotype plus noise, and a
# 0/1 value.
for t in qt cc; do
    qt=0
    [[ $t == qt ]] && qt=1
    awk -v n=200 -v n_block=4 -v block_len=10 -v n_indep=20 -v seed=8 -v miss=0.01 -v qt=$qt -v map=tmp_l$t.map -f make_blocks.awk | awk '{ if ($NF != "0") $NF = $(NF - 1) = "A"; print }' > tmp_l$t.ped
    plink --file tmp_l$t --make-bed --out tmp_l$t > /dev/null 2>&1
    plink --bfile tmp_l$t --recode A --out tmp_l$t > /dev/null 2>&1
done
awk -v s=17 'function rnd() { s = (s * 69069 + 1) % 4294967296; return s / 4294967296 }
BEGIN { print "FID IID c1 c2 c3" }
{ printf "%s %s %.3f %.3f %d\n", $1, $2, rnd() * 4, $6 + 3 * (rnd() - 0.5), (rnd() < 0.3) }' tmp_lqt.fam > tmp_lcov.txt
# (c2 last, for --lasso-select-covars c2: unpenalized covariates come first)
awk '{ print $1, $2, $3, $5, $4 }' tmp_lcov.txt > tmp_lcov_c2last.txt

# With an explicit minimum lambda, the result is deterministic and
# lasso.awk reproduces PLINK's coordinate descent.
la_ok=0
plink --bfile tmp_lqt --lasso 0.5 0.02 --out tmp_la1 > /dev/null 2>&1
awk -v minlambda=0.02 -f lasso.awk /dev/null tmp_lqt.raw > tmp_la1_exp.txt
lasso_match tmp_la1_exp.txt tmp_la1.lasso || la_ok=1
plink --bfile tmp_lcc --lasso 0.3 0.03 --out tmp_la2 > /dev/null 2>&1
awk -v cc=1 -v minlambda=0.03 -f lasso.awk /dev/null tmp_lcc.raw > tmp_la2_exp.txt
lasso_match tmp_la2_exp.txt tmp_la2.lasso || la_ok=1
# Some coefficients must be zero, and some nonzero.
n_nz=$(wc -l < tmp_la1_exp.txt)
[[ $n_nz -gt 3 && $n_nz -lt 55 ]] || la_ok=1
check "--lasso matches a coordinate-descent recalculation" $la_ok

# Covariates are unpenalized (and unreported) unless --lasso-select-covars
# names them (or is given without names).
ls_ok=0
plink --bfile tmp_lqt --lasso 0.5 0.02 --covar tmp_lcov.txt --out tmp_ls0 > /dev/null 2>&1
awk -v minlambda=0.02 -v n_unsel=3 -f lasso.awk tmp_lcov.txt tmp_lqt.raw > tmp_ls0_exp.txt
lasso_match tmp_ls0_exp.txt tmp_ls0.lasso || ls_ok=1
! grep -q '^COV' tmp_ls0.lasso || ls_ok=1
plink --bfile tmp_lqt --lasso 0.5 0.02 --covar tmp_lcov.txt --lasso-select-covars --out tmp_ls1 > /dev/null 2>&1
awk -v minlambda=0.02 -v n_unsel=0 -f lasso.awk tmp_lcov.txt tmp_lqt.raw > tmp_ls1_exp.txt
lasso_match tmp_ls1_exp.txt tmp_ls1.lasso || ls_ok=1
plink --bfile tmp_lqt --lasso 0.5 0.02 --covar tmp_lcov.txt --lasso-select-covars c2 --out tmp_ls2 > /dev/null 2>&1
awk -v minlambda=0.02 -v n_unsel=2 -f lasso.awk tmp_lcov_c2last.txt tmp_lqt.raw > tmp_ls2_exp.txt
lasso_match tmp_ls2_exp.txt tmp_ls2.lasso || ls_ok=1
grep -q '^COV	c2	NA	' tmp_ls2.lasso || ls_ok=1
# The unpenalized covariates must change the variant coefficients.
! cmp -s <(grep -v '^COV' tmp_ls1_exp.txt) <(grep -v '^COV' tmp_ls2_exp.txt) || ls_ok=1
check "--lasso covariates and --lasso-select-covars match the recalculation" $ls_ok

# 'report-zeroes' adds the zero coefficients, and NA for monomorphic
# variants.
lz_ok=0
plink --bfile tmp_lqt --lasso 0.5 0.02 report-zeroes --out tmp_lz > /dev/null 2>&1
[[ $(awk 'NR > 1' tmp_lz.lasso | wc -l) -eq $(wc -l < tmp_lqt.bim) ]] || lz_ok=1
cmp -s <(awk 'NR > 1 && $4 != 0 && $4 != "NA"' tmp_lz.lasso) <(awk 'NR > 1' tmp_la1.lasso) || lz_ok=1
[[ $(awk 'NR > 1 && $4 == "NA"' tmp_lz.lasso | wc -l) -eq 1 ]] || lz_ok=1
tail -n 1 tmp_lz.lasso | grep -q '	NA$' || lz_ok=1
check "--lasso report-zeroes" $lz_ok

# The heritability estimate only calibrates the default minimum lambda,
# which scales with sqrt(1 - h2 + 1/n).
lh_ok=0
plink --bfile tmp_lqt --lasso 0.9 0.02 --out tmp_lh > /dev/null 2>&1
cmp -s tmp_la1.lasso tmp_lh.lasso || lh_ok=1
for h in 0.2 0.6; do
    plink --bfile tmp_lqt --lasso $h --seed 3 --out tmp_lh$h > /dev/null 2>&1
done
lam1=$(sed -n 's/^--lasso: using min lambda = \(.*\)\.$/\1/p' tmp_lh0.2.log)
lam2=$(sed -n 's/^--lasso: using min lambda = \(.*\)\.$/\1/p' tmp_lh0.6.log)
awk -v a="$lam1" -v b="$lam2" 'BEGIN {
    r = a / b
    e = sqrt((1 - 0.2 + 1 / 200) / (1 - 0.6 + 1 / 200))
    d = r - e
    if (d < 0) d = -d
    exit !((a > 0) && (d < 1e-5 * e))
}' || lh_ok=1
check "--lasso heritability scales only the default minimum lambda" $lh_ok

lx_ok=0
fails_with "Invalid --lasso heritability estimate '0'" plink --bfile tmp_lqt --lasso 0 --out tmp_lx || lx_ok=1
fails_with "Invalid --lasso heritability estimate '1.5'" plink --bfile tmp_lqt --lasso 1.5 --out tmp_lx || lx_ok=1
fails_with "Invalid --lasso minimum lambda '0'" plink --bfile tmp_lqt --lasso 0.5 0 --out tmp_lx || lx_ok=1
fails_with "min lambda >= max lambda" plink --bfile tmp_lqt --lasso 0.5 100 --out tmp_lx || lx_ok=1
fails_with "lasso-select-covars must be used with --covar" plink --bfile tmp_lqt --lasso 0.5 --lasso-select-covars --out tmp_lx || lx_ok=1
check "--lasso rejects invalid parameters" $lx_ok

# --- --set-p, --set-r2, --set-test-lambda (and --set-max) ---

# Sets: each of the first three chromosome 1 blocks, the other three blocks
# together, chromosome 2, and one set spanning both chromosomes.
printf '1 1000 15000 S1\n1 16000 30000 S2\n1 31000 45000 S3\n1 46000 90000 S4\n2 91000 120000 S5\n1 1000 8000 S6\n2 91000 95000 S6\n' > tmp_sets.txt
plink --bfile tmp_bl --assoc --out tmp_bl > /dev/null 2>&1

# NSIG, ISIG and SNPS are deterministic: recompute them from the --assoc
# chi-squares and --r2 values.  (thr = 1df chi-square at p = --set-p.)
sp_ok=0
for spec in "- 3.841459 - - -" "0.001 10.827566 - - -" "0.01 6.634897 2 - -" "0.05 3.841459 1.5 0.2 3" "0.05 3.841459 - 0.8 8" "0.0001 15.136705 - 0.1 2" "1 0 - 1.5 20" "0.01 6.634897 - 0 -"; do
    set -- $spec
    flags=
    vars="-v thr=$2"
    [[ "$1" != - ]] && flags="$flags --set-p $1"
    [[ "$3" != - ]] && flags="$flags --set-test-lambda $3" && vars="$vars -v lambda=$3"
    [[ "$4" != - ]] && flags="$flags --set-r2 $4" && vars="$vars -v r2max=$4"
    [[ "$5" != - ]] && flags="$flags --set-max $5" && vars="$vars -v set_max=$5"
    # --set-r2 0 means one variant per set.
    [[ "$4" == 0 ]] && vars="-v thr=$2 -v set_max=1"
    rm -f tmp_st.assoc.set.mperm
    plink --bfile tmp_bl --assoc mperm=50 set-test --make-set tmp_sets.txt $flags --seed 1 --out tmp_st > /dev/null 2>&1
    awk $vars -f set_select.awk tmp_sets.txt tmp_bl.bim tmp_bl.assoc tmp_blg.ld > tmp_st_exp.txt || sp_ok=1
    awk 'NR > 1 { print $1, $2, $3, $4, $6 }' tmp_st.assoc.set.mperm | cmp -s - tmp_st_exp.txt || sp_ok=1
done
check "--set-p/--set-r2/--set-test-lambda/--set-max significant variant selection" $sp_ok

# 'write' dumps the LD map.
sw_ok=0
for r in 0.5 0.3; do
    rm -f tmp_sw.ldset
    plink --bfile tmp_bl --assoc mperm=20 set-test --make-set tmp_sets.txt --set-r2 $r write --seed 1 --out tmp_sw > /dev/null 2>&1
    awk -v ldset=1 -v thr=3.841459 -v r2max=$r -f set_select.awk tmp_sets.txt tmp_bl.bim tmp_bl.assoc tmp_blg.ld | sed 's/ *$//' > tmp_sw_exp.txt || sw_ok=1
    test -s tmp_sw_exp.txt && sed 's/ *$//' tmp_sw.ldset | cmp -s - tmp_sw_exp.txt || sw_ok=1
done
check "--set-r2 write matches the --r2 LD map" $sw_ok

# A lambda below 1 is raised to 1, i.e. no correction.
sl_ok=0
plink --bfile tmp_bl --assoc mperm=200 set-test --make-set tmp_sets.txt --seed 1 --out tmp_sl1 > /dev/null 2>&1
plink --bfile tmp_bl --assoc mperm=200 set-test --make-set tmp_sets.txt --set-test-lambda 0.5 --seed 1 --out tmp_sl2 > tmp_sl2.out 2>&1
grep -q 'set-test-lambda parameter set to 1' tmp_sl2.out || sl_ok=1
test -s tmp_sl1.assoc.set.mperm && cmp -s tmp_sl1.assoc.set.mperm tmp_sl2.assoc.set.mperm || sl_ok=1
check "--set-test-lambda below 1 is ignored" $sl_ok

sx_ok=0
fails_with "Invalid --set-p parameter '0'" plink --bfile tmp_bl --assoc mperm=10 set-test --make-set tmp_sets.txt --set-p 0 --out tmp_sx || sx_ok=1
fails_with "Invalid --set-p parameter '1.5'" plink --bfile tmp_bl --assoc mperm=10 set-test --make-set tmp_sets.txt --set-p 1.5 --out tmp_sx || sx_ok=1
fails_with "Invalid --set-r2 parameter '-1'" plink --bfile tmp_bl --assoc mperm=10 set-test --make-set tmp_sets.txt --set-r2 -1 --out tmp_sx || sx_ok=1
fails_with "Invalid --set-r2 parameter 'wrote'" plink --bfile tmp_bl --assoc mperm=10 set-test --make-set tmp_sets.txt --set-r2 0.5 wrote --out tmp_sx || sx_ok=1
fails_with "Invalid --set-test-lambda parameter 'x'" plink --bfile tmp_bl --assoc mperm=10 set-test --make-set tmp_sets.txt --set-test-lambda x --out tmp_sx || sx_ok=1
for f in "--set-p 0.01" "--set-r2 0.3" "--set-test-lambda 1.2"; do
    fails_with "${f%% *} must be used with --set/--make-set" plink --bfile tmp_bl --assoc $f --out tmp_sx || sx_ok=1
done
check "--set-p/--set-r2/--set-test-lambda reject invalid parameters" $sx_ok

# --- --je-cellmin ---

# On the chromosome 2 variants, every pair of variants that passes the
# genotype count prefilter is reported (with --epi1 1), with nan for pairs
# failing the cell count requirement; the others' statistics are unchanged
# from an unfiltered run.
plink --bfile tmp_bl --chr 2 --recode A --out tmp_bl2 > /dev/null 2>&1
plink --bfile tmp_bl --chr 2 --fast-epistasis joint-effects --epi1 1 --je-cellmin 0 --out tmp_je0 > /dev/null 2>&1
je_ok=0
for c in 0 1 2 3 5 8 default; do
    flag="--je-cellmin $c"
    cm=$c
    if [[ $c == default ]]; then
        flag=
        cm=5
    fi
    rm -f tmp_je.epi.cc
    plink --bfile tmp_bl --chr 2 --fast-epistasis joint-effects --epi1 1 $flag --out tmp_je > /dev/null 2>&1
    awk -v cellmin=$cm -f je_cells.awk tmp_bl2.raw > tmp_je_exp.txt
    skipped=$(sed -n 's/^--fast-epistasis: Skipping \([0-9]*\) sites.*/\1/p' tmp_je.log)
    valid=$(sed -n 's/^\([0-9]*\) valid tests\{0,1\} performed.*/\1/p' tmp_je.log)
    awk -v skipped="${skipped:-0}" -v valid="$valid" '
    FILENAME == ARGV[1] { if (FNR > 1) stat0[$2, $4] = $5; next }
    FILENAME == ARGV[2] {
        if ($1 == "SITES") { if ($2 != skipped) exit 1; next }
        ++ne
        pair[ne] = $1 " " $2
        tested[ne] = $3
        nvalid += $3
        next
    }
    FNR > 1 {
        ++no
        if (pair[no] != $2 " " $4) exit 1
        if (tested[no]) {
            if ($5 != stat0[$2, $4]) exit 1
        } else if ($5 != "nan") {
            exit 1
        }
    }
    END { exit !((no == ne) && (ne > 0) && (nvalid == valid)) }' tmp_je0.epi.cc tmp_je_exp.txt tmp_je.epi.cc || je_ok=1
done
check "--je-cellmin skips the variants and pairs with too few observations" $je_ok

jx_ok=0
fails_with "Too few cases or controls for --je-cellmin 30" plink --bfile tmp_bl --chr 2 --fast-epistasis joint-effects --je-cellmin 30 --out tmp_jx || jx_ok=1
fails_with 'je-cellmin must be used with "--fast-epistasis joint-effects"' plink --bfile tmp_bl --chr 2 --fast-epistasis --je-cellmin 2 --out tmp_jx || jx_ok=1
fails_with "Invalid --je-cellmin parameter 'x'" plink --bfile tmp_bl --chr 2 --fast-epistasis joint-effects --je-cellmin x --out tmp_jx || jx_ok=1
check "--je-cellmin rejects invalid parameters" $jx_ok

# --- --distance-wts ---

dist="--distance square allele-ct ibs 1-ibs"

# Without missing calls, the allele-count distance is the weighted sum of
# genotype differences, and 1-ibs is that divided by twice the weight sum.
awk -v n=60 -v n_block=3 -v block_len=10 -v n_indep=20 -v seed=4 -v miss=0 -v map=tmp_dn.map -f make_blocks.awk > tmp_dn.ped
plink --file tmp_dn --make-bed --out tmp_dn > /dev/null 2>&1
plink --bfile tmp_dn --recode A --out tmp_dn > /dev/null 2>&1
awk 'BEGIN { print "ID WT" } { print $2, 0.3 + (NR % 5) * 0.4 }' tmp_dn.bim > tmp_dnw.txt
plink --bfile tmp_dn $dist --distance-wts tmp_dnw.txt --out tmp_dnd > /dev/null 2>&1
awk 'NR == FNR { if (FNR > 1) { w[FNR - 1] = $2; wsum += $2 }; next }
FNR > 1 { ++n; for (j = 7; j <= NF; ++j) g[n, j - 6] = $j; m = NF - 6 }
END {
    for (a = 1; a <= n; ++a) {
        l1 = l2 = l3 = ""
        for (b = 1; b <= n; ++b) {
            d = 0
            for (j = 1; j <= m; ++j) {
                x = g[a, j] - g[b, j]
                d += w[j] * ((x < 0)? -x : x)
            }
            sep = (b > 1)? "\t" : ""
            l1 = l1 sep d
            l2 = l2 sep (1 - d / (2 * wsum))
            l3 = l3 sep (d / (2 * wsum))
        }
        print l1 > "tmp_dne.dist"
        print l2 > "tmp_dne.mibs"
        print l3 > "tmp_dne.mdist"
    }
}' tmp_dnw.txt tmp_dn.raw
dw_ok=0
for e in dist mibs mdist; do
    num_match tmp_dnd.$e tmp_dne.$e || dw_ok=1
done
check "--distance-wts weights genotype differences" $dw_ok

# With missing calls (whose correction depends on the weights):
# - all-1 weights give the unweighted matrices;
# - scaling the weights scales the allele-count distances only;
# - a weight of 0, or no weight, is the same as excluding the variant;
# - 'noheader' reads the first line as a weight.
dm_ok=0
awk 'BEGIN { print "ID WT" } { print $2, 1 }' tmp_bl.bim > tmp_w1.txt
awk 'BEGIN { print "ID WT" } { print $2, 0.3 + (NR % 5) * 0.4 }' tmp_bl.bim > tmp_wa.txt
awk 'NR == 1 { print; next } { print $1, $2 * 3 }' tmp_wa.txt > tmp_wb.txt
awk 'NR == 1 { print; next } { print $1, (NR % 9 == 0)? 0 : $2 }' tmp_wa.txt > tmp_wz.txt
awk 'NR > 1 && $2 == 0 { print $1 }' tmp_wz.txt > tmp_wz_ids.txt
awk 'NR == 1 || $2 != 0' tmp_wz.txt > tmp_wn.txt
tail -n +2 tmp_wa.txt > tmp_wa_noheader.txt
plink --bfile tmp_bl $dist --out tmp_du > /dev/null 2>&1
for w in 1 a b z n; do
    plink --bfile tmp_bl $dist --distance-wts tmp_w$w.txt --out tmp_d$w > /dev/null 2>&1
done
plink --bfile tmp_bl $dist --distance-wts tmp_wa.txt --exclude tmp_wz_ids.txt --out tmp_dx > /dev/null 2>&1
plink --bfile tmp_bl $dist --distance-wts tmp_wa_noheader.txt noheader --out tmp_dh > /dev/null 2>&1
for e in dist mibs mdist; do
    num_match tmp_du.$e tmp_d1.$e || dm_ok=1
    num_match tmp_dz.$e tmp_dx.$e || dm_ok=1
    num_match tmp_dn.$e tmp_dx.$e || dm_ok=1
    ! num_match tmp_dz.$e tmp_da.$e || dm_ok=1
    ! num_match tmp_da.$e tmp_du.$e || dm_ok=1
    cmp -s tmp_da.$e tmp_dh.$e || dm_ok=1
done
num_match tmp_da.dist tmp_db.dist 3 || dm_ok=1
num_match tmp_da.mibs tmp_db.mibs || dm_ok=1
num_match tmp_da.mdist tmp_db.mdist || dm_ok=1
check "--distance-wts with missing calls: unit, scaled, zero and absent weights; noheader" $dm_ok

# "exp=<x>" weights each variant by (2q(1-q))^(-x); with --read-freq, q is
# known exactly, so a weight file can match it.
de_ok=0
plink --bfile tmp_bl --freq --out tmp_blf > /dev/null 2>&1
awk 'NR == 1 { print; next } { $5 = sprintf("%.2f", $5); print }' tmp_blf.frq > tmp_blf2.frq
awk 'NR > 1 { q = $5; printf "%s %.17g\n", $2, (2 * q * (1 - q)) ^ (-0.7) }' tmp_blf2.frq > tmp_wexp.txt
plink --bfile tmp_bl --read-freq tmp_blf2.frq $dist --distance-wts exp=0.7 --out tmp_de > /dev/null 2>&1
plink --bfile tmp_bl --read-freq tmp_blf2.frq $dist --distance-wts tmp_wexp.txt noheader --out tmp_df > /dev/null 2>&1
for e in dist mibs mdist; do
    num_match tmp_de.$e tmp_df.$e || de_ok=1
done
! num_match tmp_de.mibs tmp_du.mibs || de_ok=1
check "--distance-wts exp=<x> matches explicit (2q(1-q))^-x weights" $de_ok

dx_ok=0
printf 'ID WT\nv1 -1\n' > tmp_wneg.txt
printf 'ID WT\nv1 abc\n' > tmp_wbad.txt
printf 'ID WT\nv1 1\nv1 2\n' > tmp_wdup.txt
printf 'ID WT\nv1 0\nv2 0\n' > tmp_wzero.txt
fails_with "Invalid weight" plink --bfile tmp_bl --distance --distance-wts tmp_wneg.txt --out tmp_dxx || dx_ok=1
fails_with "Invalid weight" plink --bfile tmp_bl --distance --distance-wts tmp_wbad.txt --out tmp_dxx || dx_ok=1
fails_with "Duplicate variant ID 'v1'" plink --bfile tmp_bl --distance --distance-wts tmp_wdup.txt --out tmp_dxx || dx_ok=1
fails_with "No valid nonzero entries" plink --bfile tmp_bl --distance --distance-wts tmp_wzero.txt --out tmp_dxx || dx_ok=1
fails_with "Invalid --distance-wts exponent" plink --bfile tmp_bl --distance --distance-wts exp=x --out tmp_dxx || dx_ok=1
fails_with "cannot be used with --distance's" plink --bfile tmp_bl --distance flat-missing --distance-wts tmp_wa.txt --out tmp_dxx || dx_ok=1
fails_with "must be used with --distance" plink --bfile tmp_bl --freq --distance-wts tmp_wa.txt --out tmp_dxx || dx_ok=1
check "--distance-wts rejects invalid weights and combinations" $dx_ok

exit $failed
