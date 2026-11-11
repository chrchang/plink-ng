#!/bin/bash

# PLINK 1.9 flags with no plink2 equivalent, checked against PLINK 1.07:
# --cluster constraints, family/stratified association (--mh --family,
# --mh2, --loop-assoc, --qfam*), --genome --min/--max, --neighbour,
# --ibs-test, --fast-epistasis case-only --gap, --annotate --border,
# --write-dosage, --with-phenotype, the --clump options plink2 lacks, and
# --simulate-label.  A few related checks that don't need PLINK 1.07 come
# first.
#
# Usage: ./run_tests.sh {plink 1.9 dir}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  The PLINK 1.07 binary is taken from $PLINK107,
# defaulting to the build in the repository's 1.07/ directory; the checks
# that need it are skipped if it isn't there.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}
plink107="${PLINK107:-../../../1.07/plink}"
compare_cols="../compare_cols.awk"
compare_tokens="../compare_tokens.awk"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}
p107() {
    "$plink107" --noweb --silent "$@" > /dev/null 2>&1
}

# --- Data ---

# 61 samples in families of 4, 500 variants: for clustering, --genome,
# --neighbour and --ibs-test.
plink --dummy 61 500 0.03 acgt --seed 23 --make-bed --out tmp_rel > /dev/null 2>&1
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_rel.bim > tmp_rel.bim.new
mv tmp_rel.bim.new tmp_rel.bim
awk '{ $1 = "f" int((NR - 1) / 4); print }' tmp_rel.fam > tmp_rel.fam.new
mv tmp_rel.fam.new tmp_rel.fam
# The LD data of TEST_19_LD_FAMILY_VS_PLINK2: 203 samples, 150 variants in LD
# blocks, for clumping.
awk -v n=203 -v m=150 -v seed=17 -f ../TEST_19_LD_FAMILY_VS_PLINK2/make_ld.awk > tmp_ld.ped
awk 'BEGIN { for (j = 1; j <= 150; ++j) print (j <= 100)? 1 : 2, "v" j, j * 0.05, 2000 * j }' > tmp_ld.map
plink --file tmp_ld --keep-allele-order --make-bed --out tmp_ld > /dev/null 2>&1

# --- Checks that don't need PLINK 1.07 ---

# PLINK 1.07 ignores --max without --min (and --min --max is compared below),
# so here the PLINK 1.9 report must hold exactly the pairs with PI_HAT <= max.
plink --bfile tmp_rel --genome --out tmp_gall > /dev/null 2>&1
plink --bfile tmp_rel --genome --max 0.02 --out tmp_gmax > /dev/null 2>&1
cmp -s <(awk 'NR == 1 || $10 <= 0.02' tmp_gall.genome) tmp_gmax.genome && test "$(wc -l < tmp_gmax.genome)" -lt "$(wc -l < tmp_gall.genome)"
check "--genome --max keeps the pairs with PI_HAT <= max" $?

# --simulate-prevalence: PLINK 1.9 treats the simulation file's odds ratios as
# odds ratios (PLINK 1.07 treats them as relative risks, which only agree at
# low prevalence).  The mean risk allele frequency over 100 variants (2000
# cases, 2000 controls, risk allele frequency 0.2, odds ratio 2.5) must be
# close to the expected value under a logistic model, at two prevalences.
printf '100 snpA 0.2 0.2 2.5 mult\n' > tmp_sim.txt
prev_ok=0
for k in 0.01 0.4; do
    rm -f tmp_simp.assoc
    plink --simulate tmp_sim.txt --simulate-ncases 2000 --simulate-ncontrols 2000 --simulate-prevalence $k --seed 1 --assoc --out tmp_simp > /dev/null 2>&1
    awk -v k=$k 'BEGIN {
        p = 0.2; odds_ratio = 2.5; lo = 1e-9; hi = 1e9
        for (it = 0; it < 200; ++it) {
            o0 = sqrt(lo * hi); m = 0
            for (g = 0; g <= 2; ++g) { f = (g == 0)? (1 - p)^2 : (g == 1)? 2 * p * (1 - p) : p^2; o = o0 * odds_ratio^g; m += f * o / (1 + o) }
            if (m > k) hi = o0; else lo = o0
        }
        for (g = 0; g <= 2; ++g) {
            f = (g == 0)? (1 - p)^2 : (g == 1)? 2 * p * (1 - p) : p^2; o = o0 * odds_ratio^g; pen = o / (1 + o)
            ca += f * pen * g; ta += f * pen; cu += f * (1 - pen) * g; tu += f * (1 - pen)
        }
        want_a = ca / ta / 2; want_u = cu / tu / 2
    }
    NR > 1 { if ($4 == "D") { fa += $5; fu += $6 } else { fa += 1 - $5; fu += 1 - $6 } ++n }
    END {
        da = fa / n - want_a; du = fu / n - want_u
        exit (n != 100 || da > 0.01 || da < -0.01 || du > 0.01 || du < -0.01)
    }' tmp_simp.assoc || prev_ok=1
done
check "--simulate-prevalence matches the expected allele frequencies" $prev_ok

# --clump-replicate keeps a clump iff its p2-significant results, including the
# index variant's own, come from at least two files, as the PLINK 1.07
# documentation describes.  (PLINK 1.07's implementation differs in one case:
# it needs at least two secondary results, so it drops a clump whose only
# support is a single result in another file.)
# The index variant is v12 (p = 1e-6) in file 1; v15 is in strong LD with it.
# Prints a results file with the given p-values ($1, e.g. "v12=1e-3") and 0.5
# elsewhere.
clump_rep() {
    awk -v spec="$1" 'BEGIN {
        n = split(spec, a, ",")
        for (i = 1; i <= n; ++i) {
            split(a[i], kv, "=")
            p[kv[1]] = kv[2]
        }
        print "CHR SNP BP P"
    }
    { print $1, $2, $4, ($2 in p)? p[$2] : 0.5 }' tmp_ld.bim
}
rep_ok=0
# file 1 results ; file 2 results ; 1 if the clump must be kept
for case in "v12=1e-6;v12=1e-3;1" "v12=1e-6;v12=1e-3,v15=1e-3;1" "v12=1e-6,v15=1e-3;v12=1e-3;1" \
    "v12=1e-6,v15=1e-3;v99=1e-3;0" "v12=1e-6;v99=1e-3;0"; do
    IFS=';' read -r spec1 spec2 want <<< "$case"
    clump_rep "$spec1" > tmp_rep1.txt
    clump_rep "$spec2" > tmp_rep2.txt
    rm -f tmp_rep.clumped tmp_norep.clumped
    plink --bfile tmp_ld --clump tmp_rep1.txt,tmp_rep2.txt --out tmp_norep > /dev/null 2>&1
    plink --bfile tmp_ld --clump tmp_rep1.txt,tmp_rep2.txt --clump-replicate --out tmp_rep > /dev/null 2>&1
    grep -q ' v12 ' tmp_norep.clumped || rep_ok=1
    kept=0
    if [[ -s tmp_rep.clumped ]] && grep -q ' v12 ' tmp_rep.clumped; then
        kept=1
    fi
    test $kept -eq "$want" || rep_ok=1
done
check "--clump-replicate needs p2-significant results from two files" $rep_ok

# --- Checks against PLINK 1.07 ---

checks=(
    "--cluster --mc matches PLINK 1.07"
    "--cluster --ppc matches PLINK 1.07"
    "--cluster --ppc --ppc-gap matches PLINK 1.07"
    "--cluster --qmatch --qt matches PLINK 1.07"
    "--cluster --ibm matches PLINK 1.07"
    "--cluster --read-genome matches PLINK 1.07"
    "--mh --family matches PLINK 1.07"
    "--mh2 matches PLINK 1.07"
    "--loop-assoc matches PLINK 1.07"
    "--qfam matches PLINK 1.07"
    "--qfam-parents matches PLINK 1.07"
    "--qfam-between matches PLINK 1.07"
    "--qfam-total matches PLINK 1.07"
    "--genome --min matches PLINK 1.07"
    "--genome --min --max matches PLINK 1.07"
    "--neighbour matches PLINK 1.07"
    "--ibs-test statistics match PLINK 1.07"
    "--fast-epistasis case-only --gap matches PLINK 1.07"
    "--annotate --border matches PLINK 1.07"
    "--write-dosage matches PLINK 1.07"
    "--write-covar --with-phenotype matches PLINK 1.07"
    "--clump-verbose matches PLINK 1.07"
    "--clump-verbose --clump-annotate matches PLINK 1.07"
    "--clump-best matches PLINK 1.07"
    "--clump-replicate matches PLINK 1.07"
    "--clump-index-first matches PLINK 1.07"
    "--clump-best with two files matches PLINK 1.07"
    "--clump-best --clump-replicate matches PLINK 1.07"
    "--simulate-label matches PLINK 1.07"
)
if [[ ! -x "$plink107" ]]; then
    echo "note: PLINK 1.07 not found at $plink107; skipping checks that need it"
    for c in "${checks[@]}"; do
        echo "skipped: $c"
    done
    exit $failed
fi
check_idx=0
# Records the next check in ${checks[@]}.
next_check() {
    check "${checks[$check_idx]}" $1
    check_idx=$((check_idx + 1))
}

# --- Clustering constraints ---

# $1: input fileset, $2: flags for both.  PLINK 1.9 uses 'old-tiebreaks', and
# all three .cluster files must match.  The constraint must also change PLINK
# 1.9's solution.
plink --bfile tmp_rel --cluster old-tiebreaks --out tmp_cl_base > /dev/null 2>&1
same_cluster() {
    rm -f tmp_cl19.cluster* tmp_cl107.cluster*
    p107 --bfile "$1" --cluster $2 --out tmp_cl107
    plink --bfile "$1" --cluster old-tiebreaks $2 --out tmp_cl19 > /dev/null 2>&1
    local ok=0
    for f in 1 2 3; do
        test -s tmp_cl19.cluster$f && cmp -s tmp_cl19.cluster$f tmp_cl107.cluster$f || ok=1
    done
    plink --bfile "$1" --cluster old-tiebreaks --out tmp_cl_unc > /dev/null 2>&1
    ! cmp -s tmp_cl19.cluster1 tmp_cl_unc.cluster1 || ok=1
    next_check $ok
}
same_cluster tmp_rel "--mc 3"
same_cluster tmp_rel "--ppc 0.5"
same_cluster tmp_rel "--ppc 0.5 --ppc-gap 100"
awk '{ print $1, $2, (NR % 7) / 2 }' tmp_rel.fam > tmp_qmatch.txt
printf '1\n' > tmp_qt.txt
same_cluster tmp_rel "--qmatch tmp_qmatch.txt --qt tmp_qt.txt"
# For --ibm, families get one of three missingness patterns (the first or
# second 200 variants missing, or neither).
plink --bfile tmp_rel --keep-allele-order --recode --out tmp_relp > /dev/null 2>&1
awk '{
    f = substr($1, 2) + 0
    if (f % 3 == 0) for (i = 7; i <= 406; ++i) $i = 0
    else if (f % 3 == 1) for (i = 407; i <= 806; ++i) $i = 0
    print
}' tmp_relp.ped > tmp_im.ped
cp tmp_relp.map tmp_im.map
plink --file tmp_im --keep-allele-order --make-bed --out tmp_im > /dev/null 2>&1
same_cluster tmp_im "--ibm 0.8"
# The .genome file is from the first half of the variants only, so it gives
# different distances.
plink --bfile tmp_rel --snps snp0-snp249 --genome --out tmp_gpart > /dev/null 2>&1
same_cluster tmp_rel "--read-genome tmp_gpart.genome"

# --- Family and stratified association ---

# 513 cases and controls in families of 9, 200 variants, and 4 strata.
plink --dummy 513 200 0.02 acgt --seed 29 --make-bed --out tmp_cc > /dev/null 2>&1
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_cc.bim > tmp_cc.bim.new
mv tmp_cc.bim.new tmp_cc.bim
awk '{ $1 = "f" int((NR - 1) / 9); print }' tmp_cc.fam > tmp_cc.fam.new
mv tmp_cc.fam.new tmp_cc.fam
awk '{ print $1, $2, "s" ((NR - 1) % 4) }' tmp_cc.fam > tmp_strata.txt

p107 --bfile tmp_cc --mh --family --out tmp_mhf107
plink --bfile tmp_cc --mh --family --out tmp_mhf19 > /dev/null 2>&1
test -s tmp_mhf19.cmh && cmp -s tmp_mhf19.cmh tmp_mhf107.cmh
next_check $?

# PLINK 1.9 adds a DF column and prints more digits.
p107 --bfile tmp_cc --within tmp_strata.txt --mh2 --out tmp_mh2_107
plink --bfile tmp_cc --within tmp_strata.txt --mh2 --out tmp_mh2_19 > /dev/null 2>&1
awk -v c1=1,2,3,5 -v c2=1,2,3,4 -f $compare_cols tmp_mh2_19.cmh2 tmp_mh2_107.cmh2
next_check $?

p107 --bfile tmp_cc --loop-assoc tmp_strata.txt --assoc --out tmp_la107
plink --bfile tmp_cc --loop-assoc tmp_strata.txt --assoc --out tmp_la19 > /dev/null 2>&1
la_ok=0
for s in s0 s1 s2 s3; do
    test -s tmp_la19.$s.assoc && cmp -s tmp_la19.$s.assoc tmp_la107.$s.assoc || la_ok=1
done
next_check $la_ok

# QFAM: 80 trios with a quantitative phenotype.  The empirical p-values
# depend on each version's permutations, so the other columns are compared.
# (PLINK 1.9 names the p-value column RAW_P.)
awk -v n_trio=80 -v n_var=60 -v seed=37 -f ../TEST_19_TDT_EPI_CLUSTER/make_trios.awk | \
    awk 'BEGIN { s = 5 } { s = (s * 69069 + 1) % 4294967296; $6 = int(s / 65536) % 1000 / 100; print }' > tmp_qf.ped
awk 'BEGIN { for (j = 1; j <= 60; ++j) print 1, "v" j, 0, 1000 * j }' > tmp_qf.map
plink --file tmp_qf --make-bed --out tmp_qf > /dev/null 2>&1
for m in qfam:within qfam-parents:parents qfam-between:between qfam-total:total; do
    flag=${m%%:*}
    ext=${m##*:}
    rm -f tmp_qf19.qfam.$ext tmp_qf107.qfam.$ext
    p107 --bfile tmp_qf --$flag --mperm 20 --out tmp_qf107
    plink --bfile tmp_qf --$flag mperm=20 --out tmp_qf19 > /dev/null 2>&1
    awk -v c1=1,2,3,4,5,6,7,8,9 -v c2=1,2,3,4,5,6,7,8,9 -f $compare_cols tmp_qf19.qfam.$ext tmp_qf107.qfam.$ext && test "$(wc -l < tmp_qf19.qfam.$ext)" -eq 61
    next_check $?
done

# --- Relationships ---

p107 --bfile tmp_rel --genome --min 0.01 --out tmp_gmin107
plink --bfile tmp_rel --genome --min 0.01 --out tmp_gmin19 > /dev/null 2>&1
test -s tmp_gmin19.genome && cmp -s tmp_gmin19.genome tmp_gmin107.genome && test "$(wc -l < tmp_gmin19.genome)" -lt "$(wc -l < tmp_gall.genome)"
next_check $?
p107 --bfile tmp_rel --genome --min 0.005 --max 0.03 --out tmp_gmm107
plink --bfile tmp_rel --genome --min 0.005 --max 0.03 --out tmp_gmm19 > /dev/null 2>&1
test -s tmp_gmm19.genome && cmp -s tmp_gmm19.genome tmp_gmm107.genome && test "$(wc -l < tmp_gmm19.genome)" -lt "$(wc -l < tmp_gmin19.genome)"
next_check $?

# --neighbour: PLINK 1.9 leaves out PLINK 1.07's PROP_DIFF column.  Ranks,
# distances and Z scores must match row by row; the neighbors themselves may
# be listed in a different order when tied, and may differ only when tied
# with the last reported distance.
p107 --bfile tmp_rel --neighbour 1 5 --out tmp_nb107
plink --bfile tmp_rel --neighbour 1 5 --out tmp_nb19 > /dev/null 2>&1
nb_ranks() {
    awk 'NR > 1 { print $1, $2, $3, $4, $5 }' "$1"
}
nb_sets() {
    awk 'NR > 1 {
        k = $1 " " $2
        row[k, $3] = $4 " " $6 " " $7
        d[k, $3] = $4
        if ($3 > maxnn) maxnn = $3
        keys[k] = 1
    }
    END { for (k in keys) for (r = 1; r <= maxnn; ++r) if (d[k, r] != d[k, maxnn]) print k, row[k, r] }' "$1" | LC_ALL=C sort
}
test -s tmp_nb19.nearest && cmp -s <(nb_ranks tmp_nb19.nearest) <(nb_ranks tmp_nb107.nearest) && cmp -s <(nb_sets tmp_nb19.nearest) <(nb_sets tmp_nb107.nearest)
next_check $?

# --ibs-test: the p-values come from permutations, so the IBS means and SDs
# (from the logs) are compared.  PLINK 1.07 labels the case group "(2)".
"$plink107" --noweb --bfile tmp_rel --ibs-test --out tmp_ibs107 > /dev/null 2>&1
plink --bfile tmp_rel --ibs-test --out tmp_ibs19 > /dev/null 2>&1
ibs_stats() {
    grep -E 'Between-group IBS|In-group \((case|2)\) IBS' "$1" | sed 's/.*= *//; s/,/ /'
}
awk -v c1=1,2 -v c2=1,2 -v skip1=0 -v skip2=0 -f $compare_cols <(ibs_stats tmp_ibs19.log) <(ibs_stats tmp_ibs107.log) && test "$(ibs_stats tmp_ibs19.log | wc -l)" -eq 2
next_check $?

# --- Epistasis, annotation, dosage and covariate output ---

# --fast-epistasis case-only skips pairs closer than --gap (default 1 Mb), so
# the variants are 20 kb apart on two chromosomes.
awk 'BEGIN { OFS = "\t" } { $1 = (NR <= 100)? 1 : 2; $4 = 20000 * NR; print }' tmp_cc.bim > tmp_co.bim
cp tmp_cc.bed tmp_co.bed
cp tmp_cc.fam tmp_co.fam
gap_ok=0
: > tmp_co_sums.txt
for gap in "" "--gap 500" "--gap 2500"; do
    rm -f tmp_co19.epi.co* tmp_co107.epi.co*
    p107 --bfile tmp_co --fast-epistasis --case-only $gap --epi1 0.05 --out tmp_co107
    plink --bfile tmp_co --fast-epistasis no-ueki case-only $gap --epi1 0.05 --out tmp_co19 > /dev/null 2>&1
    for e in epi.co epi.co.summary; do
        test -s tmp_co19.$e && cmp -s tmp_co19.$e tmp_co107.$e || gap_ok=1
    done
    cksum < tmp_co19.epi.co >> tmp_co_sums.txt
done
test "$(sort -u tmp_co_sums.txt | wc -l)" -eq 3 || gap_ok=1
next_check $gap_ok

# --annotate with GB a single position and GC between variants.
awk 'BEGIN { print "CHR SNP BP P" } { print $1, $2, $4, NR / 1000 }' tmp_rel.bim | head -41 > tmp_annot_in.txt
printf '1 3000 7000 GA\n1 12000 12000 GB\n1 20500 26000 GC\n' > tmp_annot_ranges.txt
an_ok=0
: > tmp_an_cts.txt
for border in "" "--border 1" "--border 2"; do
    rm -f tmp_an19.annot tmp_an107.annot
    p107 --annotate tmp_annot_in.txt ranges=tmp_annot_ranges.txt $border --out tmp_an107
    plink --annotate tmp_annot_in.txt ranges=tmp_annot_ranges.txt $border --out tmp_an19 > /dev/null 2>&1
    test -s tmp_an19.annot && cmp -s <(tr -s ' \t' ' ' < tmp_an19.annot) <(tr -s ' \t' ' ' < tmp_an107.annot) || an_ok=1
    grep -c 'G[ABC]' tmp_an19.annot >> tmp_an_cts.txt
done
test "$(sort -u tmp_an_cts.txt | wc -l)" -eq 3 || an_ok=1
next_check $an_ok

# --write-dosage merges two variant batches of a format=1 dosage file.
head -37 tmp_rel.fam > tmp_dose.fam
awk 'BEGIN { s = 5; ORS = "" }
  NR == FNR { fid[NR] = $1; iid[NR] = $2; n = NR; next }
  FNR == 1 { print "SNP A1 A2"; for (i = 1; i <= n; ++i) print " " fid[i] " " iid[i]; print "\n" }
  { print "v" FNR " A C"; for (i = 1; i <= n; ++i) { s = (s * 69069 + 1) % 4294967296; print " " int(s / 65536) % 2001 / 1000 } print "\n" }' tmp_dose.fam <(seq 1 30) > tmp_dose_all.txt
head -16 tmp_dose_all.txt > tmp_dose1.txt
{ head -1 tmp_dose_all.txt; tail -n +17 tmp_dose_all.txt; } > tmp_dose2.txt
printf '1 tmp_dose1.txt\n2 tmp_dose2.txt\n' > tmp_dose_list.txt
p107 --dosage tmp_dose_list.txt list format=1 --fam tmp_dose.fam --write-dosage --out tmp_wd107
plink --dosage tmp_dose_list.txt list format=1 --fam tmp_dose.fam --write-dosage --out tmp_wd19 > /dev/null 2>&1
test -s tmp_wd19.out.dosage && cmp -s tmp_wd19.out.dosage tmp_wd107.out.dosage && test "$(wc -l < tmp_wd19.out.dosage)" -eq 31
next_check $?

# --with-phenotype (PLINK 1.07 has none of PLINK 1.9's modifiers).
awk '{ if (NR % 5 == 0) { $3 = "per1"; $4 = "per2" } print }' tmp_rel.fam > tmp_wp.fam
cp tmp_rel.bed tmp_wp.bed
cp tmp_rel.bim tmp_wp.bim
awk 'BEGIN { print "FID IID CA CB" } { print $1, $2, NR / 10, NR % 3 }' tmp_wp.fam > tmp_wp_covar.txt
p107 --bfile tmp_wp --covar tmp_wp_covar.txt --write-covar --with-phenotype --out tmp_wp107
plink --bfile tmp_wp --covar tmp_wp_covar.txt --write-covar --with-phenotype --out tmp_wp19 > /dev/null 2>&1
test -s tmp_wp19.cov && cmp -s <(tr -s ' \t' ' ' < tmp_wp19.cov | sed 's/ $//') <(tr -s ' \t' ' ' < tmp_wp107.cov | sed 's/ $//') && \
    test "$(head -1 tmp_wp19.cov | wc -w)" -gt 4
next_check $?

# --- Clump options plink2 lacks ---

# The LD data and p-values of TEST_19_LD_FAMILY_VS_PLINK2, plus a BETA column
# for --clump-annotate.  The reports are compared token by token, since the
# two versions round some p-values differently.  A second results file has
# its strongest signals at other variants, for --clump-index-first; a third
# is for --clump-replicate, using data where PLINK 1.07 follows its
# documented rule (see the PLINK 1.9 --clump-replicate check above).
awk 'BEGIN { s = 23; print "CHR SNP BP P PALT MYID BETA" } {
    s = (s * 69069 + 1) % 4294967296; u = int(s / 65536) % 10000 / 10000
    k = NR % 10
    p = (k == 2 || k == 3)? u / 1e5 : (k == 4)? u / 1e3 : (NR % 7 == 0)? u / 50 : u
    s = (s * 69069 + 1) % 4294967296; q = int(s / 65536) % 10000 / 1e6
    print $1, $2, $4, p, q, $2, NR / 100
}' tmp_ld.bim > tmp_res1.txt
awk 'BEGIN { s = 3 } NR == 1 { print; next } { s = (s * 69069 + 1) % 4294967296; u = int(s / 65536) % 10000 / 10000; $4 = (NR % 10 == 7)? u * 1e-7 : (NR % 10 == 3 || NR % 10 == 4)? u * 1e-4 : u; print }' tmp_res1.txt > tmp_res2.txt
awk 'BEGIN { s = 3 } NR == 1 { print; next } { s = (s * 69069 + 1) % 4294967296; $4 = (NR % 10 == 3 || NR % 10 == 4)? 1e-4 * int(s / 65536) % 100 / 100 : int(s / 65536) % 10000 / 10000; print }' tmp_res1.txt > tmp_res3.txt
# $1: flags for both, $2: extensions to compare, $3: PLINK 1.9 report the
# flags must change (if nonempty).
same_clump() {
    local ok=0
    rm -f tmp_cl19.clumped* tmp_cl107.clumped*
    p107 --bfile tmp_ld $1 --out tmp_cl107
    plink --bfile tmp_ld $1 --out tmp_cl19 > /dev/null 2>&1
    for e in $2; do
        awk -f $compare_tokens tmp_cl19.$e tmp_cl107.$e 2> /dev/null || ok=1
    done
    if [[ -n "${3:-}" ]] && cmp -s tmp_cl19.clumped "$3"; then
        ok=1
    fi
    next_check $ok
}
plink --bfile tmp_ld --clump tmp_res1.txt --out tmp_clbase > /dev/null 2>&1
plink --bfile tmp_ld --clump tmp_res1.txt --clump-verbose --out tmp_clverb > /dev/null 2>&1
plink --bfile tmp_ld --clump tmp_res1.txt,tmp_res2.txt --out tmp_cltwo > /dev/null 2>&1
plink --bfile tmp_ld --clump tmp_res1.txt,tmp_res3.txt --out tmp_clthree > /dev/null 2>&1
same_clump "--clump tmp_res1.txt --clump-verbose" clumped tmp_clbase.clumped
same_clump "--clump tmp_res1.txt --clump-verbose --clump-annotate BETA" clumped tmp_clverb.clumped
same_clump "--clump tmp_res1.txt --clump-best" "clumped clumped.best"
# PLINK 1.07 drops clumps whose only support is a single secondary result in a
# file other than the index variant's (see above), so those are removed from
# PLINK 1.9's report before comparing; there must be some.
rm -f tmp_cl19.clumped tmp_cl107.clumped
p107 --bfile tmp_ld --clump tmp_res1.txt,tmp_res3.txt --clump-replicate --out tmp_cl107
plink --bfile tmp_ld --clump tmp_res1.txt,tmp_res3.txt --clump-replicate --out tmp_cl19 > /dev/null 2>&1
awk 'NF > 5 && NR > 1 { n = split($12, sp, ","); if (n == 1 && sp[1] ~ /\([0-9]+\)$/) { f = sp[1]; sub(/.*\(/, "", f); sub(/\)/, "", f); if (f != $2) { ++dropped; next } } } { print } END { print dropped + 0 > "tmp_rep_dropped.txt" }' tmp_cl19.clumped > tmp_cl19_as107.clumped
awk -f $compare_tokens tmp_cl19_as107.clumped tmp_cl107.clumped 2> /dev/null && test "$(cat tmp_rep_dropped.txt)" -gt 0 && ! cmp -s tmp_cl19.clumped tmp_clthree.clumped
next_check $?
same_clump "--clump tmp_res1.txt,tmp_res2.txt --clump-index-first" clumped tmp_cltwo.clumped
# With two files, --clump-best requires --clump-index-first.  Proxies come
# from either file, unless --clump-replicate is also given; then they only
# come from the second file (and clumps aren't filtered).
same_clump "--clump tmp_res1.txt,tmp_res3.txt --clump-best --clump-index-first" "clumped clumped.best"
cp tmp_cl19.clumped.best tmp_clbest2.clumped.best
same_clump "--clump tmp_res1.txt,tmp_res3.txt --clump-best --clump-replicate --clump-index-first" "clumped clumped.best"
! cmp -s tmp_cl19.clumped.best tmp_clbest2.clumped.best && grep -q ' NA ' tmp_cl19.clumped.best
check "--clump-replicate restricts --clump-best proxies to the second file" $?

# --- Simulation ---

# Simulated genotypes depend on each version's random numbers, so only the
# sample IDs and variant names are compared.
printf '20 snpA 0.1 0.4 1.5 1.5\n' > tmp_siml.txt
p107 --simulate tmp_siml.txt --simulate-ncases 30 --simulate-ncontrols 30 --simulate-label XY --make-bed --out tmp_sl107
plink --simulate tmp_siml.txt --simulate-ncases 30 --simulate-ncontrols 30 --simulate-label XY --make-bed --out tmp_sl19 > /dev/null 2>&1
test -s tmp_sl19.fam && cmp -s tmp_sl19.fam tmp_sl107.fam && cmp -s <(cut -f 1-4 tmp_sl19.bim) <(cut -f 1-4 tmp_sl107.bim) && \
    test "$(cut -d ' ' -f 1 tmp_sl19.fam | head -1)" = XY-per0
next_check $?

test $check_idx -eq ${#checks[@]}
check "every PLINK 1.07 check ran" $?

exit $failed
