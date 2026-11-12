#!/bin/bash

# PLINK 1.9 computations with no direct plink2 or PLINK 1.07 equivalent,
# checked against independent recalculations from PLINK's own matrix and
# frequency outputs: --rel-cutoff (from genotypes, --grm-bin and --grm-gz),
# --groupdist, --regress-distance, --regress-rel, --pca-cluster-names,
# --pca-clusters and --mendel-multigen.
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

# Doubles (f8) or floats (f4) from a binary matrix, one per line.
binvals() {
    od -A n -v -t "$1" "$2" | tr -s ' \t' '\n\n' | grep -v '^$'
}

# Drops the diagonal from a lower-triangular matrix that includes it.
drop_diag() {
    awk '{ if (j < i) print; if (j == i) { ++i; j = 0 } else ++j }'
}

# The rows of a .fam/.id file at the given 0-based indices, as "FID IID".
pick_ids() {
    awk 'NR == FNR { keep[$1] = 1; next } (FNR - 1) in keep { print $1, $2 }' "$1" "$2"
}

# 15 families of four plus 23 unrelated samples (83 in all) over 500
# variants; a copy with 2% missing calls; and a copy with quantitative
# phenotypes.
awk -v n_fam=15 -v n_unrel=23 -v n_var=500 -v seed=11 -f make_families.awk > tmp_fam.ped
awk 'BEGIN { for (j = 1; j <= 500; ++j) print 1 + int((j - 1) / 250), "v" j, 0, 1000 * j }' > tmp_fam.map
awk -v seed=5 '{
    for (i = 7; i < NF; i += 2) {
        seed = (seed * 69069 + 1) % 4294967296
        if (seed / 4294967296 < 0.02) $i = $(i + 1) = 0
    }
    print
}' tmp_fam.ped > tmp_fm.ped
awk -v n_fam=15 -v n_unrel=23 -v n_var=500 -v seed=11 -v qt=1 -f make_families.awk > tmp_fq.ped
cp tmp_fam.map tmp_fm.map
cp tmp_fam.map tmp_fq.map
for f in fam fm fq; do
    plink --file tmp_$f --make-bed --out tmp_$f > /dev/null 2>&1
done
n=$(wc -l < tmp_fam.fam | tr -d ' ')

# --- --rel-cutoff ---

# From genotypes, against the greedy pruning applied to --make-rel's matrix.
plink --bfile tmp_fam --make-rel triangle bin --out tmp_rel > /dev/null 2>&1
binvals f8 tmp_rel.rel.bin | drop_diag > tmp_rel_vals.txt
rc_ok=0
prev=
for cutoff in 0.025 0.1 0.2 0.35; do
    rm -f tmp_rc.rel.id
    plink --bfile tmp_fam --rel-cutoff $cutoff --out tmp_rc > /dev/null 2>&1
    awk -v n=$n -v cutoff=$cutoff -f rel_cutoff.awk tmp_rel_vals.txt > tmp_rc_idx.txt
    pick_ids tmp_rc_idx.txt tmp_fam.fam > tmp_rc_exp.txt
    test -s tmp_rc_exp.txt && cmp -s tmp_rc_exp.txt <(awk '{ print $1, $2 }' tmp_rc.rel.id) || rc_ok=1
    # Something must be pruned at the lower cutoffs.
    if [[ $cutoff != 0.35 ]]; then
        test "$(wc -l < tmp_rc_exp.txt)" -lt $n || rc_ok=1
    fi
    if [[ -n "$prev" ]]; then
        ! cmp -s "$prev" tmp_rc_exp.txt || [[ $cutoff == 0.35 ]] || rc_ok=1
    fi
    cp tmp_rc_exp.txt tmp_rc_prev.txt
    prev=tmp_rc_prev.txt
done
check "--rel-cutoff keeps the samples greedy pruning keeps" $rc_ok

# With --grm-bin and --grm-gz input, against the same pruning applied to the
# GRM's own values; the genotype path must agree.
plink --bfile tmp_fam --make-grm-bin --out tmp_grm > /dev/null 2>&1
plink --bfile tmp_fam --make-grm-gz --out tmp_grmgz > /dev/null 2>&1
binvals f4 tmp_grm.grm.bin | drop_diag > tmp_grm_vals.txt
gunzip -c tmp_grmgz.grm.gz | awk '$1 != $2 { print $4 }' > tmp_grmgz_vals.txt
grm_ok=0
for cutoff in 0.025 0.1; do
    rm -f tmp_rcb.grm.id tmp_rcg.grm.id tmp_rc.rel.id
    plink --grm-bin tmp_grm --rel-cutoff $cutoff --out tmp_rcb > /dev/null 2>&1
    plink --grm-gz tmp_grmgz --rel-cutoff $cutoff --out tmp_rcg > /dev/null 2>&1
    plink --bfile tmp_fam --rel-cutoff $cutoff --out tmp_rc > /dev/null 2>&1
    awk -v n=$n -v cutoff=$cutoff -f rel_cutoff.awk tmp_grm_vals.txt > tmp_rcb_idx.txt
    awk -v n=$n -v cutoff=$cutoff -f rel_cutoff.awk tmp_grmgz_vals.txt > tmp_rcg_idx.txt
    pick_ids tmp_rcb_idx.txt tmp_grm.grm.id > tmp_rcb_exp.txt
    pick_ids tmp_rcg_idx.txt tmp_grmgz.grm.id > tmp_rcg_exp.txt
    test -s tmp_rcb_exp.txt && test "$(wc -l < tmp_rcb_exp.txt)" -lt $n && \
        cmp -s tmp_rcb_exp.txt <(awk '{ print $1, $2 }' tmp_rcb.grm.id) && \
        cmp -s tmp_rcg_exp.txt <(awk '{ print $1, $2 }' tmp_rcg.grm.id) && \
        cmp -s tmp_rcb_exp.txt <(awk '{ print $1, $2 }' tmp_rc.rel.id) || grm_ok=1
done
check "--rel-cutoff on --grm-bin/--grm-gz input keeps the expected samples" $grm_ok

# --- --groupdist ---

# Mean, standard error ("sem") and median of the AA, AU and UU
# distances, against --distance's matrix; with and without missing calls.
# (The jackknife differences are random, so they aren't checked.)
gd_ok=0
for f in fam fm; do
    plink --bfile tmp_$f --distance triangle --out tmp_gd_$f > /dev/null 2>&1
    tr -s ' \t' '\n\n' < tmp_gd_$f.dist | grep -v '^$' > tmp_gd_vals.txt
    awk -f groupdist.awk tmp_$f.fam tmp_gd_vals.txt > tmp_gd_exp.txt
    plink --bfile tmp_$f --groupdist 100 --seed 1 --out tmp_gd_$f > /dev/null 2>&1
    awk -F'[:,()]' '/Mean \(sem\), median dists between/ {
        g = ($0 ~ /2x affected/)? "AA" : ($0 ~ /aff\. and unaff\./)? "AU" : "UU"
        print g, $5 + 0, $6 + 0, $8 + 0
    }' tmp_gd_$f.log > tmp_gd_got.txt
    awk -v tol=1e-5 -f "$compare_tokens" tmp_gd_exp.txt tmp_gd_got.txt || gd_ok=1
done
check "--groupdist means, standard errors and medians" $gd_ok

# --- --regress-distance and --regress-rel ---

# Both regression slopes, against --distance's and --make-rel's matrices.
plink --bfile tmp_fq --distance triangle bin --out tmp_rd > /dev/null 2>&1
binvals f8 tmp_rd.dist.bin > tmp_rd_vals.txt
awk -f pair_regress.awk tmp_fq.fam tmp_rd_vals.txt > tmp_rd_exp.txt
plink --bfile tmp_fq --regress-distance 100 --seed 1 --out tmp_rd > /dev/null 2>&1
awk -F': ' '/^Regression slope \(y = genomic distance/ { a = $2 } /^Regression slope \(y = avg phenotype/ { b = $2 } END { print a, b }' tmp_rd.log > tmp_rd_got.txt
awk -v tol=1e-5 -f "$compare_tokens" tmp_rd_exp.txt tmp_rd_got.txt
check "--regress-distance slopes" $?

plink --bfile tmp_fq --make-rel triangle bin --out tmp_rr > /dev/null 2>&1
binvals f8 tmp_rr.rel.bin > tmp_rr_vals.txt
awk -v diag=1 -f pair_regress.awk tmp_fq.fam tmp_rr_vals.txt > tmp_rr_exp.txt
plink --bfile tmp_fq --regress-rel 100 --seed 1 --out tmp_rr > /dev/null 2>&1
awk -F': ' '/^Regression slope \(y = genomic relationship/ { a = $2 } /^ *\(y = avg phenotype, x = genomic relationship\)/ { b = $2 } END { print a, b }' tmp_rr.log > tmp_rr_got.txt
awk -v tol=1e-5 -f "$compare_tokens" tmp_rr_exp.txt tmp_rr_got.txt
check "--regress-rel slopes" $?

# --- --pca-cluster-names and --pca-clusters ---

# Samples cycle through clusters cA, cC and cB; the PCA uses the cA and cB
# samples and projects the rest.  Cluster samples' eigenvectors must equal a
# --keep run's (with frequencies from all founders, as --pca uses), and every
# sample's must match its projection onto that run's variant weights.
awk '{ print $1, $2, substr("cAcCcB", 2 * ((NR - 1) % 3) + 1, 2) }' tmp_fam.fam > tmp_clus.txt
awk '$3 != "cC" { print $1, $2 }' tmp_clus.txt > tmp_clus_keep.txt
printf 'cA\ncB\n' > tmp_clus_names.txt
plink --bfile tmp_fam --freqx --out tmp_pcf > /dev/null 2>&1
plink --bfile tmp_fam --recode A --out tmp_pcf > /dev/null 2>&1
plink --bfile tmp_fam --keep tmp_clus_keep.txt --read-freq tmp_pcf.frqx --pca 4 var-wts --out tmp_pck > /dev/null 2>&1
plink --bfile tmp_fam --within tmp_clus.txt --pca 4 --pca-cluster-names cA cB --out tmp_pcn > /dev/null 2>&1
plink --bfile tmp_fam --within tmp_clus.txt --pca 4 --pca-clusters tmp_clus_names.txt --out tmp_pcc > /dev/null 2>&1
test -s tmp_pcn.eigenvec && test "$(wc -l < tmp_pcn.eigenvec)" -eq $n && \
    cmp -s tmp_pck.eigenval tmp_pcn.eigenval && \
    cmp -s tmp_pck.eigenvec <(awk 'NR == FNR { keep[$1, $2] = 1; next } ($1, $2) in keep' tmp_clus_keep.txt tmp_pcn.eigenvec)
check "--pca-cluster-names matches a PCA of the cluster samples" $?
awk -v k=4 -f pca_project.awk tmp_pcf.frqx tmp_pck.eigenvec.var tmp_pck.eigenvec tmp_pcf.raw > tmp_pc_proj.txt
awk -v abstol=2e-4 -f "$compare_tokens" tmp_pc_proj.txt tmp_pcn.eigenvec
check "--pca-cluster-names projects other samples onto the cluster PCs" $?
test -s tmp_pcc.eigenvec && cmp -s tmp_pcc.eigenvec tmp_pcn.eigenvec && cmp -s tmp_pcc.eigenval tmp_pcn.eigenval
check "--pca-clusters matches --pca-cluster-names" $?

# --- --mendel-multigen ---

# Child ch is C/C at every variant, with an A/C mother and a father, fa, who
# is missing at v1-v3 and A/A at v4.  Without --mendel-multigen only v4 is an
# error.  With it, fa's missing calls are inferred from his ancestors: at v1
# his parents are both A/A; at v3 his father gf is also missing, but gf's
# parents are A/A; at v2 his mother gm is A/C, so there's no error.
plink --file mendel_multigen --mendel --out tmp_mg0 > /dev/null 2>&1
plink --file mendel_multigen --mendel --mendel-multigen --out tmp_mg1 > /dev/null 2>&1
mendel_rows() {
    awk 'NR > 1 { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }' "$1"
}
cmp -s <(mendel_rows tmp_mg0.mendel) <(echo 'F ch 1 v4 3 A/A x */* -> C/C') && \
    cmp -s <(mendel_rows tmp_mg1.mendel) <(printf 'F ch 1 v%s 3 A/A x */* -> C/C\n' 1 3 4)
check "--mendel-multigen infers missing parental calls from ancestors" $?

exit $failed
