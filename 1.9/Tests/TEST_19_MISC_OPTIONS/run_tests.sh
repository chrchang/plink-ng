#!/bin/bash

# PLINK 1.9 options with no plink2 or PLINK 1.07 equivalent, checked against
# expected results computed here, or against runs without the option:
# --meta-analysis-*-field, --meta-analysis-report-dups,
# --gene-report-snp-field, --missing-genotype2, --allow-no-vars,
# --list-duplicate-vars, --make-perm-pheno, --pool-size, --homozyg-match,
# --ld-xchr and --read-dists (with --groupdist, --regress-distance,
# --ibs-test, --cluster, --neighbour and --mds-plot).
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

# --- Meta-analysis and gene-report field names ---

# Two association reports over 40 variants.
for k in 1 2; do
    awk -v k=$k 'BEGIN {
        s = 7 * k
        print "CHR SNP BP A1 A2 OR SE P NMISS"
        for (j = 1; j <= 40; ++j) {
            s = (s * 69069 + 1) % 4294967296
            u = int(s / 65536) % 1000 / 1000
            print 1, "v" j, 1000 * j, "A", "C", exp(u - 0.5), 0.1 + u / 10, u, 300 + j
        }
    }' > tmp_meta$k.txt
    awk 'NR == 1 { print "CC SS BB XA XB OR SEX PV NN"; next } { print }' tmp_meta$k.txt > tmp_meta_renamed$k.txt
done
# All eight --meta-analysis-*-field flags, on copies with every such column
# renamed; the report must match the default one, and the renamed copies
# must fail without the flags.
plink --meta-analysis tmp_meta1.txt tmp_meta2.txt + weighted-z --out tmp_ma_def > /dev/null 2>&1
plink --meta-analysis tmp_meta_renamed1.txt tmp_meta_renamed2.txt + weighted-z \
    --meta-analysis-chr-field CC --meta-analysis-snp-field SS --meta-analysis-bp-field BB \
    --meta-analysis-a1-field XA --meta-analysis-a2-field XB --meta-analysis-se-field SEX \
    --meta-analysis-p-field PV --meta-analysis-ess-field NN --out tmp_ma_ren > /dev/null 2>&1
plink --meta-analysis tmp_meta_renamed1.txt tmp_meta_renamed2.txt + weighted-z --out tmp_ma_bad > /dev/null 2>&1
bad_rc=$?
test -s tmp_ma_def.meta && cmp -s tmp_ma_def.meta tmp_ma_ren.meta && test $bad_rc -ne 0 && \
    grep -q 'WEIGHTED_Z' tmp_ma_def.meta && test "$(wc -l < tmp_ma_def.meta)" -eq 41
check "--meta-analysis-*-field flags match the default field names" $?

# --meta-analysis-report-dups lists the variants a file repeats (in .prob).
{ cat tmp_meta1.txt; tail -3 tmp_meta1.txt; } > tmp_meta_dup.txt
rm -f tmp_md_on.prob tmp_md_off.prob
plink --meta-analysis tmp_meta_dup.txt tmp_meta2.txt --meta-analysis-report-dups --out tmp_md_on > /dev/null 2>&1
plink --meta-analysis tmp_meta_dup.txt tmp_meta2.txt --out tmp_md_off > /dev/null 2>&1
cmp -s tmp_md_on.prob <(printf 'tmp_meta_dup.txt\tv38\tDUPLICATE\ntmp_meta_dup.txt\tv39\tDUPLICATE\ntmp_meta_dup.txt\tv40\tDUPLICATE\n') && test ! -e tmp_md_off.prob
check "--meta-analysis-report-dups lists repeated variants" $?

# --gene-report-snp-field (only relevant with --extract): the report from a
# copy with the SNP column renamed must match, apart from the header name.
printf '1 2500 9500 GA\n1 15000 30000 GB\n' > tmp_genes.txt
awk 'NR == 1 || NR % 3 != 0' tmp_meta1.txt > tmp_gr.txt
awk 'NR == 1 { $2 = "MARKER" } { print }' tmp_gr.txt > tmp_gr_renamed.txt
awk 'NR > 1 && NR % 2 == 0 { print $2 }' tmp_meta1.txt > tmp_gr_extract.txt
plink --gene-report tmp_gr.txt tmp_genes.txt --extract tmp_gr_extract.txt --out tmp_gr_def > /dev/null 2>&1
plink --gene-report tmp_gr_renamed.txt tmp_genes.txt --extract tmp_gr_extract.txt --gene-report-snp-field MARKER --out tmp_gr_ren > /dev/null 2>&1
plink --gene-report tmp_gr_renamed.txt tmp_genes.txt --extract tmp_gr_extract.txt --out tmp_gr_bad > /dev/null 2>&1
bad_rc=$?
test -s tmp_gr_def.range.report && cmp -s tmp_gr_def.range.report <(sed 's/ MARKER / SNP /' tmp_gr_ren.range.report) && test $bad_rc -ne 0
check "--gene-report-snp-field matches the default field name" $?

# --- Input handling ---

plink --dummy 21 12 0.05 acgt --seed 3 --make-bed --out tmp_small > /dev/null 2>&1
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_small.bim > tmp_small.bim.new
mv tmp_small.bim.new tmp_small.bim

# --missing-genotype2: an N allele code with the flag must be read exactly
# like the default '.' code.
for code in dot N; do
    awk -v c=$code 'BEGIN { OFS = "\t" } NR % 4 == 0 { $6 = (c == "dot")? "." : c } { print }' tmp_small.bim > tmp_mg_$code.bim
    cp tmp_small.bed tmp_mg_$code.bed
    cp tmp_small.fam tmp_mg_$code.fam
done
plink --bfile tmp_mg_dot --make-bed --out tmp_mg_out_dot > /dev/null 2>&1
plink --bfile tmp_mg_N --missing-genotype2 N --make-bed --out tmp_mg_out_N > /dev/null 2>&1
plink --bfile tmp_mg_N --make-bed --out tmp_mg_out_plain > /dev/null 2>&1
test -s tmp_mg_out_dot.bim && cmp -s tmp_mg_out_dot.bim tmp_mg_out_N.bim && cmp -s tmp_mg_out_dot.bed tmp_mg_out_N.bed && \
    ! cmp -s tmp_mg_out_dot.bim tmp_mg_out_plain.bim
check "--missing-genotype2 sets an extra missing allele code" $?

# --allow-no-vars: a fileset without variants is rejected without the flag,
# and with it, --make-bed writes an empty variant set and the same samples.
: > tmp_nv.bim
cp tmp_small.fam tmp_nv.fam
printf '\154\033\001' > tmp_nv.bed
plink --bfile tmp_nv --freq --out tmp_nv_bad > /dev/null 2>&1
bad_rc=$?
rm -f tmp_nv_out.*
plink --bfile tmp_nv --allow-no-vars --make-bed --out tmp_nv_out > /dev/null 2>&1
test $bad_rc -ne 0 && test -e tmp_nv_out.bim && test ! -s tmp_nv_out.bim && test "$(wc -c < tmp_nv_out.bed)" -eq 3 && \
    cmp -s tmp_nv_out.fam tmp_small.fam
check "--allow-no-vars accepts a fileset without variants" $?

# --list-duplicate-vars, with each modifier, against dupvar.awk.  snp1/snp2
# and snp5-snp7 share positions and allele sets; snp5's alleles are swapped
# relative to snp6/snp7, so it drops out with 'require-same-ref'; snp9/snp10
# share a position but not alleles.
awk 'BEGIN { OFS = "\t" } {
    if (NR == 2 || NR == 3) { $4 = 2000; $5 = "X"; $6 = "Y" }
    if (NR >= 6 && NR <= 8) { $4 = 6000; $5 = (NR == 6)? "T" : "G"; $6 = (NR == 6)? "G" : "T" }
    if (NR == 10 || NR == 11) { $4 = 10000; $5 = (NR == 10)? "C" : "A"; $6 = (NR == 10)? "T" : "C" }
    print
}' tmp_small.bim > tmp_dv.bim
cp tmp_small.bed tmp_dv.bed
cp tmp_small.fam tmp_dv.fam
dv_ok=0
for mode in "" require-same-ref ids-only suppress-first "ids-only suppress-first"; do
    vars=""
    for m in $mode; do
        case $m in
            require-same-ref) vars="$vars -v same_ref=1" ;;
            ids-only) vars="$vars -v ids_only=1" ;;
            suppress-first) vars="$vars -v suppress_first=1" ;;
        esac
    done
    rm -f tmp_dv_out.dupvar
    plink --bfile tmp_dv --list-duplicate-vars $mode --out tmp_dv_out > /dev/null 2>&1
    test -s tmp_dv_out.dupvar && awk $vars -f dupvar.awk tmp_dv.bim | cmp -s - tmp_dv_out.dupvar || dv_ok=1
done
check "--list-duplicate-vars (each modifier) matches the expected groups" $dv_ok

# --make-perm-pheno: each column must keep missing phenotypes missing and be a
# permutation of the other phenotypes (within each cluster, with --within),
# and the permutations must not all be the identity.
plink --dummy 31 10 0.02 acgt --seed 3 --make-bed --out tmp_pp > /dev/null 2>&1
awk '{ $6 = (NR % 7 == 0)? -9 : 1 + (NR % 3 == 0); print }' tmp_pp.fam > tmp_pp.fam.new
mv tmp_pp.fam.new tmp_pp.fam
awk '{ print $1, $2, "c" (NR % 2) }' tmp_pp.fam > tmp_pp_within.txt
# $1: .pphe file, $2: cluster file (or empty).
perm_ok() {
    awk -v clfile="$2" 'BEGIN { if (clfile != "") while ((getline line < clfile) > 0) { split(line, f, " "); cl[f[1] " " f[2]] = f[3] } }
      NR == FNR { orig[$1 " " $2] = $6; next }
      {
        k = $1 " " $2
        c = (k in cl)? cl[k] : "all"
        if (orig[k] == -9) {
            for (i = 3; i <= NF; ++i) if ($i != -9) bad = 1
            next
        }
        want[c, orig[k]]++
        for (i = 3; i <= NF; ++i) {
            got[c, $i, i]++
            if ($i != orig[k]) moved[i] = 1
            if (i > maxcol) maxcol = i
        }
      }
      END {
        for (key in want) {
            split(key, kk, SUBSEP)
            for (i = 3; i <= maxcol; ++i) if (got[kk[1], kk[2], i] != want[key]) bad = 1
        }
        for (i = 3; i <= maxcol; ++i) if (!moved[i]) bad = 1
        exit (bad || maxcol < 4)
      }' tmp_pp.fam "$1"
}
plink --bfile tmp_pp --make-perm-pheno 5 --seed 2 --out tmp_pp_all > /dev/null 2>&1
plink --bfile tmp_pp --within tmp_pp_within.txt --make-perm-pheno 4 --seed 2 --out tmp_pp_cl > /dev/null 2>&1
test "$(wc -l < tmp_pp_all.pphe)" -eq 31 && perm_ok tmp_pp_all.pphe "" && perm_ok tmp_pp_cl.pphe tmp_pp_within.txt
check "--make-perm-pheno writes permutations of the phenotypes" $?

# --- Homozygosity pools ---

# 41 samples, 400 variants 5 kb apart.  Samples 1-12 share a homozygous
# haplotype over variants 80-260, samples 13-20 another over 150-330 (with
# allele differences in samples 17-20), and samples 21-23 a third over
# 335-400, so --homozyg group forms a pool of 20 and a pool of 3.
awk 'BEGIN {
    s = 77; n = 41; m = 400
    for (j = 1; j <= m; ++j) {
        s = (s * 69069 + 1) % 4294967296; shared1[j] = int(s / 65536) % 2
        s = (s * 69069 + 1) % 4294967296; shared2[j] = int(s / 65536) % 2
        s = (s * 69069 + 1) % 4294967296; shared3[j] = int(s / 65536) % 2
    }
    for (i = 1; i <= n; ++i) {
        line = "f" i " i" i " 0 0 " (1 + i % 2) " " (1 + (i % 3 == 0))
        for (j = 1; j <= m; ++j) {
            if (i <= 12 && j >= 80 && j <= 260) a = b = shared1[j]
            else if (i > 12 && i <= 20 && j >= 150 && j <= 330) a = b = (i <= 16 || j % 25)? shared2[j] : 1 - shared2[j]
            else if (i > 20 && i <= 23 && j >= 335) a = b = shared3[j]
            else {
                s = (s * 69069 + 1) % 4294967296; a = int(s / 65536) % 2
                s = (s * 69069 + 1) % 4294967296; b = int(s / 65536) % 2
            }
            line = line " " (a? "C" : "A") " " (b? "C" : "A")
        }
        print line
    }
}' > tmp_hz.ped
awk 'BEGIN { for (j = 1; j <= 400; ++j) print 1, "v" j, 0, 5000 * j }' > tmp_hz.map
plink --file tmp_hz --make-bed --out tmp_hz > /dev/null 2>&1
hz="--homozyg-snp 50 --homozyg-kb 200 --homozyg-window-snp 50 --homozyg-window-het 1"
plink --bfile tmp_hz $hz --homozyg group --out tmp_hz_def > /dev/null 2>&1
# Pool sizes and allelic-match group count of a .hom.overlap report.
pool_sizes() {
    awk 'NR > 1 && $1 ~ /^S/ && $2 != "CON" && $2 != "UNION" { n[$1]++ } END { for (k in n) print k, n[k] }' "$1" | sort
}
match_groups() {
    # (counted on insertion: mawk 1.3.3, Debian 10's default awk, has no
    # length(array))
    awk 'NR > 1 && $1 ~ /^S/ && $2 != "CON" && $2 != "UNION" && !(($1 " " $13) in g) { g[$1 " " $13] = 1; ++ct } END { print ct + 0 }' "$1"
}
# --pool-size 5 must give the default report without the 3-sample pool
# (ignoring blank separator lines).
plink --bfile tmp_hz $hz --homozyg group --pool-size 5 --out tmp_hz_ps > /dev/null 2>&1
test "$(pool_sizes tmp_hz_def.hom.overlap | tr '\n' ' ')" = "S1 20 S2 3 " && \
    cmp -s <(awk 'NR == FNR { if (FNR > 1 && $1 ~ /^S/ && $2 != "CON" && $2 != "UNION") n[$1]++; next } (FNR == 1 || !($1 ~ /^S/) || n[$1] >= 5) && NF' tmp_hz_def.hom.overlap tmp_hz_def.hom.overlap) \
      <(awk 'NF' tmp_hz_ps.hom.overlap)
check "--pool-size drops smaller --homozyg group pools" $?
# A stricter --homozyg-match splits pools into more allelic-match groups.
plink --bfile tmp_hz $hz --homozyg group --homozyg-match 0.99 --out tmp_hz_hm > /dev/null 2>&1
test "$(match_groups tmp_hz_hm.hom.overlap)" -gt "$(match_groups tmp_hz_def.hom.overlap)" && \
    test "$(pool_sizes tmp_hz_hm.hom.overlap)" = "$(pool_sizes tmp_hz_def.hom.overlap)"
check "--homozyg-match changes the allelic-match groups" $?

# --- --ld-xchr ---

# 81 samples (alternating sex), 30 chrX variants; males' heterozygous calls
# are set to missing.  --r2 under modes 1 and 2 must match xchr_r2.awk.
plink --dummy 81 30 0.03 acgt --seed 61 --make-bed --out tmp_xraw > /dev/null 2>&1
awk 'BEGIN { OFS = "\t" } { $1 = 23; $4 = 1000 * NR; print }' tmp_xraw.bim > tmp_xraw.bim.new
mv tmp_xraw.bim.new tmp_xraw.bim
awk '{ $5 = 1 + NR % 2; print }' tmp_xraw.fam > tmp_xraw.fam.new
mv tmp_xraw.fam.new tmp_xraw.fam
plink --bfile tmp_xraw --set-hh-missing --keep-allele-order --make-bed --out tmp_x > /dev/null 2>&1
plink --bfile tmp_x --recode A --out tmp_x > /dev/null 2>&1
xc_ok=0
for mode in 1 2; do
    plink --bfile tmp_x --r2 --ld-xchr $mode --ld-window-r2 0 --out tmp_x_r2_$mode > /dev/null 2>&1
    awk -v c1=3,6,7 -v c2=1,2,3 -v skip2=0 -f $compare_cols tmp_x_r2_$mode.ld <(awk -v mode=$mode -f xchr_r2.awk tmp_x.raw) || xc_ok=1
done
! cmp -s tmp_x_r2_1.ld tmp_x_r2_2.ld || xc_ok=1
check "--ld-xchr 1/2 --r2 match the recomputed values" $xc_ok

# Mode 3 (--indep[-pairwise] only) also gives females double weight.  For the
# adjacent pair whose weighted and unweighted (mode 2) r^2 differ most,
# --indep-pairwise 2 1 <t>, with t halfway between them, must prune one
# variant in exactly one of modes 2 and 3, as the r^2 values predict.
read -r xa xb r2_mode2 r2_mode3 <<< "$(awk 'NR == 1 {
    for (i = 7; i <= NF; ++i) { n = $i; sub(/_[^_]*$/, "", n); id[i - 6] = n }
    nv = NF - 6
    next
}
{ ++ns; male[ns] = ($5 == 1); for (j = 1; j <= nv; ++j) g[ns, j] = $(j + 6) }
END {
    for (a = 1; a < nv; ++a) {
        for (w = 0; w <= 1; ++w) {
            n = sx = sy = sxx = syy = sxy = 0
            for (i = 1; i <= ns; ++i) {
                x = g[i, a]; y = g[i, a + 1]
                if (x == "NA" || y == "NA") continue
                wt = (w && !male[i])? 2 : 1
                n += wt; sx += wt * x; sy += wt * y; sxx += wt * x * x; syy += wt * y * y; sxy += wt * x * y
            }
            cov = sxy - sx * sy / n
            r[w] = cov * cov / ((sxx - sx * sx / n) * (syy - sy * sy / n))
        }
        d = r[1] - r[0]
        if (d < 0) d = -d
        if (d > best) { best = d; line = id[a] " " id[a + 1] " " r[0] " " r[1] }
    }
    print line
}' tmp_x.raw)"
thresh=$(awk -v a="$r2_mode2" -v b="$r2_mode3" 'BEGIN { print (a + b) / 2 }')
ld3_ok=0
for mode in 2 3; do
    rm -f tmp_x_ip$mode.prune.out
    plink --bfile tmp_x --snps "$xa-$xb" --ld-xchr $mode --indep-pairwise 2 1 "$thresh" --out tmp_x_ip$mode > /dev/null 2>&1
    want=$(awk -v r=$([[ $mode == 2 ]] && echo "$r2_mode2" || echo "$r2_mode3") -v t="$thresh" 'BEGIN { print (r > t)? 1 : 0 }')
    test "$(grep -c . tmp_x_ip$mode.prune.out)" -eq "$want" || ld3_ok=1
done
! cmp -s tmp_x_ip2.prune.out tmp_x_ip3.prune.out || ld3_ok=1
check "--ld-xchr 3 weights females in --indep-pairwise" $ld3_ok

# --- --read-dists ---

# --read-dists loads the allele-count distances "--distance triangle bin"
# writes.  With a matrix written from the same data, every command that
# accepts it must give the same report as calculating from the genotypes.
# With missing calls, that takes the matching missingness treatment:
# --groupdist, --regress-distance and --ibs-test use --distance's default
# weighted correction, while --cluster, --neighbour and --mds-plot use none,
# i.e. a 'flat-missing' matrix.  (--cluster, --neighbour and --mds-plot
# convert the distances to IBS similarities; before 11 Nov 2026, they used
# the distances directly.)
# Report lines of a log (leaving out the "Options in effect" block, the run
# details, and the lines that differ between loading and calculating).
report_lines() {
    grep -vE '^(PLINK|Logging|Options|  --|Hostname|Working directory|Start time|End time|Random number seed|[0-9]+ MB RAM)' "$1" | grep -v -e '--read-dists: ' -e 'Distance matrix calculation complete' | grep .
}
rd_ok=0
cl_ok=0
nb_ok=0
for missing in 0 0.03; do
    plink --dummy 41 400 $missing acgt --seed 5 --make-bed --out tmp_rd > /dev/null 2>&1
    awk '{ $6 = 1 + NR % 2; print }' tmp_rd.fam > tmp_rd.fam.new
    mv tmp_rd.fam.new tmp_rd.fam
    awk 'BEGIN { s = 3 } { s = (s * 69069 + 1) % 4294967296; print $1, $2, int(s / 65536) % 1000 / 100 }' tmp_rd.fam > tmp_rd_qt.txt
    rm -f tmp_rd_dist.* tmp_rd_flat.*
    plink --bfile tmp_rd --distance triangle bin --out tmp_rd_dist > /dev/null 2>&1
    plink --bfile tmp_rd --distance triangle bin flat-missing --out tmp_rd_flat > /dev/null 2>&1
    # --groupdist, --regress-distance, --ibs-test (report in the log).
    for cmd in "--groupdist 50:dists between" "--pheno tmp_rd_qt.txt --regress-distance 50:Regression slope" "--ibs-test 2000:Between-group IBS"; do
        plink --bfile tmp_rd ${cmd%%:*} --seed 3 --out tmp_rd_geno > /dev/null 2>&1
        plink --bfile tmp_rd ${cmd%%:*} --seed 3 --read-dists tmp_rd_dist.dist.bin tmp_rd_dist.dist.id --out tmp_rd_read > /dev/null 2>&1
        grep -q "${cmd##*:}" tmp_rd_geno.log && cmp -s <(report_lines tmp_rd_geno.log) <(report_lines tmp_rd_read.log) || rd_ok=1
    done
    # --cluster (each linkage mode; also without the ID file).
    for cmd in "--cluster" "--cluster old-tiebreaks" "--cluster group-avg" "--cluster old-tiebreaks --K 3"; do
        rm -f tmp_rdc_geno.cluster* tmp_rdc_read.cluster*
        plink --bfile tmp_rd $cmd --out tmp_rdc_geno > /dev/null 2>&1
        plink --bfile tmp_rd $cmd --read-dists tmp_rd_flat.dist.bin tmp_rd_flat.dist.id --out tmp_rdc_read > /dev/null 2>&1
        for f in 1 2; do
            test -s tmp_rdc_geno.cluster$f && cmp -s tmp_rdc_geno.cluster$f tmp_rdc_read.cluster$f || cl_ok=1
        done
    done
    rm -f tmp_rdc_read.cluster1
    plink --bfile tmp_rd --cluster old-tiebreaks --read-dists tmp_rd_flat.dist.bin --out tmp_rdc_read > /dev/null 2>&1
    plink --bfile tmp_rd --cluster old-tiebreaks --out tmp_rdc_geno > /dev/null 2>&1
    cmp -s tmp_rdc_geno.cluster1 tmp_rdc_read.cluster1 || cl_ok=1
    # --neighbour and --mds-plot.  With missing calls, the default matrix
    # must give a different --neighbour report, so the flat-missing matrix
    # is the one that matters.
    for cmd in "--neighbour 1 3" "--cluster old-tiebreaks --mds-plot 2" "--cluster old-tiebreaks --K 3 --mds-plot 2 by-cluster"; do
        rm -f tmp_rdn_geno.* tmp_rdn_read.*
        plink --bfile tmp_rd $cmd --out tmp_rdn_geno > /dev/null 2>&1
        plink --bfile tmp_rd $cmd --read-dists tmp_rd_flat.dist.bin tmp_rd_flat.dist.id --out tmp_rdn_read > /dev/null 2>&1
        for e in nearest mds; do
            if [[ -e tmp_rdn_geno.$e ]]; then
                test -s tmp_rdn_geno.$e && cmp -s tmp_rdn_geno.$e tmp_rdn_read.$e || nb_ok=1
            fi
        done
    done
    if [[ $missing != 0 ]]; then
        plink --bfile tmp_rd --neighbour 1 3 --read-dists tmp_rd_dist.dist.bin tmp_rd_dist.dist.id --out tmp_rdn_dflt > /dev/null 2>&1
        ! cmp -s tmp_rdn_geno.nearest tmp_rdn_dflt.nearest || nb_ok=1
    fi
done
check "--read-dists gives the same --groupdist/--regress-distance/--ibs-test reports" $rd_ok
check "--read-dists gives the same --cluster solutions" $cl_ok
check "--read-dists gives the same --neighbour/--mds-plot reports" $nb_ok

exit $failed
