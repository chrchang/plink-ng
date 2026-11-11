#!/bin/bash

# PLINK 1.9 commands with no plink2 equivalent, checked without a reference
# binary: --tdt (against transmissions recounted by check_tdt.awk),
# --within/--filter-females/--write-cluster, and --fast-epistasis --parallel
# with --epistasis-summary-merge (against a single run).  (Ported from the
# old 1.9/tests_old suite, which compared against PLINK 1.07.)
#
# Usage: ./run_tests.sh
# Requires plink (1.9) to be in the system PATH.

set -uo pipefail

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}

# --tdt: 121 trios with affected children, plus one unrelated sample (364
# samples), with missing calls and some Mendel errors.
awk -v n_trio=121 -v n_var=200 -v seed=19 -f make_trios.awk > tmp_trio.ped
awk 'BEGIN { for (j = 1; j <= 200; ++j) print 1, "v" j, 0, 1000 * j }' > tmp_trio.map
plink --file tmp_trio --tdt --out tmp_tdt > /dev/null
awk -f check_tdt.awk tmp_tdt.tdt tmp_trio.ped
check "--tdt T/U counts match the transmissions" $?
plink --file tmp_trio --make-bed --out tmp_trio > /dev/null
plink --bfile tmp_trio --tdt --out tmp_tdt_b > /dev/null
cmp -s tmp_tdt.tdt tmp_tdt_b.tdt
check "--tdt gives the same report from a .bed" $?

# --within --filter-females --write-cluster.  Females are the odd-numbered
# samples; every fifth sample has no cluster.
plink --dummy 201 300 0.02 acgt --seed 17 --make-bed --out tmp_epi > /dev/null
awk '{ $5 = 1 + NR % 2; print }' tmp_epi.fam > tmp_epi.fam.new
mv tmp_epi.fam.new tmp_epi.fam
awk 'NR % 5 != 0 { print $1, $2, "c" (NR % 4) }' tmp_epi.fam > tmp_within.txt
plink --bfile tmp_epi --within tmp_within.txt --filter-females --write-cluster --out tmp_wc > /dev/null
plink --bfile tmp_epi --within tmp_within.txt --filter-females --write-cluster omit-unassigned --out tmp_wco > /dev/null
awk 'NR == FNR { c[$1 " " $2] = $3; next } $5 == 2 { k = $1 " " $2; print k, (k in c)? c[k] : "NA" }' tmp_within.txt tmp_epi.fam > tmp_wc_expected.txt
cmp -s tmp_wc.clst tmp_wc_expected.txt && cmp -s tmp_wco.clst <(grep -v ' NA$' tmp_wc_expected.txt)
check "--within --filter-females --write-cluster" $?

# --fast-epistasis --parallel 3 ways, then --epistasis-summary-merge, must
# reproduce a single run.  Only the first piece has a header line.
plink --bfile tmp_epi --fast-epistasis no-ueki --epi1 0.01 --out tmp_epi1 > /dev/null
for i in 1 2 3; do
    plink --bfile tmp_epi --fast-epistasis no-ueki --epi1 0.01 --parallel $i 3 --out tmp_epi3 > /dev/null
done
cat tmp_epi3.epi.cc.1 tmp_epi3.epi.cc.2 tmp_epi3.epi.cc.3 > tmp_epi3.epi.cc
plink --epistasis-summary-merge tmp_epi3.epi.cc 3 --out tmp_epi3m > /dev/null
cmp -s tmp_epi1.epi.cc tmp_epi3.epi.cc && cmp -s tmp_epi1.epi.cc.summary tmp_epi3m.epi.cc.summary && test "$(wc -l < tmp_epi1.epi.cc)" -gt 100
check "--fast-epistasis --parallel + --epistasis-summary-merge matches a single run" $?

exit $failed
