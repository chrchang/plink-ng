#!/bin/bash

# PLINK 1.9 with contig names longer than its filename buffer.
#
# Usage: ./run_tests.sh {plink 1.9 dir}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.

set -exo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}

long_contig=$(printf 'c%.0s' $(seq 1 15900))
awk -v OFS='\t' 'BEGIN { for (i = 1; i <= 40; ++i) { print "s" i, "s" i, 0, 0, 1 + i % 2, 1 } }' > tmp.fam
awk -v OFS='\t' -v c="$long_contig" 'BEGIN { for (i = 1; i <= 400; ++i) { print c, "rs" i, 0, i, "A", "G" } }' > tmp_long.bim
sed "s/^c*/chrZ/" tmp_long.bim > tmp_short.bim
python3 -c "
import random
random.seed(1)
open('tmp.bed', 'wb').write(bytes([0x6c, 0x1b, 1]) + bytes(random.randrange(256) for _ in range(4000)))
"
awk '{print $1, $2, "c" ($5)}' tmp.fam > tmp_clusters.txt

# --freq --within output must not depend on the contig name's length.
for f in long short; do
    plink --bed tmp.bed --bim tmp_$f.bim --fam tmp.fam --allow-extra-chr --within tmp_clusters.txt --freq --out freq_$f
    awk '{$1 = ""; print}' freq_$f.frq.strat > freq_$f.txt
done
cmp freq_long.txt freq_short.txt

# --recode must reject a name that does not fit, rather than overflow.
for mode in beagle fastphase HV; do
    set +e
    plink --bed tmp.bed --bim tmp_long.bim --fam tmp.fam --allow-extra-chr --allow-no-sex --recode $mode --out recode_$mode > tmp.log 2>&1
    rc=$?
    set -e
    [[ $rc -ne 0 && $rc -lt 128 ]] || exit 1
    grep -q 'too long for --recode' tmp.log
done
