#!/bin/bash

# --gwas-ssf rsid=infer when the first rsID appears late in a large file.
#
# The rsID column is guessed from the second line; on encountering a later
# rsID, --gwas-ssf restarts with the column enabled.  The input stream wasn't
# rewound for the restart, which only matters once the file is larger than
# the stream's buffer.

set -exo pipefail

plink2="$1/plink2 $2 $3"

python3 - <<'PY'
n = 20000
with open("tmp_assoc.glm.linear", "w") as f:
    f.write("#CHROM\tPOS\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tOMITTED\tA1_FREQ\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\tERRCODE\n")
    for i in range(n):
        vid = "rs%d" % (i + 1) if i == n - 10 else "v%d" % (i + 1)
        f.write("1\t%d\t%s\tA\tG\tN\tG\tA\t0.3125\tADD\t40\t0.0512047\t0.263066\t0.194646\t0.846707\t.\n" % (i + 1, vid))
PY
$plink2 --gwas-ssf rsid=infer file=tmp_assoc.glm.linear --out tmp_ssf > /dev/null
out=tmp_assoc.glm.linear.ssf.tsv
[[ $(wc -l < $out) -eq 20001 ]]
head -n 1 $out | tr '\t' '\n' | grep -qx rsid
[[ $(grep -c $'\trs19991\t' $out) -eq 1 ]]
