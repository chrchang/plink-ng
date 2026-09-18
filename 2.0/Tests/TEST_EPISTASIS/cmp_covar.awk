# Compares oracle_covar.py's output against a plink2 .epi.cc from a
# covariate-adjusted run.  DF has to agree exactly; STAT is compared with a
# relative tolerance, since plink2 prints six significant digits.
#
# awk still runs END after an exit from a main rule, so a failure sets failed
# and END stops at once instead of reporting on a file read only part way.
function abs(x) { return (x < 0)? -x : x }
function fail(msg) { print msg; failed = 1; exit 1 }
FNR == NR { if (FNR > 1 && NF) { stat[$1 "|" $2] = $3 + 0; df[$1 "|" $2] = $4 }; next }
/^#/ { next }
{
    k = $2 "|" $4
    if (!(k in stat)) { fail("pair missing from the oracle: " k) }
    seen[k] = 1
    if (df[k] != $6) { fail("DF differs on " k ": " df[k] " vs " $6) }
    if (abs(stat[k] - ($5 + 0)) > 1e-4 * (1 + abs(stat[k]))) {
        fail("STAT differs on " k ": " stat[k] " vs " $5)
    }
    ++matched
}
END {
    if (failed) { exit 1 }
    for (k in stat) { if (!(k in seen)) { print "pair missing from the plink2 report: " k; exit 1 } }
    if (matched == 0) { print "no pairs compared"; exit 1 }
    print matched " covariate-adjusted pairs matched the oracle"
}
