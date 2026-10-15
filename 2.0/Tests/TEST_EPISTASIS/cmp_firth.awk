# Compares a covariate-adjusted .epi.cc against firth/oracle.tsv, whose
# FIRTH_PLR column (3) is the penalized likelihood ratio from direct
# maximization of the Jeffreys-penalized log-likelihood in R, and whose LR
# column (4) is the unpenalized likelihood ratio from R's glm().  The
# full-model penalized log-likelihood agrees with logistf's.  logistf's own
# constrained null fit is not used, since it stops short of the maximum on
# some of these sparse pairs.
#
#   -v col=3  every pair of the oracle must be in the report (Firth)
#   -v col=4  only the pairs the report has are compared (logistic, which
#             leaves out the pairs it cannot fit)
#
# STAT is compared with a relative tolerance, since plink2 prints six
# significant digits.
function abs(x) { return (x < 0)? -x : x }
function fail(msg) { print msg; failed = 1; exit 1 }
FNR == NR { if (FNR > 1 && NF) { want[$1 "|" $2] = $col + 0 }; next }
/^#/ { next }
{
    k = $2 "|" $4
    if (!(k in want)) { fail("pair missing from the oracle: " k) }
    seen[k] = 1
    if (abs(want[k] - ($5 + 0)) > 1e-4 * (1 + abs(want[k]))) {
        fail("STAT differs on " k ": " want[k] " vs " $5)
    }
    ++matched
}
END {
    if (failed) { exit 1 }
    if (col == 3) {
        for (k in want) { if (!(k in seen)) { print "pair missing from the plink2 report: " k; exit 1 } }
    }
    if (matched < 100) { print "only " matched + 0 " pairs compared"; exit 1 }
    print matched " pairs matched the oracle"
}
