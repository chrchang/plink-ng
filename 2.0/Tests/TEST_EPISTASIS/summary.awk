# Compares a PLINK 1.9 .epi.cc.summary against a plink2 one.  N_TOT is the
# count of pairs that were actually tested, so it is also what checks that both
# programs threw out the same degenerate pairs.
#
# awk still runs END after an exit from a main rule, so a failure sets failed
# and END stops at once instead of reporting on a file read only part way.
function abs(x) { return (x < 0)? -x : x }
function fail(msg) { print msg; failed = 1; exit 1 }
FNR == NR { if (FNR > 1) { nsig[$2] = $3; ntot[$2] = $4; best[$2] = $6 + 0; ++n1 }; next }
/^#/ { next }
{
    ++n2
    if (!($2 in nsig)) { fail("variant missing from the PLINK 1.9 summary: " $2) }
    if (nsig[$2] != $3) { fail("N_SIG differs on " $2 ": " nsig[$2] " vs " $3) }
    if (ntot[$2] != $4) { fail("N_TOT differs on " $2 ": " ntot[$2] " vs " $4) }
    if (abs(best[$2] - $6) > 1e-3 * (1 + abs(best[$2]))) {
        fail("BEST_CHISQ differs on " $2 ": " best[$2] " vs " $6)
    }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "summary row count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "empty summary, so nothing was compared"; exit 1 }
    print n1 " summary rows matched"
}
