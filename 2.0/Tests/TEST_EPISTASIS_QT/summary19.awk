# Compares a PLINK 1.9 .epi.qt.summary against a plink2 one.
#
# N_SIG is not compared: 1.9 counts a pair as significant by the normal
# approximation to the t-statistic's p-value, and this reports the t p-value
# itself, so the two disagree on pairs between the thresholds.  N_TOT (the
# pairs that were fit at all) and BEST_CHISQ (the largest squared
# t-statistic, which 1.9 prints to four significant digits) are directly
# comparable.
function abs(x) { return (x < 0)? -x : x }
FNR == NR { if (FNR > 1) { ntot[$2] = $4; best[$2] = $6 + 0; ++n1 }; next }
/^#/ { next }
{
    ++n2
    if (!($2 in ntot)) { print "variant missing from the PLINK 1.9 summary: " $2; exit 1 }
    if (ntot[$2] != $4) { print "N_TOT differs on " $2 ": " ntot[$2] " vs " $4; exit 1 }
    if (abs(best[$2] - $6) > 5.5e-4 * (1 + abs(best[$2]))) {
        print "BEST_CHISQ differs on " $2 ": " best[$2] " vs " $6; exit 1
    }
}
END {
    if (n1 != n2) { print "summary row count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "empty summary, so nothing was compared"; exit 1 }
    print n1 " summary rows matched PLINK 1.9"
}
