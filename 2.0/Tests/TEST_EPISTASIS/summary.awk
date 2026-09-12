# Compares a PLINK 1.9 .epi.cc.summary against a plink2 one.  N_TOT is the
# count of pairs that were actually tested, so it is also what checks that both
# programs threw out the same degenerate pairs.
function abs(x) { return (x < 0)? -x : x }
FNR == NR { if (FNR > 1) { nsig[$2] = $3; ntot[$2] = $4; best[$2] = $6 + 0; ++n1 }; next }
/^#/ { next }
{
    ++n2
    if (!($2 in nsig)) { print "variant missing from the PLINK 1.9 summary: " $2; exit 1 }
    if (nsig[$2] != $3) { print "N_SIG differs on " $2 ": " nsig[$2] " vs " $3; exit 1 }
    if (ntot[$2] != $4) { print "N_TOT differs on " $2 ": " ntot[$2] " vs " $4; exit 1 }
    if (abs(best[$2] - $6) > 1e-3 * (1 + abs(best[$2]))) {
        print "BEST_CHISQ differs on " $2 ": " best[$2] " vs " $6; exit 1
    }
}
END {
    if (n1 != n2) { print "summary row count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "empty summary, so nothing was compared"; exit 1 }
    print n1 " summary rows matched"
}
