# Compares a PLINK 1.9 .epi.qt report against a plink2 one.
#
# 1.9 reports the squared t-statistic as a 1-df chi-square, so the comparison
# is against the square root of it; and it prints six significant digits, so
# the tolerance is half of the last printed place.  BETA_INT is comparable
# directly: flipping a variant's allele coding negates the interaction
# coefficient, and both programs count the A1 allele of the .bim here.
function abs(x) { return (x < 0)? -x : x }
function key(a, b) { return (a < b)? (a "|" b) : (b "|" a) }
FNR == NR { if (FNR > 1) { k = key($2, $4); beta[k] = $5 + 0; stat[k] = $6 + 0; ++n1 }; next }
/^#/ { next }
{
    ++n2
    k = key($2, $4)
    if (!(k in beta)) { print "pair missing from the PLINK 1.9 report: " k; exit 1 }
    if (abs(beta[k] - ($5 + 0)) > 5.5e-6 * (1 + abs(beta[k]))) {
        print "BETA_INT differs on " k ": " beta[k] " vs " $5; exit 1
    }
    want_t = sqrt(stat[k])
    if (abs(want_t - abs($7 + 0)) > 5.5e-6 * (1 + want_t)) {
        print "T_STAT differs on " k ": sqrt(" stat[k] ") vs " $7; exit 1
    }
}
END {
    if (n1 != n2) { print "pair count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no pairs tested, so nothing was compared"; exit 1 }
    print n1 " pairs matched PLINK 1.9"
}
