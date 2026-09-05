# Compares the .flipscan.verbose files.  Both list one line per
# above-threshold neighbor pair of a flagged variant; the allele columns use
# PLINK 1.9's A1 and plink2's ALT.  Recoding one variant of a pair flips the
# sign of both correlations, so what is comparable is their magnitudes and
# whether they agree in sign, which is the quantity the scan is about.
function abs(x) { return (x < 0)? -x : x }
function same_value(a, b) { return abs(a - b) <= 1e-3 + 1e-3 * abs(a) }
function sgn(x) { return (x > 0)? 1 : ((x < 0)? -1 : 0) }
FNR == NR {
    if (FNR > 1) { key = $2 "|" $5; rcase[key] = $8; rctrl[key] = $9; ++n1 }
    next
}
/^#/ { next }
{
    ++n2;
    key = $2 "|" $5;
    if (!(key in rcase)) { print "pair missing from the PLINK 1.9 report: " key; failed = 1; exit 1 }
    if (!same_value(abs(rcase[key]), abs($8)) || !same_value(abs(rctrl[key]), abs($9)) ||
        (sgn(rcase[key]) * sgn(rctrl[key]) != sgn($8) * sgn($9))) {
        print "correlation differs on " key ": " rcase[key] "/" rctrl[key] " vs " $8 "/" $9; failed = 1; exit 1
    }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "pair count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no pairs compared"; exit 1 }
    print n1 " neighbor pairs matched"
}
