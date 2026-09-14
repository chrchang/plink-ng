# Compares a PLINK 1.9 .epi.cc report against a plink2 one.
#
# The pair key is normalized, since the two programs emit a pair in whatever
# order their scans reach it.  PLINK 1.9 prints STAT to four significant
# digits, so the comparison allows half of its last printed place rather than a
# fixed tolerance.
#
# A statistic that is mathematically zero is the exception: it comes back as
# rounding noise, and the noise is absolute rather than relative, since the
# statistic is a sum over 2N cell terms.  Such a row is accepted when both
# sides are within zero_tol of zero, and a near-zero row present on one side
# alone is not a disagreement either: under a permissive --epi1 the screening
# threshold is zero, so which side of it such a pair falls on is noise as
# well.  A chi-square statistic that small is a p-value of 0.999
# either way.
function abs(x) { return (x < 0)? -x : x }
function sigdigits(t,   s) {
    s = t; sub(/^-/, "", s); sub(/[eE].*$/, "", s); gsub(/\./, "", s)
    sub(/^0+/, "", s); sub(/0+$/, "", s)
    return (length(s) > 0)? length(s) : 1
}
function same(txt, a, b,   nd, mag, tol) {
    nd = sigdigits(txt)
    mag = (a == 0)? 0 : int(log(abs(a)) / log(10))
    if (abs(a) < 1) { mag = mag - 1 }
    tol = 0.55 * (10 ^ (mag - nd + 1))
    return (abs(a - b) <= tol)
}
function pairkey(x, y) { return (x < y)? (x "|" y) : (y "|" x) }
BEGIN { zero_tol = 1e-6 }
# PLINK 1.9 emits a row with STAT=nan for a pair whose statistic is undefined
# (a table the fit stage cannot handle), but only when --epi1 is permissive
# enough, and its own .summary excludes those
# pairs from N_TOT.  The port omits the row instead, so those rows are skipped
# here rather than being counted as a disagreement.
FNR == NR {
    if (FNR == 1) { has_df = ($6 == "DF"); next }
    if ($5 == "nan" || $5 == "-nan" || $5 == "inf" || $5 == "-inf") { ++nan1; next }
    k = pairkey($2, $4); stat[k] = $5 + 0; txt[k] = $5; df1[k] = $6
    if (abs($5 + 0) < zero_tol) { near_zero[k] = 1 } else { ++n1 }
    next
}
/^#/ { next }
{
    k = pairkey($2, $4)
    if (!(k in stat)) {
        if (abs($5 + 0) < zero_tol) { ++zero_only2; next }
        print "pair missing from the PLINK 1.9 report: " k; exit 1
    }
    seen[k] = 1
    if (!same(txt[k], stat[k], $5 + 0)) {
        if ((k in near_zero) && (abs($5 + 0) < zero_tol)) { ++loose; next }
        print "STAT differs on " k ": " txt[k] " vs " $5; exit 1
    }
    ++matched
    if (has_df && (df1[k] != $6)) {
        print "DF differs on " k ": " df1[k] " vs " $6; exit 1
    }
}
END {
    for (k in stat) {
        if ((!(k in seen)) && (!(k in near_zero))) {
            print "pair missing from the plink2 report: " k; exit 1
        }
    }
    if (matched + loose == 0) { print "no pairs tested, so nothing was compared"; exit 1 }
    msg = matched + 0 " pairs matched"
    if (loose > 0) { msg = msg ", " loose " zero-statistic rows within the noise floor" }
    if (zero_only2 > 0) { msg = msg ", " zero_only2 " zero-statistic rows unmatched" }
    if (nan1 > 0) { msg = msg ", " nan1 " undefined rows skipped" }
    print msg
}
