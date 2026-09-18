# Compares a PLINK 1.9 .epi.cc report against a plink2 one.
#
# The pair key is normalized, since the two programs emit a pair in whatever
# order their scans reach it.  Both programs print STAT to six significant
# digits, and both printed values are rounded, so two statistics that agree to
# far more digits than that can still print a whole unit of the last place
# apart when they straddle a rounding boundary: 0.65577450000006 prints as
# 0.655775 and 0.65577449999999 as 0.655774.  The comparison therefore allows
# half of the last printed place on each side, plus a relative slack far below
# a printed digit for the difference between the two fits themselves.
#
# A statistic that is mathematically zero is the exception: it comes back as
# rounding noise, and the noise is absolute rather than relative, since the
# statistic is a sum over 2N cell terms.  Such a row is accepted when both
# sides are within zero_tol of zero, and a near-zero row present on one side
# alone is not a disagreement either: under a permissive --epi1 the screening
# threshold is zero, so which side of it such a pair falls on is noise as
# well.  A chi-square statistic that small is a p-value of 0.999
# either way.
#
# awk still runs END after an exit from a main rule, so a failure sets failed
# and END stops at once; otherwise END's missing-pair scan would run on a
# report that was only read part of the way, and name a pair that is there.
function abs(x) { return (x < 0)? -x : x }
function sigdigits(t,   s) {
    s = t; sub(/^-/, "", s); sub(/[eE].*$/, "", s); gsub(/\./, "", s)
    sub(/^0+/, "", s); sub(/0+$/, "", s)
    return (length(s) > 0)? length(s) : 1
}
# The place value of the last of stat_digits significant digits in a printed
# number, read off the text so that no log() rounding enters; 0 for a zero.
function last_place(t,   m, e, pt, lead) {
    m = t; sub(/^-/, "", m); e = 0
    if (match(m, /[eE]/)) { e = substr(m, RSTART + 1) + 0; m = substr(m, 1, RSTART - 1) }
    pt = index(m, ".")
    if (pt == 0) { m = m "."; pt = length(m) }
    if (match(m, /[1-9]/) == 0) { return 0 }
    # lead is the decimal exponent of the leading nonzero digit.
    lead = (RSTART < pt)? (pt - RSTART - 1) : (pt - RSTART)
    return 10 ^ (lead + e - stat_digits + 1)
}
function same(txt1, txt2, a, b,   tol) {
    tol = 0.5 * (last_place(txt1) + last_place(txt2)) + rel_slack * (abs(a) + abs(b))
    return (abs(a - b) <= tol)
}
function fail(msg) { print msg; failed = 1; exit 1 }
function pairkey(x, y) { return (x < y)? (x "|" y) : (y "|" x) }
BEGIN { zero_tol = 1e-6; stat_digits = 6; rel_slack = 1e-9 }
# PLINK 1.9 emits a row with STAT=nan for a pair whose statistic is undefined
# (a table the fit stage cannot handle), but only when --epi1 is permissive
# enough, and its own .summary excludes those
# pairs from N_TOT.  The port omits the row instead, so those rows are skipped
# here rather than being counted as a disagreement.
FNR == NR {
    if (FNR == 1) { has_df = ($6 == "DF"); next }
    if ($5 == "nan" || $5 == "-nan" || $5 == "inf" || $5 == "-inf") { ++nan1; next }
    if (sigdigits($5) > stat_digits) { fail("PLINK 1.9 STAT has more than " stat_digits " significant digits: " $5) }
    k = pairkey($2, $4); stat[k] = $5 + 0; txt[k] = $5; df1[k] = $6
    if (abs($5 + 0) < zero_tol) { near_zero[k] = 1 } else { ++n1 }
    next
}
/^#/ { next }
{
    if (sigdigits($5) > stat_digits) { fail("plink2 STAT has more than " stat_digits " significant digits: " $5) }
    k = pairkey($2, $4)
    if (!(k in stat)) {
        if (abs($5 + 0) < zero_tol) { ++zero_only2; next }
        fail("pair missing from the PLINK 1.9 report: " k)
    }
    seen[k] = 1
    if (!same(txt[k], $5, stat[k], $5 + 0)) {
        if ((k in near_zero) && (abs($5 + 0) < zero_tol)) { ++loose; next }
        fail("STAT differs on " k ": " txt[k] " vs " $5)
    }
    ++matched
    if (has_df && (df1[k] != $6)) {
        fail("DF differs on " k ": " df1[k] " vs " $6)
    }
}
END {
    if (failed) { exit 1 }
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
