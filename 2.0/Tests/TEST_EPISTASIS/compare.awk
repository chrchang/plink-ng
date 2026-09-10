# Compares a PLINK 1.9 .epi.cc/.epi.co report against a plink2 one.
#
# The pair key is normalized, since the two programs emit a pair in whatever
# order their scans reach it.  PLINK 1.9 prints STAT to four significant
# digits, so the comparison allows half of its last printed place rather than a
# fixed tolerance.  A statistic that is mathematically zero comes back as
# log(1) rounding noise, around 1e-30, from both programs independently; those
# are treated as equal.
function abs(x) { return (x < 0)? -x : x }
function sigdigits(t,   s) {
    s = t; sub(/^-/, "", s); sub(/[eE].*$/, "", s); gsub(/\./, "", s)
    sub(/^0+/, "", s); sub(/0+$/, "", s)
    return (length(s) > 0)? length(s) : 1
}
function same(txt, a, b,   nd, mag, tol) {
    if (abs(a) < 1e-20 && abs(b) < 1e-20) { return 1 }
    nd = sigdigits(txt)
    mag = (a == 0)? 0 : int(log(abs(a)) / log(10))
    if (abs(a) < 1) { mag = mag - 1 }
    tol = 0.55 * (10 ^ (mag - nd + 1))
    return (abs(a - b) <= tol)
}
function pairkey(x, y) { return (x < y)? (x "|" y) : (y "|" x) }
# PLINK 1.9 emits a row with STAT=nan for a pair whose statistic is undefined
# (an empty allele-table cell with 'no-ueki'), but only when --epi1 is
# permissive enough, and its own .summary excludes those pairs from N_TOT.  The
# port omits the row instead, so those rows are skipped here rather than being
# counted as a disagreement.
FNR == NR {
    if (FNR > 1) {
        if ($5 == "nan" || $5 == "-nan" || $5 == "inf" || $5 == "-inf") { ++nan1; next }
        k = pairkey($2, $4); stat[k] = $5 + 0; txt[k] = $5; ++n1
    }
    next
}
/^#/ { next }
{
    ++n2
    k = pairkey($2, $4)
    if (!(k in stat)) { print "pair missing from the PLINK 1.9 report: " k; exit 1 }
    if (!same(txt[k], stat[k], $5 + 0)) {
        print "STAT differs on " k ": " txt[k] " vs " $5; exit 1
    }
}
END {
    if (n1 != n2) { print "pair count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no pairs tested, so nothing was compared"; exit 1 }
    print n1 " pairs matched" ((nan1 > 0)? (", " nan1 " undefined rows skipped") : "")
}
