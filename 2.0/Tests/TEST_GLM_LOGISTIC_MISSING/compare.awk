# Compares a --glm report against an oracle report of the same variants,
# row by row, keyed on ID and TEST.
#
# Text fields must match exactly.  Numbers are printed to six significant
# digits and both sides are rounded, so two values that agree to far more
# digits than that can still print a whole unit of the last place apart when
# they straddle a rounding boundary.  A numeric field is therefore allowed half
# of the last printed place on each side, plus a relative slack far below a
# printed digit for the two regressions having been computed differently.
#
# One text difference is allowed.  A genotype column with no variance (e.g. a
# recessive coding of a variant with no hom-ALT call) is reported as
# VIF_INFINITE by the regression path that handles missing calls, and as
# RANK_DEFICIENT by the one the oracle takes, which has no missing calls to
# handle.  Either way the row has no estimates.
#
# awk still runs END after an exit from a main rule, so a failure sets failed
# and END stops at once.
function abs(x) { return (x < 0)? -x : x }
function is_num(t) { return (t ~ /^-?[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?$/) }
# The place value of the last printed digit, read off the text.
function last_place(t,   m, e, pt) {
    m = t; sub(/^-/, "", m); e = 0
    if (match(m, /[eE]/)) { e = substr(m, RSTART + 1) + 0; m = substr(m, 1, RSTART - 1) }
    pt = index(m, ".")
    if (pt == 0) { return 10 ^ e }
    return 10 ^ (e - (length(m) - pt))
}
# Integers (OBS_CT) have no rounding to allow for.
function same(t1, t2,   a, b, tol) {
    if (t1 == t2) { return 1 }
    if (!(is_num(t1) && is_num(t2))) { return 0 }
    if ((t1 ~ /^-?[0-9]+$/) || (t2 ~ /^-?[0-9]+$/)) { return 0 }
    a = t1 + 0; b = t2 + 0
    tol = 0.5 * (last_place(t1) + last_place(t2)) + rel_slack * (abs(a) + abs(b))
    return (abs(a - b) <= tol)
}
function same_errcode(e1, e2) {
    if (e1 == e2) { return 1 }
    return ((e1 == "VIF_INFINITE") && (e2 == "RANK_DEFICIENT")) || ((e1 == "RANK_DEFICIENT") && (e2 == "VIF_INFINITE"))
}
function fail(msg) { print msg; failed = 1; exit 1 }
BEGIN { FS = "\t"; rel_slack = 1e-8 }
FNR == 1 {
    if (NR == 1) { header = $0 }
    else if ($0 != header) { fail("headers differ: " header " vs " $0) }
    for (i = 1; i <= NF; ++i) {
        if ($i == "ID") { id_col = i }
        if ($i == "TEST") { test_col = i }
        if ($i == "ERRCODE") { errcode_col = i }
    }
    next
}
FNR == NR {
    k = $id_col "|" $test_col
    if (k in oracle) { fail("duplicate oracle row: " k) }
    oracle[k] = $0
    next
}
{
    k = $id_col "|" $test_col
    if (!(k in oracle)) { fail("row missing from the oracle: " k) }
    n = split(oracle[k], o, "\t")
    if (n != NF) { fail("field count differs on " k) }
    for (i = 1; i <= NF; ++i) {
        if ((i == errcode_col)? (!same_errcode(o[i], $i)) : (!same(o[i], $i))) { fail("field " i " differs on " k ": " o[i] " vs " $i) }
    }
    seen[k] = 1
    ++matched
}
END {
    if (failed) { exit 1 }
    for (k in oracle) {
        if (!(k in seen)) { print "row missing from the report: " k; exit 1 }
    }
    if (matched == 0) { print "no rows compared"; exit 1 }
    print matched " rows matched"
}
