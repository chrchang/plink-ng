# Compares a two-category --glm multinomial report (second file) against a
# --glm no-firth logistic report (first file) on the same data, keyed on ID and
# TEST.  Every per-category row of the multinomial report must have a logistic
# counterpart with the same A1, OBS_CT, odds ratio, standard error, Z statistic,
# p-value and error code; the LRT row has no counterpart and is only counted.
#
# The two fits stop on different convergence rules (the logistic one imitates
# R's glm.fit, on the deviance), so beyond the half-unit rounding allowance on
# each side, a relative slack of 1e-6 is allowed.
function abs(x) { return (x < 0)? -x : x }
function is_num(t) { return (t ~ /^-?[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?$/) }
function last_place(t,   m, e, pt) {
    m = t; sub(/^-/, "", m); e = 0
    if (match(m, /[eE]/)) { e = substr(m, RSTART + 1) + 0; m = substr(m, 1, RSTART - 1) }
    pt = index(m, ".")
    if (pt == 0) { return 10 ^ e }
    return 10 ^ (e - (length(m) - pt))
}
function same(t1, t2,   a, b, tol) {
    if (t1 == t2) { return 1 }
    if (!(is_num(t1) && is_num(t2))) { return 0 }
    a = t1 + 0; b = t2 + 0
    tol = 0.5 * (last_place(t1) + last_place(t2)) + rel_slack * (abs(a) + abs(b))
    return (abs(a - b) <= tol)
}
function fail(msg) { print msg; failed = 1; exit 1 }
BEGIN { FS = "\t"; rel_slack = 1e-6; field_ct = split("A1 OBS_CT OR LOG(OR)_SE P ERRCODE", fields, " ") }
FNR == 1 {
    split("", c)
    for (i = 1; i <= NF; ++i) { h = $i; sub(/^#/, "", h); c[h] = i }
    if (NR == 1) {
        for (f = 1; f <= field_ct; ++f) { lcol[f] = c[fields[f]] }
        lcol[field_ct + 1] = c["Z_STAT"]; lid = c["ID"]; ltest = c["TEST"]
    } else {
        for (f = 1; f <= field_ct; ++f) { mcol[f] = c[fields[f]] }
        mcol[field_ct + 1] = c["Z_OR_CHISQ_STAT"]; mid = c["ID"]; mtest = c["TEST"]; mcat = c["CATEGORY"]
    }
    next
}
FNR == NR { logistic[$lid "|" $ltest] = $0; next }
{
    if ($mtest ~ /^LRT_/) {
        if ($mtest != "LRT_1DF") { fail("unexpected joint test " $mtest) }
        ++lrt_ct
        next
    }
    if ($mcat != "hi") { fail("unexpected category " $mcat) }
    k = $mid "|" $mtest
    if (!(k in logistic)) { fail("no logistic row for " k) }
    split(logistic[k], o, "\t")
    for (f = 1; f <= field_ct + 1; ++f) {
        if (!same(o[lcol[f]], $mcol[f])) { fail("field " f " differs on " k ": " o[lcol[f]] " vs " $mcol[f]) }
    }
    ++matched
}
END {
    if (failed) { exit 1 }
    if ((matched == 0) || (lrt_ct == 0)) { print "no rows compared"; exit 1 }
    print matched " rows matched, " lrt_ct " LRT rows"
}
