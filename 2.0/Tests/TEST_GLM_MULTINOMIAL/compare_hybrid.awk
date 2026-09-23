# Checks a firth-fallback (.glm.multinomial.hybrid) report, the third file,
# against the 'firth' report (first file) and the 'no-firth' report (second
# file) of the same run.  With the FIRTH? column removed, every row of a
# variant flagged FIRTH?=Y must be identical to its 'firth' row, and the
# 'no-firth' run must have reported that variant as NA; every row of a
# variant flagged N must be identical to its 'no-firth' row.  Both kinds of
# variant must be present.
function fail(msg) { print msg; failed = 1; exit 1 }
function strip(line, col,   n, t, i, out) {
    n = split(line, t, "\t")
    out = ""
    for (i = 1; i <= n; ++i) {
        if (i != col) { out = out ((out == "")? "" : "\t") t[i] }
    }
    return out
}
BEGIN { FS = "\t" }
FNR == 1 {
    ++file_idx
    for (i = 1; i <= NF; ++i) {
        if ($i == "ID") { id_col = i }
        if ($i == "TEST") { test_col = i }
        if ($i == "CATEGORY") { cat_col = i }
        if ($i == "ERRCODE") { err_col = i }
        if ($i == "FIRTH?") { firth_col = i }
    }
    if (file_idx < 3) {
        header[file_idx] = $0
    } else {
        if (!firth_col) { fail("no FIRTH? column in the hybrid report") }
        h = strip($0, firth_col)
        if ((h != header[1]) || (h != header[2])) { fail("headers differ") }
    }
    next
}
file_idx < 3 {
    k = $id_col "|" $test_col "|" $cat_col
    rows[file_idx, k] = $0
    if ((file_idx == 2) && ($err_col != ".")) { nofirth_na[$id_col] = 1 }
    next
}
{
    k = $id_col "|" $test_col "|" $cat_col
    line = strip($0, firth_col)
    if ($firth_col == "Y") {
        if (!($id_col in nofirth_na)) { fail("Firth fallback on " $id_col ", which 'no-firth' could fit") }
        if (line != rows[1, k]) { fail("FIRTH?=Y row differs from 'firth': " k) }
        y_vars[$id_col] = 1
    } else if ($firth_col == "N") {
        if (line != rows[2, k]) { fail("FIRTH?=N row differs from 'no-firth': " k) }
        n_vars[$id_col] = 1
    } else {
        fail("bad FIRTH? value on " k)
    }
    ++matched
}
END {
    if (failed) { exit 1 }
    for (v in y_vars) { ++y_ct }
    for (v in n_vars) { ++n_ct }
    if ((!y_ct) || (!n_ct)) { print "expected both FIRTH?=Y and FIRTH?=N variants"; exit 1 }
    print matched " rows matched (" y_ct " Firth-fallback variants, " n_ct " others)"
}
