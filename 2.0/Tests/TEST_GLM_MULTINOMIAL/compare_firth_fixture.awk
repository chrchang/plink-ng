# Compares a --glm firth multinomial report (second file) against the expected
# values in expected_firth.txt (first file), keyed on ID, TEST and CATEGORY.
# The expected values come from an independent reference (see run_tests.sh)
# and carry 10 significant digits, so this checks plink2's 6-digit output
# itself: a field may differ from the expected value by at most max_units
# units of its own last printed digit (0.5 for a correctly rounded value).
# P values go through plink2's float-precision p-value routines, which can be
# a little further off, so they are allowed p_units.
function abs(x) { return (x < 0)? -x : x }
function last_place(t,   m, e, pt) {
    m = t; sub(/^-/, "", m); e = 0
    if (match(m, /[eE]/)) { e = substr(m, RSTART + 1) + 0; m = substr(m, 1, RSTART - 1) }
    pt = index(m, ".")
    if (pt == 0) { return 10 ^ e }
    return 10 ^ (e - (length(m) - pt))
}
function fail(msg) { print msg; failed = 1; exit 1 }
BEGIN {
    FS = "\t"; max_units = 0.51; p_units = 0.75
    field_ct = split("OBS_CT OR LOG(OR)_SE Z_OR_CHISQ_STAT P", fields, " ")
}
FNR == 1 {
    split("", c)
    for (i = 1; i <= NF; ++i) { h = $i; sub(/^#/, "", h); c[h] = i }
    for (f = 1; f <= field_ct; ++f) {
        if (!(fields[f] in c)) { fail("column " fields[f] " missing") }
        if (NR == 1) { ecol[f] = c[fields[f]] } else { pcol[f] = c[fields[f]] }
    }
    if (NR == 1) { eid = c["ID"]; etest = c["TEST"]; ecat = c["CATEGORY"] }
    else { pid = c["ID"]; ptest = c["TEST"]; pcat = c["CATEGORY"]; perr = c["ERRCODE"] }
    next
}
FNR == NR { expected[$eid "|" $etest "|" $ecat] = $0; next }
{
    k = $pid "|" $ptest "|" $pcat
    if (!(k in expected)) { fail("unexpected row " k) }
    if ($perr != ".") { fail("error code " $perr " on " k) }
    split(expected[k], e, "\t")
    for (f = 1; f <= field_ct; ++f) {
        t = $pcol[f]; x = e[ecol[f]]
        if ((t == "NA") || (x == "NA")) {
            if (t != x) { fail(fields[f] " differs on " k ": " t " vs expected " x) }
            continue
        }
        u = abs(t - x) / last_place(t)
        if (u > ((fields[f] == "P")? p_units : max_units)) { fail(fields[f] " differs on " k ": " t " vs expected " x " (" u " units)") }
        if (u > worst[fields[f]]) { worst[fields[f]] = u }
    }
    seen[k] = 1
    ++matched
}
END {
    if (failed) { exit 1 }
    for (k in expected) {
        if (!(k in seen)) { print "row missing from the report: " k; exit 1 }
    }
    printf("%d rows matched; largest differences, in units of the last printed digit:", matched)
    for (f = 2; f <= field_ct; ++f) { printf(" %s %.3f", fields[f], worst[fields[f]]) }
    printf("\n")
}
