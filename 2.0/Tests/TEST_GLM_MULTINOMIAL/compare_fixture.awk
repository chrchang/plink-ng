# Compares a --glm multinomial report (second file) against expected values
# computed with statsmodels (first file, see run_tests.sh), keyed on ID, TEST
# and CATEGORY, at 5 significant digits: a numeric field may differ from the
# expected value by 0.55 units of its 5th significant digit (half a unit, plus
# the rounding of plink2's own 6-digit output).
function abs(x) { return (x < 0)? -x : x }
function close5(t, e,   a, b, unit) {
    if (t == e) { return 1 }
    if ((t == "NA") || (e == "NA")) { return 0 }
    a = t + 0; b = e + 0
    if (b == 0) { return (a == 0) }
    unit = 10 ^ (int(log(abs(b)) / log(10) + 100) - 100 - 4)
    return (abs(a - b) <= 0.55 * unit)
}
function fail(msg) { print msg; failed = 1; exit 1 }
BEGIN { FS = "\t"; field_ct = split("A1 OBS_CT OR LOG(OR)_SE Z_OR_CHISQ_STAT P", fields, " ") }
FNR == 1 {
    split("", c)
    for (i = 1; i <= NF; ++i) { h = $i; sub(/^#/, "", h); c[h] = i }
    for (f = 1; f <= field_ct; ++f) {
        if (!(fields[f] in c)) { fail("column " fields[f] " missing") }
        if (NR == 1) { ecol[f] = c[fields[f]] } else { pcol[f] = c[fields[f]] }
    }
    if (NR == 1) { eid = c["ID"]; etest = c["TEST"]; ecat = c["CATEGORY"] }
    else { pid = c["ID"]; ptest = c["TEST"]; pcat = c["CATEGORY"] }
    next
}
FNR == NR { expected[$eid "|" $etest "|" $ecat] = $0; next }
{
    k = $pid "|" $ptest "|" $pcat
    if (!(k in expected)) { fail("unexpected row " k) }
    split(expected[k], e, "\t")
    for (f = 1; f <= field_ct; ++f) {
        if (!close5($pcol[f], e[ecol[f]])) { fail(fields[f] " differs on " k ": " $pcol[f] " vs expected " e[ecol[f]]) }
    }
    seen[k] = 1
    ++matched
}
END {
    if (failed) { exit 1 }
    for (k in expected) {
        if (!(k in seen)) { print "row missing from the report: " k; exit 1 }
    }
    print matched " rows matched"
}
