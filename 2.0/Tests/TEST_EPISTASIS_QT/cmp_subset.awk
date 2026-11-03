# Checks that every row of a set-mode .epi.qt report (second file) is the same
# test as the matching pair in the full all-pairs report (first file), and that
# the set-mode report has exactly expected_n rows.  A pair is looked up with
# its variants in either order, since the set scans reach it from whichever
# side their row set puts first.  Columns 5-8 (BETA_INT, SE, T_STAT, P) are
# compared to a relative 1e-9, which is far below a printed digit but allows
# for the row and column roles swapping.
function abs(x) { return (x < 0)? -x : x }
function close_enough(a, b) { return (abs(a - b) <= 1e-9 * (abs(a) + abs(b)) + 1e-300) }
FNR == NR {
    if (/^#/) { next }
    full[$2 "|" $4] = $0
    next
}
/^#/ { next }
{
    ++n
    if (($2 "|" $4) in full) { line = full[$2 "|" $4] } else if (($4 "|" $2) in full) { line = full[$4 "|" $2] } else {
        print "pair missing from the full report: " $2 " " $4; exit 1
    }
    split(line, f, "\t")
    split($0, g, "\t")
    for (i = 5; i <= 8; ++i) {
        if (!close_enough(f[i] + 0, g[i] + 0)) {
            print "column " i " differs on " $2 " " $4 ": " f[i] " vs " g[i]; exit 1
        }
    }
}
END {
    if (n != expected_n) { print "expected " expected_n " rows, found " n; exit 1 }
    print n " rows matched the full report"
}
