# Compares two whitespace-delimited text files token by token, with a
# relative tolerance for numbers, for reports whose lines don't share a
# column layout.  Exits 1, printing the first few mismatches to stderr, if the
# files differ in line count, any line differs in token count, or any token
# differs; also if both files are empty.
#
# Usage: awk [-v tol=<x>] [-v abstol=<x>] -f compare_tokens.awk file1 file2
#   tol: relative tolerance (default 1e-3).
#   abstol: values this close always match (default 1e-6).
# Non-numeric tokens must match exactly.

function isnum(x) {
    return x ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/
}

function same(a, b,    d, m) {
    if (isnum(a) && isnum(b)) {
        a += 0
        b += 0
        d = a - b
        if (d < 0) d = -d
        if (d <= abstol) return 1
        m = (a < 0)? -a : a
        if (b > m) m = b
        if (-b > m) m = -b
        return d <= tol * m
    }
    return a == b
}

function fail(msg) {
    if (++bad <= 5) print msg > "/dev/stderr"
}

BEGIN {
    if (tol == "") tol = 1e-3
    if (abstol == "") abstol = 1e-6
}

FNR == 1 { ++file_idx }

file_idx == 1 {
    line1[++ct1] = $0
    next
}

{
    ++ct2
    if (ct2 > ct1) next
    n1 = split(line1[ct2], f1)
    if (n1 != NF) {
        fail("line " ct2 ": " n1 " vs " NF " tokens")
        next
    }
    for (i = 1; i <= NF; ++i) {
        if (!same(f1[i], $i)) fail("line " ct2 ": '" f1[i] "' vs '" $i "'")
    }
}

END {
    if (ct1 != ct2) fail("line counts differ: " ct1 " vs " ct2)
    if (!ct1) fail("nothing to compare")
    exit (bad > 0)
}
