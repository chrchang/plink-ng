# Compares columns of two whitespace-delimited files, line by line, with a
# relative tolerance for numbers (PLINK 1.9 often prints only 4 significant
# digits).  Exits 1, printing the first few mismatches to stderr, if any
# compared field differs or the files have different line counts.
#
# Usage: awk -v c1=<cols> -v c2=<cols> [-v skip1=<n>] [-v skip2=<n>]
#            [-v tol=<x>] [-v abstol=<x>] -f compare_cols.awk file1 file2
#   c1, c2: comma-separated 1-based column numbers, compared pairwise.
#   skip1, skip2: header lines to skip in each file (default 1).
#   tol: relative tolerance (default 1e-3).
#   abstol: values this close always match (default 1e-6).  Use this for
#           output printed with a fixed number of decimal places.
# Non-numeric fields (including NA/nan/inf) must match exactly, ignoring case.

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
    return tolower(a) == tolower(b)
}

BEGIN {
    if (tol == "") tol = 1e-3
    if (abstol == "") abstol = 1e-6
    if (skip1 == "") skip1 = 1
    if (skip2 == "") skip2 = 1
    n = split(c1, col1, ",")
    if (split(c2, col2, ",") != n) {
        print "compare_cols.awk: c1 and c2 have different lengths" > "/dev/stderr"
        exit 2
    }
}

FNR == 1 { ++file_idx }

file_idx == 1 {
    if (FNR > skip1) line1[++ct1] = $0
    next
}

FNR > skip2 {
    ++ct2
    if (ct2 > ct1) next
    split(line1[ct2], f1)
    for (i = 1; i <= n; ++i) {
        a = f1[col1[i]]
        b = $(col2[i])
        if (!same(a, b)) {
            if (++bad <= 5) {
                printf "line %d: column %d '%s' vs column %d '%s'\n", ct2, col1[i], a, col2[i], b > "/dev/stderr"
            }
        }
    }
}

END {
    if (ct1 != ct2) {
        printf "line counts differ: %d vs %d\n", ct1, ct2 > "/dev/stderr"
        exit 1
    }
    exit (bad > 0)
}
