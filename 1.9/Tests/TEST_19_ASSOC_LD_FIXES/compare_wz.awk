# Compares the P(WZ) column of a PLINK 1.9 --meta-analysis 'weighted-z'
# report (second file) with plink2's P_WZ (first file), matching rows by
# variant ID.  Values are compared as log10s parsed from the text, since
# plink2 prints p-values below the double range (e.g. 7e-846), which PLINK
# 1.9 must print as 0 or another value below 1e-300.  PLINK 1.9 prints 4
# significant digits.
#
# Exits nonzero on any mismatch, or unless at least 5 rows have
# 1e-300 < p < 1e-16 (the range where the old formula printed 0).

function lg10(x,    a, n, m) {
    n = split(x, a, /[eE]/)
    m = a[1] + 0
    if (m <= 0) return -1e9
    return log(m) / log(10) + ((n > 1)? a[2] + 0 : 0)
}

NR == FNR {
    if (FNR > 1) ref[$3] = $14
    next
}

FNR > 1 {
    if (!($3 in ref)) exit 1
    lr = lg10(ref[$3])
    l19 = lg10($14)
    if (lr < -307) {
        if (l19 > -300) exit 1
        next
    }
    d = l19 - lr
    if (d < 0) d = -d
    if (d > 2.2e-4) exit 1
    if ((lr < -16) && (lr > -300)) ++small
    ++rows
}

END { exit !((rows > 0) && (small >= 5)) }
