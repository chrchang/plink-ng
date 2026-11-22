# Recomputes chrX --r2 values from a --recode A file (A1 dosages; males'
# chrX calls are 0 or 2 there) under an --ld-xchr coding: with mode=1, males
# are coded 0/1, and with mode=2, 0/2.  r^2 is the squared Pearson
# correlation of the codes over samples genotyped at both variants.  Prints
# "SNP_A SNP_B R2" for each pair at most `window` variants apart (default 9,
# matching --ld-window 10).
#
# Usage: awk -v mode=<1|2> [-v window=<n>] -f xchr_r2.awk x.raw

NR == 1 {
    if (window == "") window = 9
    for (i = 7; i <= NF; ++i) {
        n = $i
        sub(/_[^_]*$/, "", n)
        id[i - 6] = n
    }
    nv = NF - 6
    next
}
{
    ++ns
    for (j = 1; j <= nv; ++j) {
        v = $(j + 6)
        if (v != "NA" && $5 == 1 && mode == 1) v /= 2
        g[ns, j] = v
    }
}
END {
    for (a = 1; a < nv; ++a) {
        for (b = a + 1; b <= nv && b - a <= window; ++b) {
            n = sx = sy = sxx = syy = sxy = 0
            for (i = 1; i <= ns; ++i) {
                x = g[i, a]
                y = g[i, b]
                if (x == "NA" || y == "NA") continue
                ++n
                sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y
            }
            cov = sxy - sx * sy / n
            vx = sxx - sx * sx / n
            vy = syy - sy * sy / n
            print id[a], id[b], (vx > 0 && vy > 0)? cov * cov / (vx * vy) : "nan"
        }
    }
}
