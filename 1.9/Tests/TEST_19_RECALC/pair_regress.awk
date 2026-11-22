# Recomputes --regress-distance/--regress-rel slopes: over all sample pairs,
# y = the pair's matrix value and x = the pair's average phenotype.  Prints
# the slope of y on x and of x on y.
#
# Usage: awk [-v diag=1] -f pair_regress.awk x.fam values.txt
#   x.fam: the samples, in matrix order (phenotypes in column 6).
#   values.txt: the lower-triangular matrix, one value per line, row by row
#               ((2,1), (3,1), (3,2), ...), or with diag=1 including the
#               diagonal ((1,1), (2,1), (2,2), ...).

NR == FNR {
    ph[FNR - 1] = $6
    n = FNR
    next
}
{
    v[FNR] = $1
}
END {
    k = 0
    for (i = 0; i < n; ++i) {
        for (j = 0; j < i + (diag? 1 : 0); ++j) {
            ++k
            if (j == i) continue
            x = (ph[i] + ph[j]) / 2
            y = v[k]
            ++m
            sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y
        }
    }
    cov = sxy - sx * sy / m
    printf "%.10g %.10g\n", cov / (sxx - sx * sx / m), cov / (syy - sy * sy / m)
}
