# Recomputes the --groupdist summary: for affected-affected (AA),
# affected-unaffected (AU) and unaffected-unaffected (UU) sample pairs, the
# mean distance, its standard error (reported by PLINK 1.9 as "sem"),
# and the median (the average of the two middle values for an
# even count).  Prints "AA mean se median", then AU and UU.
#
# Usage: awk -f groupdist.awk x.fam values.txt
#   values.txt: the lower-triangular distance matrix without the diagonal,
#               one value per line, row by row ((2,1), (3,1), (3,2), ...).

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
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            ++k
            g = (ph[i] == 2 && ph[j] == 2)? "AA" : (ph[i] == 2 || ph[j] == 2)? "AU" : "UU"
            vals[g, ++ct[g]] = v[k]
            s[g] += v[k]
            ss[g] += v[k] * v[k]
        }
    }
    split("AA AU UU", groups, " ")
    for (gi = 1; gi <= 3; ++gi) {
        g = groups[gi]
        m = ct[g]
        for (a = 1; a <= m; ++a) sorted[a] = vals[g, a]
        # insertion sort; the counts here are small
        for (a = 2; a <= m; ++a) {
            x = sorted[a]
            for (b = a - 1; b >= 1 && sorted[b] > x; --b) sorted[b + 1] = sorted[b]
            sorted[b + 1] = x
        }
        median = (m % 2)? sorted[(m + 1) / 2] : (sorted[m / 2] + sorted[m / 2 + 1]) / 2
        mean = s[g] / m
        printf "%s %.10g %.10g %.10g\n", g, mean, sqrt((ss[g] - m * mean * mean) / (m - 1) / m), median
    }
}
