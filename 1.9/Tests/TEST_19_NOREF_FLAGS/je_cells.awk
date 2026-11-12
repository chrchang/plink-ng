# Recomputes which variants and pairs "--fast-epistasis joint-effects
# --je-cellmin <cellmin>" tests, from a --recode A (.raw) file whose
# PHENOTYPE column is case/control (2 = case).
#
# A variant is skipped unless each of its three genotypes is seen in at least
# 3 * cellmin cases and 3 * cellmin controls.  A pair of remaining variants
# is tested only if each of the 18 two-locus genotype x status cells (samples
# with both calls) holds at least cellmin samples.  Prints "SITES <skipped
# variant count>", then "<ID1> <ID2> <1 if tested, 0 if not>" for every pair
# of remaining variants, in order.

FNR == 1 {
    for (k = 7; k <= NF; ++k) {
        id[k - 6] = $k
        sub(/_[^_]*$/, "", id[k - 6])
    }
    m = NF - 6
    next
}

$6 == 1 || $6 == 2 {
    ++n
    st[n] = $6 - 1
    for (k = 7; k <= NF; ++k) g[n, k - 6] = $k
}

END {
    skipped = 0
    kept_ct = 0
    for (j = 1; j <= m; ++j) {
        split("", c)
        for (i = 1; i <= n; ++i) {
            if (g[i, j] != "NA") ++c[st[i], g[i, j]]
        }
        ok = 1
        for (s = 0; s <= 1; ++s) {
            for (x = 0; x <= 2; ++x) {
                if (c[s, x] < 3 * cellmin) ok = 0
            }
        }
        if (ok) {
            kept[++kept_ct] = j
        } else {
            ++skipped
        }
    }
    print "SITES", skipped
    for (a = 1; a < kept_ct; ++a) {
        for (b = a + 1; b <= kept_ct; ++b) {
            j1 = kept[a]
            j2 = kept[b]
            split("", c)
            for (i = 1; i <= n; ++i) {
                if ((g[i, j1] != "NA") && (g[i, j2] != "NA")) ++c[st[i], g[i, j1], g[i, j2]]
            }
            ok = 1
            for (s = 0; s <= 1; ++s) {
                for (x = 0; x <= 2; ++x) {
                    for (y = 0; y <= 2; ++y) {
                        if (c[s, x, y] < cellmin) ok = 0
                    }
                }
            }
            print id[j1], id[j2], ok
        }
    }
}
