# Recomputes PLINK 1.9 --model statistics from the genotype counts in the
# GENO rows, and checks the TREND/ALLELIC/DOM/REC counts and every CHISQ
# against them.  With a second file (an --assoc counts or --assoc fisher
# report for the same data), also checks its counts, CHISQ and OR, or its
# Fisher's exact test P, against the allelic table.
#
# Usage: awk -f check_model.awk <.model> [<.assoc>|<.assoc.fisher>]
# Exits 1, printing the first few mismatches to stderr, on any mismatch.

function close_enough(a, b,    d, m) {
    if (a == "NA" || b == "NA") return a == b
    d = a - b
    if (d < 0) d = -d
    m = (a < 0)? -a : a
    if (b > m) m = b
    return d <= 1e-3 * m + 1e-6
}

function fail(msg) {
    if (++bad <= 5) print msg > "/dev/stderr"
}

# Pearson chi-square of a 2 x k table (rows r1[], r2[]).
function pearson(k,    i, t1, t2, n, col, e, x) {
    t1 = t2 = 0
    for (i = 1; i <= k; ++i) {
        t1 += r1[i]
        t2 += r2[i]
    }
    n = t1 + t2
    x = 0
    for (i = 1; i <= k; ++i) {
        col = r1[i] + r2[i]
        if (col == 0) continue
        e = t1 * col / n
        x += (r1[i] - e) ^ 2 / e
        e = t2 * col / n
        x += (r2[i] - e) ^ 2 / e
    }
    return x
}

# Cochran-Armitage trend chi-square on allele counts, as PLINK computes it.
# a[] = case genotype counts and u[] = control genotype counts, for
# (A1A1, A1A2, A2A2).
function trend(a, u,    n, r, s, ca, cu, t, v, cn, i, ct) {
    r = a[1] + a[2] + a[3]
    s = u[1] + u[2] + u[3]
    n = r + s
    t = 0
    v = 0
    cn = 0
    for (i = 1; i <= 3; ++i) {
        ct = a[i] + u[i]
        # score = number of A1 alleles
        t += (3 - i) * (a[i] * s - u[i] * r)
        v += (3 - i) ^ 2 * ct
        cn += (3 - i) * ct
    }
    v = r * s * (n * v - cn * cn) / n
    return t * t / v
}

function lfact(k) {
    return lf[k]
}

# Two-sided Fisher's exact test P for the 2 x 2 table (a b / c d).
function fisher(a, b, c, d,    r1t, c1t, n, lo, hi, x, lp, lpo, p) {
    r1t = a + b
    c1t = a + c
    n = a + b + c + d
    lo = c1t - (n - r1t)
    if (lo < 0) lo = 0
    hi = (r1t < c1t)? r1t : c1t
    lpo = lfact(r1t) + lfact(n - r1t) + lfact(c1t) + lfact(n - c1t) - lfact(n) - lfact(a) - lfact(b) - lfact(c) - lfact(d)
    p = 0
    for (x = lo; x <= hi; ++x) {
        lp = lfact(r1t) + lfact(n - r1t) + lfact(c1t) + lfact(n - c1t) - lfact(n) - lfact(x) - lfact(r1t - x) - lfact(c1t - x) - lfact(n - r1t - c1t + x)
        if (lp <= lpo + 1e-7) p += exp(lp)
    }
    return (p > 1)? 1 : p
}

BEGIN {
    lf[0] = 0
    for (i = 1; i <= 20000; ++i) lf[i] = lf[i - 1] + log(i)
}

FNR == 1 { ++file_idx; next }

file_idx == 1 {
    snp = $2
    test = $5
    if (test == "GENO") {
        split($6, ga, "/")
        split($7, gu, "/")
        geno_a[snp] = $6
        geno_u[snp] = $7
        al_a1[snp] = 2 * ga[1] + ga[2]
        al_a2[snp] = 2 * ga[3] + ga[2]
        al_u1[snp] = 2 * gu[1] + gu[2]
        al_u2[snp] = 2 * gu[3] + gu[2]
        r1[1] = ga[1]; r1[2] = ga[2]; r1[3] = ga[3]
        r2[1] = gu[1]; r2[2] = gu[2]; r2[3] = gu[3]
        if (!close_enough($8, pearson(3))) fail(snp " GENO: " $8 " vs " pearson(3))
        ++geno_ct
        next
    }
    if (!(snp in geno_a)) {
        fail(snp " " test " before GENO")
        next
    }
    split(geno_a[snp], ga, "/")
    split(geno_u[snp], gu, "/")
    if (test == "TREND" || test == "ALLELIC") {
        want_a = al_a1[snp] "/" al_a2[snp]
        want_u = al_u1[snp] "/" al_u2[snp]
        if ($6 != want_a || $7 != want_u) fail(snp " " test " counts: " $6 " " $7 " vs " want_a " " want_u)
        if (test == "TREND") {
            x = trend(ga, gu)
        } else {
            r1[1] = al_a1[snp]; r1[2] = al_a2[snp]
            r2[1] = al_u1[snp]; r2[2] = al_u2[snp]
            x = pearson(2)
        }
    } else if (test == "DOM" || test == "REC") {
        if (test == "DOM") {
            r1[1] = ga[1] + ga[2]; r1[2] = ga[3]
            r2[1] = gu[1] + gu[2]; r2[2] = gu[3]
        } else {
            r1[1] = ga[1]; r1[2] = ga[2] + ga[3]
            r2[1] = gu[1]; r2[2] = gu[2] + gu[3]
        }
        want_a = r1[1] "/" r1[2]
        want_u = r2[1] "/" r2[2]
        if ($6 != want_a || $7 != want_u) fail(snp " " test " counts: " $6 " " $7 " vs " want_a " " want_u)
        x = pearson(2)
    } else {
        next
    }
    if ($8 != "NA" && !close_enough($8, x)) fail(snp " " test ": " $8 " vs " x)
    ++row_ct
    next
}

# --assoc counts (CHISQ column) or --assoc fisher (P column, no CHISQ).
FNR == 1 { next }
{
    snp = $2
    if (!(snp in al_a1)) {
        fail(snp ": not in the --model report")
        next
    }
    a = al_a1[snp]; b = al_a2[snp]; c = al_u1[snp]; d = al_u2[snp]
    ++assoc_ct
    if (NF == 10) {
        # CHR SNP BP A1 C_A C_U A2 CHISQ P OR
        if ($5 != a || $6 != c) fail(snp " --assoc counts: " $5 " " $6 " vs " a " " c)
        r1[1] = a; r1[2] = b
        r2[1] = c; r2[2] = d
        if (!close_enough($8, pearson(2))) fail(snp " --assoc CHISQ: " $8 " vs " pearson(2))
        if (b * c > 0 && !close_enough($10, (a * d) / (b * c))) fail(snp " --assoc OR: " $10 " vs " (a * d) / (b * c))
    } else {
        # CHR SNP BP A1 F_A F_U A2 P OR
        if (!close_enough($8, fisher(a, b, c, d))) fail(snp " --assoc fisher P: " $8 " vs " fisher(a, b, c, d))
    }
}

END {
    if (geno_ct == 0 || row_ct == 0) fail("no --model rows checked")
    if (ARGC > 2 && assoc_ct == 0) fail("no --assoc rows checked")
    exit (bad > 0)
}
