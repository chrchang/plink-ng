# Recomputes PLINK 1.9 --lasso (with an explicit minimum lambda), printing
# "<ID> <effect>" for every nonzero penalized coefficient: selected
# covariates first (as COV:<name>), then variants, in input order.
#
# Inputs, in order:
#   1. covariate file (FID IID cov1 cov2 ..., with header), or /dev/null;
#   2. --recode A output (.raw), whose PHENOTYPE column is the phenotype.
# Variables: minlambda; cc=1 for a case/control phenotype (2 = case);
# n_unsel = number of leading covariate columns that are unpenalized (as
# without --lasso-select-covars); the remaining covariates are penalized.
# All samples must have nonmissing phenotypes and covariates.
#
# PLINK standardizes every column (variants: A1 counts, with missing calls
# set to the mean) and the phenotype to mean 0 and norm 1, then runs
# coordinate descent over 100 log-spaced lambdas from the largest |x_j . y|
# down to minlambda.  Coefficients start at x_j . y and are carried over
# between lambdas; within a lambda, a coefficient that reaches 0 drops out
# until the next lambda; a lambda ends when the objective changes by a
# relative 1e-4 or less between sweeps.  The arithmetic below follows PLINK's
# operation order, so the results agree to rounding error.

FILENAME == ARGV[1] {
    if (FNR == 1) {
        for (k = 3; k <= NF; ++k) cov_name[k - 2] = $k
        cov_ct = NF - 2
        next
    }
    for (k = 3; k <= NF; ++k) cov_val[$1 SUBSEP $2, k - 2] = $k
    next
}

FNR == 1 {
    for (k = 7; k <= NF; ++k) {
        var_id[k - 6] = $k
        sub(/_[^_]*$/, "", var_id[k - 6])
    }
    var_ct = NF - 6
    next
}

{
    ++n
    key[n] = $1 SUBSEP $2
    y[n] = cc? ($6 == 2) : $6
    for (k = 7; k <= NF; ++k) g[n, k - 6] = $k
}

function standardize_cov(c, col,    i, sum, ssq, sub_by, mult, x) {
    sum = 0
    ssq = 0
    for (i = 1; i <= n; ++i) {
        x = cov_val[key[i], c]
        sum += x
        ssq += x * x
    }
    sub_by = sum * (1.0 / n)
    mult = sqrt_n_recip * sqrt((n - 1) / (ssq - sum * sub_by))
    for (i = 1; i <= n; ++i) X[col, i] = (cov_val[key[i], c] - sub_by) * mult
    label[col] = "COV:" cov_name[c]
}

END {
    sqrt_n_recip = sqrt(1.0 / n)
    # phenotype
    sum = 0
    ssq = 0
    for (i = 1; i <= n; ++i) {
        sum += y[i]
        ssq += y[i] * y[i]
    }
    mean = sum / n
    mult = sqrt_n_recip * sqrt((n - 1) / (ssq - sum * mean))
    for (i = 1; i <= n; ++i) y[i] = (y[i] - mean) * mult
    # columns: unpenalized covariates, penalized covariates, variants
    col_ct = 0
    for (c = 1; c <= cov_ct; ++c) standardize_cov(c, ++col_ct)
    for (j = 1; j <= var_ct; ++j) {
        cnt[0] = cnt[1] = cnt[2] = 0
        nm = 0
        for (i = 1; i <= n; ++i) {
            if (g[i, j] != "NA") {
                ++cnt[g[i, j]]
                ++nm
            }
        }
        if ((!cnt[2] && (!cnt[1] || !cnt[0])) || (!cnt[1] && !cnt[0])) continue
        s = 2 * cnt[2] + cnt[1]
        mean = s / nm
        mult = sqrt_n_recip * sqrt((nm - 1) / (4 * cnt[2] + cnt[1] - s * mean))
        val[0] = (-mean) * mult
        val[1] = (1 - mean) * mult
        val[2] = (2 - mean) * mult
        ++col_ct
        for (i = 1; i <= n; ++i) X[col_ct, i] = (g[i, j] == "NA")? 0 : val[g[i, j]]
        label[col_ct] = var_id[j]
    }
    lambda_max = 0
    for (c = 1; c <= col_ct; ++c) {
        d = 0
        for (i = 1; i <= n; ++i) d += X[c, i] * y[i]
        xhat[c] = d
        if (d < 0) d = -d
        if (d > lambda_max) lambda_max = d
    }
    loghi = log(lambda_max)
    loglo = log(minlambda)
    logdelta = (loghi - loglo) / 99.0
    for (lambi = 0; lambi < 100; ++lambi) {
        lambda = exp(loghi - logdelta * lambi)
        for (i = 1; i <= n; ++i) r[i] = y[i]
        for (c = 1; c <= col_ct; ++c) {
            d = -xhat[c]
            for (i = 1; i <= n; ++i) r[i] += X[c, i] * d
        }
        for (c = 1; c <= col_ct; ++c) active[c] = 1
        iter = 0
        while (1) {
            for (c = 1; c <= col_ct; ++c) {
                if (!active[c]) continue
                old = xhat[c]
                d = old
                for (i = 1; i <= n; ++i) d += X[c, i] * r[i]
                if (c > n_unsel) {
                    if (d > 0.0) {
                        d = (d - lambda > 0.0)? d - lambda : 0.0
                    } else {
                        d = (d + lambda < 0.0)? d + lambda : 0.0
                    }
                }
                xhat[c] = d
                if (d == 0.0) active[c] = 0
                d -= old
                for (i = 1; i <= n; ++i) r[i] -= X[c, i] * d
            }
            err_last = err_cur
            err_cur = 0.0
            for (c = 1; c <= col_ct; ++c) {
                if (active[c]) err_cur += (xhat[c] < 0)? -xhat[c] : xhat[c]
            }
            err_cur *= lambda
            for (i = 1; i <= n; ++i) err_cur += r[i] * r[i]
            if (iter++) {
                lo = (err_last < err_cur)? err_last : err_cur
                hi = (err_last < err_cur)? err_cur : err_last
                if (1.0 - lo / hi < 0.0001) break
            }
        }
    }
    for (c = n_unsel + 1; c <= col_ct; ++c) {
        if (xhat[c] != 0) printf "%s %.10g\n", label[c], xhat[c]
    }
}
