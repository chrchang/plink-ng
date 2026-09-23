#!/usr/bin/env python3
# statsmodels MNLogit reference values for TEST_GLM_MULTINOMIAL (needs numpy
# and statsmodels; not run by the test itself, its outputs are committed).
#
# Mirrors what plink2 --glm multinomial omit-ref does:
# * A1 = the ALT alleles, dosages at plink2's 1/16384 resolution, samples with
#   a missing call dropped per variant.
# * On chrX, the SEX covariate is added (1 = male, 2 = female) and male
#   haploid calls are coded 0/2.  The non-additive models skip chrX (the
#   fixture has males).
# * Additive model: one row per variant, testing every ALT allele jointly
#   (constant ones left out).  Other models: one row per ALT allele A1, whose
#   recoded column(s) are tested, with the other non-constant ALT alleles as
#   additive nuisance predictors in both fits.
# * When a predictor's values in some level all lie at or below (or at or
#   above) its values in every other level, the maximum-likelihood estimate
#   does not exist: SEP (plink2's SEPARATION) for the likelihood ratio and
#   Wald tests, and for the score test too when the column is a nuisance
#   one.  Like plink2, this also checks two columns that are not predictors
#   but lie in their span: the omitted allele's count for a multiallelic
#   variant in the additive model, and whichever of ADD and HOM the
#   genotypic/hethom basis leaves out.  A constant tested column
#   (non-additive models), or no nonconstant
#   ALT allele (additive), gives CONST (CONST_ALLELE).
# * With 'firth', every fit maximizes Firth's penalized log-likelihood
#   l*(b) = l(b) + 0.5 log det I(b) instead (--glm multinomial firth), so
#   there is no separation check.  LRT is the penalized likelihood ratio
#   statistic, whose restricted fit maximizes the same l* with the tested
#   coefficients at zero; WALD and the standard errors use the inverse of the
#   leverage-augmented information matrix, logistf's covariance.  This does
#   not share plink2's algebra: I and each dI/db are built from explicit
#   per-sample Kronecker products, and l* is maximized with scipy's BFGS
#   followed by Newton steps on a finite-difference Hessian.
#
# Usage: oracle.py <phenotype column> <use covariates: 0/1> <model> <output> [firth]
# where <model> is one of add, dominant, recessive, hetonly, genotypic, hethom.
import gzip
import sys
import warnings
import numpy as np
import scipy.optimize
import statsmodels.api as sm

warnings.simplefilter("ignore")
pheno_name, use_covar, model, out_fname = sys.argv[1], sys.argv[2] == "1", sys.argv[3], sys.argv[4]
use_firth = (len(sys.argv) > 5) and (sys.argv[5] == "firth")

phe_rows = [l.split() for l in open("pheno.txt")]
col = phe_rows[0].index(pheno_name)
labels = [r[col] for r in phe_rows[1:]]
n = len(labels)
level_names = sorted(set(labels))
y_all = np.array([level_names.index(l) for l in labels])
cov = np.loadtxt("covar.txt", skiprows=1, usecols=range(1, 6))
sex = np.loadtxt("sex.txt", skiprows=1, usecols=1)

# variants: (chrom, id, allele names, per-allele count matrix with NaN rows
# for missing calls)
variants = []
for vcf_fname in ("mn.vcf.gz", "mn_multi.vcf.gz"):
    with gzip.open(vcf_fname, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            alleles = [t[3]] + t[4].split(",")
            fmt = t[8].split(":")
            counts = np.zeros((len(alleles), n))
            for i, cell in enumerate(t[9:]):
                parts = cell.split(":")
                if parts[0] in ("./.", "."):
                    counts[:, i] = np.nan
                elif "DS" in fmt:
                    ds = round(float(parts[1]) * 16384) / 16384
                    counts[1, i] = ds
                    counts[0, i] = 2 - ds
                elif "/" not in parts[0]:
                    counts[int(parts[0]), i] += 2
                else:
                    for a in parts[0].split("/"):
                        counts[int(a), i] += 1
            variants.append((t[0], t[2], alleles, counts))


def separated(x, yv, present):
    return any((x[yv == l].max() <= x[yv != l].min()) or (x[yv == l].min() >= x[yv != l].max()) for l in present)


def recode(x):
    domdev = np.where(x > 1, 2 - x, x)
    if model == "dominant":
        return [np.minimum(x, 1)]
    if model == "recessive":
        return [np.where(x < 1, 0, x - 1)]
    if model == "hetonly":
        return [domdev]
    if model == "genotypic":
        return [x, domdev]
    if model == "hethom":
        return [np.where(x < 1, 0, x - 1), domdev]
    return [x]


def fit(yv, X, start=None):
    return sm.MNLogit(yv, X).fit(method="newton", maxiter=200, disp=0, tol=1e-14, start_params=start)


class FirthModel:
    # yv: class indices (0 = reference); X: n x p (column 0 the intercept);
    # parameters are class-major, J x p.
    def __init__(self, yv, X):
        self.X = X
        self.n, self.p = X.shape
        self.J = yv.max()
        self.Y = np.zeros((self.n, self.J))
        for c in range(1, self.J + 1):
            self.Y[yv == c, c - 1] = 1.0

    def probs(self, theta):
        eta = np.hstack([np.zeros((self.n, 1)), self.X @ theta.reshape(self.J, self.p).T])
        m = eta.max(1, keepdims=True)
        lse = m[:, 0] + np.log(np.exp(eta - m).sum(1))
        ll = (eta[:, 1:] * self.Y).sum() - lse.sum()
        return ll, np.exp(eta[:, 1:] - lse[:, None])

    def info(self, P, wts=None):
        W = np.einsum("ia,ab->iab", P, np.eye(self.J)) - np.einsum("ia,ib->iab", P, P)
        if wts is not None:
            W = W * wts[:, None, None]
        n, J, p = self.n, self.J, self.p
        return np.einsum("iab,ik,il->akbl", W, self.X, self.X).reshape(J * p, J * p), W

    def lstar(self, theta):
        ll, P = self.probs(theta)
        sign, logdet = np.linalg.slogdet(self.info(P)[0])
        return ll + 0.5 * logdet if sign > 0 else -np.inf

    def grad(self, theta):
        ll, P = self.probs(theta)
        I, W = self.info(P)
        Iinv = np.linalg.inv(I)
        J, p = self.J, self.p
        g = ((self.Y - P).T @ self.X).ravel()
        for c in range(J):
            # dW_i/deta_ic, with dP_a/deta_c = P_a (1[a = c] - P_c)
            dP = P * ((np.arange(J) == c)[None, :] - P[:, c:c + 1])
            dW = np.einsum("ia,ab->iab", dP, np.eye(J)) - np.einsum("ia,ib->iab", dP, P) - np.einsum("ia,ib->iab", P, dP)
            for j in range(p):
                dI = np.einsum("i,iab,ik,il->akbl", self.X[:, j], dW, self.X, self.X).reshape(J * p, J * p)
                g[c * p + j] += 0.5 * np.sum(Iinv * dI)
        return g

    def aug_cov(self, theta):
        ll, P = self.probs(theta)
        I, W = self.info(P)
        Iinv = np.linalg.inv(I)
        J, p = self.J, self.p
        # h_i^(ab) = x_i' (I^-1)_(ab) x_i; leverage tr(W_i h_i)
        H = np.einsum("ik,akbl,il->iab", self.X, Iinv.reshape(J, p, J, p), self.X)
        lev = np.einsum("iab,iba->i", W, H)
        return np.linalg.inv(self.info(P, 1.0 + lev)[0])

    def fit(self, start, free):
        # maximizes l* over the parameters flagged in free, the others held
        # at their starting values
        theta = start.copy()

        def full(t):
            v = theta.copy()
            v[free] = t
            return v
        res = scipy.optimize.minimize(lambda t: -self.lstar(full(t)), theta[free], jac=lambda t: -self.grad(full(t))[free], method="BFGS", options={"gtol": 1e-9, "maxiter": 2000})
        t = res.x
        for _ in range(20):
            g = self.grad(full(t))[free]
            if np.max(np.abs(g)) < 1e-11:
                break
            h = 1e-5
            k = len(t)
            hess = np.zeros((k, k))
            for idx in range(k):
                e = np.zeros(k)
                e[idx] = h
                hess[:, idx] = (self.grad(full(t + e))[free] - self.grad(full(t - e))[free]) / (2 * h)
            t = t - np.linalg.solve(0.5 * (hess + hess.T), g)
        return full(t), self.lstar(full(t))


def fmt_list(vals):
    return ",".join("NA" if v is None else ("%.10g" % v) for v in vals)


with open(out_fname, "w") as fo:
    fo.write("#ID\tA1\tOBS_CT\tDF\tLRT\tSCORE\tWALD" + "".join("\tBETA_%s\tSE_%s" % (l, l) for l in level_names[1:]) + "\n")
    for chrom, vid, alleles, counts in variants:
        if chrom == "X" and model != "add":
            continue
        nm = ~np.isnan(counts[0])
        yv = y_all[nm]
        present = np.unique(yv)
        yc = np.searchsorted(present, yv)
        base = [np.ones(nm.sum())]
        z = []
        if use_covar:
            z.append(cov[nm])
        if chrom == "X":
            z.append(sex[nm, None])
        if z:
            zz = np.column_stack(z)
            base.append((zz - zz.mean(0)) / zz.std(0))
        X0base = np.column_stack(base)
        alt_counts = counts[1:, nm]
        if model == "add":
            rows = [(",".join(alleles[1:]), [], list(range(len(alleles) - 1)))]
        else:
            rows = [(alleles[a + 1], [b for b in range(len(alleles) - 1) if b != a], [a]) for a in range(len(alleles) - 1)]
        for a1_label, nuisance_idxs, a1_idxs in rows:
            df_slots = len(a1_idxs) if model == "add" else len(recode(alt_counts[0]))
            na_cells = "\t".join(fmt_list([None] * df_slots) + "\t" + fmt_list([None] * df_slots) for _ in level_names[1:])

            def emit(df, lr, score, wald, coef_cells):
                fo.write("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n" % (vid, a1_label, nm.sum(), df, lr, score, wald, coef_cells))

            nuisance = [alt_counts[b] for b in nuisance_idxs if alt_counts[b].min() != alt_counts[b].max()]
            if model == "add":
                tested = [(slot, alt_counts[a]) for slot, a in enumerate(a1_idxs) if alt_counts[a].min() != alt_counts[a].max()]
                if not tested:
                    emit("NA", "CONST", "CONST", "CONST", na_cells)
                    continue
            else:
                tested = list(enumerate(recode(alt_counts[a1_idxs[0]])))
                if any(c.min() == c.max() for _, c in tested):
                    emit("NA", "CONST", "CONST", "CONST", na_cells)
                    continue
            # plink2's correlation / VIF check on the non-intercept predictors
            preds = ([X0base[:, 1:]] if X0base.shape[1] > 1 else []) + [c[:, None] for c in nuisance] + [c[:, None] for _, c in tested]
            corr = np.corrcoef(np.column_stack(preds), rowvar=False)
            corr = np.atleast_2d(corr)
            if np.any(np.abs(corr[np.triu_indices_from(corr, 1)]) > 0.999):
                emit("NA", "CORR", "CORR", "CORR", na_cells)
                continue
            if corr.shape[0] > 1 and (np.linalg.cond(corr) > 1e14 or np.max(np.diag(np.linalg.inv(corr))) > 50):
                emit("NA", "VIF", "VIF", "VIF", na_cells)
                continue
            if (not use_firth) and any(separated(c, yv, present) for c in nuisance):
                emit("NA", "SEP", "SEP", "SEP", na_cells)
                continue
            separation_cols = [c for _, c in tested]
            if model == "add" and len(alleles) > 2:
                # the omitted (REF) allele's count is in the span of the tested
                # ones; plink2 checks it too
                separation_cols.append(counts[0, nm])
            if model in ("genotypic", "hethom"):
                # (ADD, DOMDEV) and (HOM, HET) span the same space; plink2 also
                # checks the column its basis leaves out
                x = alt_counts[a1_idxs[0]]
                separation_cols.append(np.where(x < 1, 0, x - 1) if model == "genotypic" else x)
            tested_sep = any(separated(c, yv, present) for c in separation_cols)
            X0 = np.column_stack([X0base] + [c - c.mean() for c in nuisance]) if nuisance else X0base
            X1 = np.column_stack([X0] + [c - c.mean() for _, c in tested])
            if use_firth:
                fm = FirthModel(yc, X1)
                k0 = X0.shape[1]
                k1 = X1.shape[1]
                nt = len(tested)
                ncls = len(present) - 1
                df = ncls * nt
                start = np.zeros((ncls, k1))
                for c in range(1, ncls + 1):
                    start[c - 1, 0] = np.log(np.mean(yc == c) / np.mean(yc == 0))
                start = start.ravel()
                free0 = np.array([(i % k1) < k0 for i in range(ncls * k1)])
                theta0, lstar0 = fm.fit(start, free0)
                theta1, lstar1 = fm.fit(theta0, np.ones(ncls * k1, dtype=bool))
                lr = 2 * (lstar1 - lstar0)
                covm = fm.aug_cov(theta1)
                idx = [c * k1 + k0 + j for c in range(ncls) for j in range(nt)]
                g = theta1[idx]
                wald = g @ np.linalg.solve(covm[np.ix_(idx, idx)], g)
                cells = []
                for lidx in range(1, len(level_names)):
                    betas = [None] * df_slots
                    ses = [None] * df_slots
                    if present[0] == 0 and lidx in present:
                        c = list(present).index(lidx) - 1
                        for j, (slot, _) in enumerate(tested):
                            betas[slot] = theta1[c * k1 + k0 + j]
                            ses[slot] = np.sqrt(covm[c * k1 + k0 + j, c * k1 + k0 + j])
                    cells.append(fmt_list(betas) + "\t" + fmt_list(ses))
                emit(df, "%.10g" % lr, "NA", "%.10g" % wald, "\t".join(cells))
                continue
            r0 = fit(yc, X0)
            m1 = sm.MNLogit(yc, X1)
            k0 = X0.shape[1]
            k1 = X1.shape[1]
            nt = len(tested)
            start = np.vstack([r0.params, np.zeros((nt, r0.params.shape[1]))]).ravel(order="F")
            sc = m1.score(start)
            score = sc @ np.linalg.solve(-m1.hessian(start), sc)
            ncls = len(present) - 1
            df = ncls * nt
            if tested_sep:
                emit(df, "SEP", "%.10g" % score, "SEP", na_cells)
                continue
            r1 = m1.fit(method="newton", maxiter=200, disp=0, tol=1e-14, start_params=start)
            lr = 2 * (r1.llf - r0.llf)
            params = np.asarray(r1.params)
            covm = np.asarray(r1.cov_params())
            idx = [c * k1 + k0 + j for c in range(ncls) for j in range(nt)]
            g = params[k0:, :].T.ravel()
            V = covm[np.ix_(idx, idx)]
            wald = g @ np.linalg.solve(V, g)
            cells = []
            for lidx in range(1, len(level_names)):
                betas = [None] * df_slots
                ses = [None] * df_slots
                if present[0] == 0 and lidx in present:
                    c = list(present).index(lidx) - 1
                    for j, (slot, _) in enumerate(tested):
                        betas[slot] = params[k0 + j, c]
                        ses[slot] = np.sqrt(covm[c * k1 + k0 + j, c * k1 + k0 + j])
                cells.append(fmt_list(betas) + "\t" + fmt_list(ses))
            emit(df, "%.10g" % lr, "%.10g" % score, "%.10g" % wald, "\t".join(cells))
