#!/usr/bin/env python3
# statsmodels MNLogit reference values for TEST_GLM_MULTINOMIAL (needs numpy
# and statsmodels; not run by the test itself, its outputs are committed).
#
# Mirrors what plink2 --glm multinomial omit-ref does: A1 = ALT, dosages at
# plink2's 1/16384 resolution, samples with a missing call dropped per
# variant, and on chrX the SEX covariate added (1 = male, 2 = female) and male
# haploid calls coded 0/2.  When the dosages of some level all lie at or
# below (or at or above) every dosage of the other levels, the
# maximum-likelihood estimate does not exist; the variant is reported as SEP
# (plink2's SEPARATION) for the likelihood ratio and Wald tests, and the score
# test is still computed.
#
# Usage: oracle.py <phenotype column> <use covariates: 0/1> <output file>
import gzip
import sys
import warnings
import numpy as np
import statsmodels.api as sm

warnings.simplefilter("ignore")
pheno_name, use_covar, out_fname = sys.argv[1], sys.argv[2] == "1", sys.argv[3]

phe_rows = [l.split() for l in open("pheno.txt")]
col = phe_rows[0].index(pheno_name)
labels = [r[col] for r in phe_rows[1:]]
n = len(labels)
level_names = sorted(set(labels))
y_all = np.array([level_names.index(l) for l in labels])
cov = np.loadtxt("covar.txt", skiprows=1, usecols=range(1, 6))
sex = np.loadtxt("sex.txt", skiprows=1, usecols=1)

variants = []
with gzip.open("mn.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        t = line.rstrip("\n").split("\t")
        fmt = t[8].split(":")
        g = np.empty(n)
        for i, cell in enumerate(t[9:]):
            parts = cell.split(":")
            if parts[0] in ("./.", "."):
                g[i] = np.nan
            elif "DS" in fmt:
                g[i] = round(float(parts[1]) * 16384) / 16384
            elif t[0] == "X" and "/" not in parts[0]:
                g[i] = 2.0 * int(parts[0])
            else:
                g[i] = parts[0].count("1")
        variants.append((t[0], t[2], g))

def fit(yv, X, start=None):
    return sm.MNLogit(yv, X).fit(method="newton", maxiter=200, disp=0, tol=1e-14, start_params=start)

with open(out_fname, "w") as fo:
    fo.write("#ID\tOBS_CT\tDF\tLRT\tSCORE\tWALD" + "".join("\tBETA_%s\tSE_%s" % (l, l) for l in level_names[1:]) + "\n")
    for chrom, vid, g in variants:
        nm = ~np.isnan(g)
        yv = y_all[nm]
        present = np.unique(yv)
        yc = np.searchsorted(present, yv)
        cols = [np.ones(nm.sum())]
        z = []
        if use_covar:
            z.append(cov[nm])
        if chrom == "X":
            z.append(sex[nm, None])
        if z:
            zz = np.column_stack(z)
            cols.append((zz - zz.mean(0)) / zz.std(0))
        X0 = np.column_stack(cols)
        gv = g[nm]
        X1 = np.column_stack([X0, gv - gv.mean()])
        sep = any((gv[yv == l].max() <= gv[yv != l].min()) or (gv[yv == l].min() >= gv[yv != l].max()) for l in present)
        r0 = fit(yc, X0)
        m1 = sm.MNLogit(yc, X1)
        start = np.vstack([r0.params, np.zeros((1, r0.params.shape[1]))]).ravel(order="F")
        sc = m1.score(start)
        score = sc @ np.linalg.solve(-m1.hessian(start), sc)
        df = len(present) - 1
        if sep:
            fo.write("%s\t%d\t%d\tSEP\t%.10g\tSEP" % (vid, nm.sum(), df, score) + "\tNA\tNA" * (len(level_names) - 1) + "\n")
            continue
        r1 = m1.fit(method="newton", maxiter=200, disp=0, tol=1e-14, start_params=start)
        lr = 2 * (r1.llf - r0.llf)
        k = X1.shape[1]
        idx = [j * k + (k - 1) for j in range(df)]
        b = np.asarray(r1.params)[-1]
        V = np.asarray(r1.cov_params())[np.ix_(idx, idx)]
        wald = b @ np.linalg.solve(V, b)
        se = np.sqrt(np.diag(V))
        coef_cells = []
        for lidx in range(1, len(level_names)):
            if present[0] == 0 and lidx in present:
                j = list(present).index(lidx) - 1
                coef_cells.append("%.10g\t%.10g" % (b[j], se[j]))
            else:
                coef_cells.append("NA\tNA")
        fo.write("%s\t%d\t%d\t%.10g\t%.10g\t%.10g\t%s\n" % (vid, nm.sum(), df, lr, score, wald, "\t".join(coef_cells)))
