#!/usr/bin/env python3
"""--glm numerical accuracy harness.

Usage:
  ./build.sh                      # builds ./plink2_fullprec
  ./harness.py [--plink2 PATH] [--only SUBSTRING ...] [--keep]

Generates small adversarial datasets (offset phenotypes, uncentered
covariates, high-VIF genotypes, missing calls, ...), runs --glm on them with a
plink2 binary that prints 17 significant digits, and compares BETA/SE against
a high-precision reference computed here:
  * linear regression: exact rational arithmetic, starting from the exact
    double values plink2 parses;
  * logistic and Firth regression: Newton iterations in 60-digit decimal
    arithmetic, run to full convergence.

Reported per scenario, separately for ADD rows and covariate/intercept rows
(max over variants):
  dB/SE   |beta - beta_ref| / SE_ref   (error in units of the standard error;
                                         this is what perturbs Z/T and P)
  relB    |beta - beta_ref| / |beta_ref|
  relSE   |SE - SE_ref| / SE_ref
All are printed as log10 values; -16 is roughly perfect double precision.
Logistic/Firth beta errors also include plink2's convergence tolerance, which
is not a numerical-linear-algebra issue.

Pure Python (no numpy) so it runs anywhere.
"""

import argparse
import math
import os
import random
import shutil
import subprocess
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 60
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Reference implementations
# ---------------------------------------------------------------------------

def gauss_jordan_inverse(mat, zero, one):
    """Inverse of a symmetric positive definite k x k matrix (list of lists)
    in exact/high-precision arithmetic."""
    k = len(mat)
    aug = [list(mat[i]) + [one if j == i else zero for j in range(k)] for i in range(k)]
    for col in range(k):
        piv = max(range(col, k), key=lambda r: abs(aug[r][col]))
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        if pv == 0:
            raise ZeroDivisionError("singular")
        inv_pv = one / pv
        aug[col] = [v * inv_pv for v in aug[col]]
        for r in range(k):
            if r != col and aug[r][col] != 0:
                f = aug[r][col]
                rowc = aug[col]
                aug[r] = [a - f * b for a, b in zip(aug[r], rowc)]
    return [row[k:] for row in aug]


def to_scaled_ints(col):
    """Exact representation of a list of doubles as integers / 2^shift."""
    ratios = [float(v).as_integer_ratio() for v in col]
    shift = max(d.bit_length() - 1 for _, d in ratios)
    return [n << (shift - (d.bit_length() - 1)) for n, d in ratios], shift


def linear_reference(cols, y):
    """cols: predictor columns (intercept first) as lists of doubles.
    Returns (betas, ses) as floats, computed exactly."""
    n = len(y)
    k = len(cols)
    icols = [to_scaled_ints(c) for c in cols]
    iy, yshift = to_scaled_ints(y)
    xtx = [[None] * k for _ in range(k)]
    for a in range(k):
        ca, sa = icols[a]
        for b in range(a + 1):
            cb, sb = icols[b]
            xtx[a][b] = xtx[b][a] = Fraction(sum(p * q for p, q in zip(ca, cb)), 1 << (sa + sb))
    xty = []
    for a in range(k):
        ca, sa = icols[a]
        xty.append(Fraction(sum(p * q for p, q in zip(ca, iy)), 1 << (sa + yshift)))
    yty = Fraction(sum(q * q for q in iy), 1 << (2 * yshift))
    inv = gauss_jordan_inverse(xtx, Fraction(0), Fraction(1))
    beta = [sum(inv[a][b] * xty[b] for b in range(k)) for a in range(k)]
    rss = yty - sum(b * t for b, t in zip(beta, xty))
    sigma2 = rss / (n - k)
    ses = []
    for a in range(k):
        v = sigma2 * inv[a][a]
        ses.append(float((Decimal(v.numerator) / Decimal(v.denominator)).sqrt()))
    return [float(b) for b in beta], ses


def logistic_reference(cols, y, firth):
    """Newton (or Firth-penalized Newton) in high-precision decimal arithmetic.
    y: list of 0/1.  Returns (betas, ses) as floats."""
    n = len(y)
    k = len(cols)
    X = [[Decimal(float(cols[j][i])) for j in range(k)] for i in range(n)]
    Y = [Decimal(v) for v in y]
    zero = Decimal(0)
    one = Decimal(1)
    half = Decimal("0.5")
    beta = [zero] * k
    for it in range(200):
        eta = [sum(xi[j] * beta[j] for j in range(k)) for xi in X]
        p = [one / (one + (-e).exp()) for e in eta]
        v = [pi * (one - pi) for pi in p]
        info = [[zero] * k for _ in range(k)]
        for i in range(n):
            xi = X[i]
            vi = v[i]
            for a in range(k):
                t = xi[a] * vi
                row = info[a]
                for b in range(a + 1):
                    row[b] += t * xi[b]
        for a in range(k):
            for b in range(a):
                info[b][a] = info[a][b]
        inv = gauss_jordan_inverse(info, zero, one)
        resid = [Y[i] - p[i] for i in range(n)]
        if firth:
            for i in range(n):
                xi = X[i]
                q = sum(xi[a] * sum(inv[a][b] * xi[b] for b in range(k)) for a in range(k))
                h = v[i] * q
                resid[i] += h * (half - p[i])
        score = [sum(X[i][a] * resid[i] for i in range(n)) for a in range(k)]
        step = [sum(inv[a][b] * score[b] for b in range(k)) for a in range(k)]
        beta = [b + s for b, s in zip(beta, step)]
        if max(abs(s) for s in step) < Decimal("1e-40"):
            break
    else:
        raise RuntimeError("reference Newton iterations did not converge")
    # information matrix at the final beta.  For Firth, plink2 (following
    # its logistf port) reports the inverse of X' diag(v * (1 + h)) X rather
    # than of X' diag(v) X, so we do the same.
    eta = [sum(xi[j] * beta[j] for j in range(k)) for xi in X]
    p = [one / (one + (-e).exp()) for e in eta]
    v = [pi * (one - pi) for pi in p]
    info = [[zero] * k for _ in range(k)]
    for i in range(n):
        for a in range(k):
            for b in range(k):
                info[a][b] += X[i][a] * v[i] * X[i][b]
    inv = gauss_jordan_inverse(info, zero, one)
    if firth:
        w = []
        for i in range(n):
            xi = X[i]
            q = sum(xi[a] * sum(inv[a][b] * xi[b] for b in range(k)) for a in range(k))
            w.append(v[i] * (one + v[i] * q))
        info = [[sum(X[i][a] * w[i] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
        inv = gauss_jordan_inverse(info, zero, one)
    return [float(b) for b in beta], [float(inv[a][a].sqrt()) for a in range(k)]


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------

class Scenario:
    def __init__(self, name, kind, n, variant_ct=5, seed=1, y_offset=0.0,
                 covar_means=(0.0, 0.0, 0.0), covar_sds=(1.0, 1.0, 1.0),
                 vif_link=False, missing_rate=0.0, glm_mods=(), y_scale=1.0):
        self.name = name
        self.kind = kind  # "linear", "logistic", "firth"
        self.n = n
        self.variant_ct = variant_ct
        self.seed = seed
        self.y_offset = y_offset
        self.y_scale = y_scale
        self.covar_means = covar_means
        self.covar_sds = covar_sds
        self.vif_link = vif_link
        self.missing_rate = missing_rate
        self.glm_mods = list(glm_mods)

    def generate(self):
        rng = random.Random(self.seed)
        n = self.n
        genos = []
        base = None
        for v in range(self.variant_ct):
            maf = 0.05 + 0.4 * v / max(1, self.variant_ct - 1)
            if self.vif_link and base is not None:
                # LD with the first variant: copy, then rerandomize ~10%
                g = [b if rng.random() > 0.1 else (rng.random() < maf) + (rng.random() < maf) for b in base]
            else:
                g = [(rng.random() < maf) + (rng.random() < maf) for _ in range(n)]
            if base is None:
                base = g
            gm = [None if rng.random() < self.missing_rate else x for x in g]
            genos.append(gm)
        covars = []
        for c, (m, sd) in enumerate(zip(self.covar_means, self.covar_sds)):
            if self.vif_link and c == 0:
                # strongly correlated with the first variant's genotypes
                col = [m + sd * (base[i] - 0.5 + rng.gauss(0, 0.25)) for i in range(n)]
            else:
                col = [m + sd * rng.gauss(0, 1) for _ in range(n)]
            covars.append([float(repr(x)) for x in col])
        lin = []
        for i in range(n):
            z = 0.15 * base[i]
            for c, (m, sd) in enumerate(zip(self.covar_means, self.covar_sds)):
                z += 0.3 * (covars[c][i] - m) / sd
            lin.append(z)
        if self.kind == "linear":
            pheno = [float(repr(self.y_offset + self.y_scale * (z + rng.gauss(0, 1)))) for z in lin]
        else:
            pheno = [1 if rng.random() < 1.0 / (1.0 + math.exp(-(z - 0.5))) else 0 for z in lin]
        self.genos = genos
        self.covars = covars
        self.pheno = pheno


def fmt_num(x):
    return "%d" % x if float(x).is_integer() else repr(x)


def scenarios():
    unc = (100.0, 100.0, -100.0)  # |mean|/sd = 100, just below the rescale trigger
    out = [
        Scenario("lin_baseline", "linear", 2000),
        Scenario("lin_offset_1e4", "linear", 2000, y_offset=1e4),
        Scenario("lin_offset_1e7", "linear", 2000, y_offset=1e7),
        Scenario("lin_uncentered_covars", "linear", 2000, covar_means=unc),
        Scenario("lin_high_vif", "linear", 2000, vif_link=True),
        Scenario("lin_combo", "linear", 2000, y_offset=1e4, covar_means=unc, vif_link=True, missing_rate=0.05),
        Scenario("lin_combo_intercept", "linear", 2000, y_offset=1e4, covar_means=unc, vif_link=True, missing_rate=0.05, glm_mods=["intercept"]),
        Scenario("lin_bign_offset", "linear", 50000, variant_ct=3, y_offset=1e3),
        Scenario("log_baseline", "logistic", 2000, glm_mods=["no-firth"]),
        Scenario("log_uncentered_covars", "logistic", 2000, covar_means=unc, glm_mods=["no-firth"]),
        Scenario("log_combo", "logistic", 2000, covar_means=unc, vif_link=True, missing_rate=0.05, glm_mods=["no-firth"]),
        Scenario("firth_baseline", "firth", 1000, variant_ct=3, glm_mods=["firth"]),
        Scenario("firth_uncentered_covars", "firth", 1000, variant_ct=3, covar_means=unc, glm_mods=["firth"]),
        Scenario("logsp_baseline", "logistic", 4000, variant_ct=3, glm_mods=["no-firth", "single-prec-cc"]),
        Scenario("logsp_uncentered_covars", "logistic", 4000, variant_ct=3, covar_means=unc, glm_mods=["no-firth", "single-prec-cc"]),
    ]
    return out


# ---------------------------------------------------------------------------
# Running plink2
# ---------------------------------------------------------------------------

def write_inputs(sc, wdir):
    n = sc.n
    ids = ["s%d" % i for i in range(n)]
    with open(os.path.join(wdir, "input.vcf"), "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=10000000>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
        for v, g in enumerate(sc.genos):
            gts = ["./." if x is None else ("0/0", "0/1", "1/1")[x] for x in g]
            f.write("1\t%d\tv%d\tA\tG\t.\t.\t.\tGT\t%s\n" % (1000 + 100 * v, v, "\t".join(gts)))
    with open(os.path.join(wdir, "input.psam"), "w") as f:
        f.write("#IID\tSEX\tY\t" + "\t".join("C%d" % (c + 1) for c in range(len(sc.covars))) + "\n")
        for i in range(n):
            y = fmt_num(sc.pheno[i]) if sc.kind == "linear" else str(sc.pheno[i] + 1)
            f.write("%s\tNA\t%s\t%s\n" % (ids[i], y, "\t".join(fmt_num(c[i]) for c in sc.covars)))


def run_plink2(plink2, sc, wdir):
    def run(args):
        r = subprocess.run([plink2] + args, cwd=wdir, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode:
            sys.stderr.write(r.stderr)
            raise RuntimeError("plink2 failed: " + " ".join(args))
    run(["--vcf", "input.vcf", "--psam", "input.psam", "--make-pgen", "--out", "data"])
    # plink2_fullprec writes .psam phenotypes with 17 significant digits, so
    # this recovers the exact doubles plink2 parsed (its decimal parser can be
    # 1-2 ulps off; that input error must not be charged to the regression).
    # --glm below parses the same text with the same parser.
    with open(os.path.join(wdir, "data.psam")) as f:
        header = f.readline().rstrip("\n").split("\t")
        hidx = {h: i for i, h in enumerate(header)}
        parsed = [line.rstrip("\n").split("\t") for line in f]
    assert len(parsed) == sc.n
    for c in range(len(sc.covars)):
        sc.covars[c] = [float(t[hidx["C%d" % (c + 1)]]) for t in parsed]
    if sc.kind == "linear":
        sc.pheno = [float(t[hidx["Y"]]) for t in parsed]
    covar_names = ",".join("C%d" % (c + 1) for c in range(len(sc.covars)))
    run(["--pfile", "data", "--covar", "input.psam", "--pheno", "input.psam", "--covar-name", covar_names, "--pheno-name", "Y",
         "--glm", "cols=+beta"] + sc.glm_mods + ["--out", "out"])
    fn = [f for f in os.listdir(wdir) if f.startswith("out.Y.glm.")]
    assert len(fn) == 1, fn
    rows = {}
    with open(os.path.join(wdir, fn[0])) as f:
        header = f.readline().rstrip("\n").split("\t")
        hidx = {h: i for i, h in enumerate(header)}
        se_col = "SE" if "SE" in hidx else "LOG(OR)_SE"
        for line in f:
            t = line.rstrip("\n").split("\t")
            beta = t[hidx["BETA"]]
            se = t[hidx[se_col]]
            if beta == "NA":
                rows[(t[hidx["ID"]], t[hidx["TEST"]])] = None
            else:
                rows[(t[hidx["ID"]], t[hidx["TEST"]])] = (float(beta), float(se))
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def log10_or_floor(x):
    return -99.0 if x == 0 else math.log10(x)


def evaluate(sc, rows):
    covar_names = ["C%d" % (c + 1) for c in range(len(sc.covars))]
    stats = {"ADD": [-99.0] * 3, "covar": [-99.0] * 3}
    na_ct = 0
    for v, g in enumerate(sc.genos):
        keep = [i for i in range(sc.n) if g[i] is not None]
        cols = [[1.0] * len(keep), [float(g[i]) for i in keep]] + [[c[i] for i in keep] for c in sc.covars]
        y = [sc.pheno[i] for i in keep]
        if sc.kind == "linear":
            ref_b, ref_se = linear_reference(cols, y)
        else:
            ref_b, ref_se = logistic_reference(cols, y, sc.kind == "firth")
        tests = ["INTERCEPT", "ADD"] + covar_names
        for j, test in enumerate(tests):
            key = ("v%d" % v, test)
            if key not in rows:
                continue
            got = rows[key]
            if got is None:
                na_ct += 1
                continue
            b, se = got
            d = abs(b - ref_b[j])
            cur = [d / ref_se[j], d / abs(ref_b[j]) if ref_b[j] else 0.0, abs(se - ref_se[j]) / ref_se[j]]
            slot = stats["ADD" if test == "ADD" else "covar"]
            for m in range(3):
                slot[m] = max(slot[m], log10_or_floor(cur[m]))
    return stats, na_ct


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plink2", default=os.path.join(HERE, "plink2_fullprec"))
    ap.add_argument("--only", nargs="*", default=None, help="only run scenarios whose names contain one of these")
    ap.add_argument("--keep", action="store_true", help="keep work directories")
    ap.add_argument("--workdir", default=os.path.join(HERE, "work"), help="scratch directory (default ./work)")
    args = ap.parse_args()
    plink2 = os.path.abspath(args.plink2)
    root = os.path.abspath(args.workdir)
    os.makedirs(root, exist_ok=True)
    print("%-26s | %-22s | %-22s | %s" % ("scenario", "ADD dB/SE relB relSE", "covar dB/SE relB relSE", "NA"))
    for sc in scenarios():
        if args.only and not any(s in sc.name for s in args.only):
            continue
        wdir = os.path.join(root, sc.name)
        shutil.rmtree(wdir, ignore_errors=True)
        os.makedirs(wdir)
        sc.generate()
        write_inputs(sc, wdir)
        rows = run_plink2(plink2, sc, wdir)
        stats, na_ct = evaluate(sc, rows)
        fmt = lambda s: "%6.1f %6.1f %6.1f" % tuple(s)
        print("%-26s | %-22s | %-22s | %d" % (sc.name, fmt(stats["ADD"]), fmt(stats["covar"]), na_ct), flush=True)
        if not args.keep:
            shutil.rmtree(wdir)


if __name__ == "__main__":
    main()
