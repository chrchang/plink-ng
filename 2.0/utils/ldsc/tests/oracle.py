#!/usr/bin/env python3
"""An independent implementation of the LD Score regressions, for testing.

This is a second implementation of the same estimator, in a different
language, written from the same description: the regression weights of
Bulik-Sullivan et al. (2015), two iterations of IRWLS, the two-step estimator,
and the block jackknife.  Its purpose is to disagree with ldsc.cc when
ldsc.cc has a coding error, so it deliberately shares no code with it.

Deliberately numpy-free: the test runners do not have numpy.
"""
import sys


# ***** linear algebra (at most 2 parameters, so this is all closed-form) *****

def solve(mat, vec):
    """Solves mat * out = vec for a 1x1 or 2x2 system."""
    p = len(vec)
    if p == 1:
        return [vec[0] / mat[0][0]]
    det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
    return [(vec[0] * mat[1][1] - vec[1] * mat[0][1]) / det,
            (vec[1] * mat[0][0] - vec[0] * mat[1][0]) / det]


def xtx_xty(x, y, w2, lo, hi):
    """Accumulates X'WX and X'Wy over rows [lo, hi)."""
    p = len(x[0])
    mat = [[0.0] * p for _ in range(p)]
    vec = [0.0] * p
    for i in range(lo, hi):
        for j in range(p):
            vec[j] += w2[i] * x[i][j] * y[i]
            for k in range(p):
                mat[j][k] += w2[i] * x[i][j] * x[i][k]
    return mat, vec


def wls(x, y, w, lo=0, hi=None):
    if hi is None:
        hi = len(y)
    w2 = [wi * wi for wi in w]
    mat, vec = xtx_xty(x, y, w2, lo, hi)
    return solve(mat, vec)


# ***** block jackknife *****

def separators(n, n_blocks):
    import math
    sep = [int(math.floor(i * n / n_blocks)) for i in range(n_blocks + 1)]
    sep[0] = 0
    sep[n_blocks] = n
    return sep


def mat_add(a, b, sign):
    return [[a[i][j] + sign * b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def jknife_from_delete(est, delete_values):
    """est/delete values -> jackknife estimate and covariance."""
    n_blocks = len(delete_values)
    p = len(est)
    pseudo = [[n_blocks * est[j] - (n_blocks - 1) * delete_values[b][j]
               for j in range(p)] for b in range(n_blocks)]
    mean = [sum(pseudo[b][j] for b in range(n_blocks)) / n_blocks
            for j in range(p)]
    cov = [[sum((pseudo[b][j] - mean[j]) * (pseudo[b][k] - mean[k])
                for b in range(n_blocks)) / ((n_blocks - 1) * n_blocks)
            for k in range(p)] for j in range(p)]
    se = [cov[j][j] ** 0.5 for j in range(p)]
    return mean, cov, se, pseudo


def lstsq_jknife(x, y, w, sep):
    """Weighted regression plus its block jackknife delete values."""
    n_blocks = len(sep) - 1
    p = len(x[0])
    w2 = [wi * wi for wi in w]
    blocks = [xtx_xty(x, y, w2, sep[b], sep[b + 1]) for b in range(n_blocks)]
    tot_mat = [[0.0] * p for _ in range(p)]
    tot_vec = [0.0] * p
    for mat, vec in blocks:
        tot_mat = mat_add(tot_mat, mat, 1)
        tot_vec = [tot_vec[j] + vec[j] for j in range(p)]
    est = solve(tot_mat, tot_vec)
    delete_values = []
    for mat, vec in blocks:
        del_mat = mat_add(tot_mat, mat, -1)
        del_vec = [tot_vec[j] - vec[j] for j in range(p)]
        delete_values.append(solve(del_mat, del_vec))
    return est, delete_values


# ***** weights *****

def clamp(val, lo, hi):
    return lo if val < lo else (hi if val > hi else val)


def hsq_weights(ld, w_ld, n, m, hsq, intercept, idxs):
    hsq = clamp(hsq, 0.0, 1.0)
    out = []
    for i in idxs:
        cur_ld = max(ld[i], 1.0)
        cur_w = max(w_ld[i], 1.0)
        c = hsq * n[i] / m
        denom = intercept + c * cur_ld
        out.append(1.0 / (2 * denom * denom * cur_w))
    return out


def gencov_weights(ld, w_ld, n1, n2, m, h1, h2, rho_g, intercept,
                   intercept_h1, intercept_h2, idxs):
    h1 = clamp(h1, 0.0, 1.0)
    h2 = clamp(h2, 0.0, 1.0)
    rho_g = clamp(rho_g, -1.0, 1.0)
    out = []
    for i in idxs:
        cur_ld = max(ld[i], 1.0)
        cur_w = max(w_ld[i], 1.0)
        a = n1[i] * h1 * cur_ld / m + intercept_h1
        b = n2[i] * h2 * cur_ld / m + intercept_h2
        c = (n1[i] * n2[i]) ** 0.5 * rho_g * cur_ld / m + intercept
        out.append(1.0 / ((a * b + c * c) * cur_w))
    return out


# ***** the regressions *****

class Reg(object):
    """Shared LD Score regression body."""

    def __init__(self, kind, y, ld, w_ld, n, m, n_blocks, intercept=None,
                 twostep=None, step1_filter=None, gencov_args=None):
        self.kind = kind
        self.gencov_args = gencov_args
        n_snp = len(y)
        self.nbar = sum(n) / n_snp
        null_intercept = 1.0 if kind == 'hsq' else 0.0
        agg_int = intercept if intercept is not None else null_intercept
        mean_y = sum(y) / n_snp
        mean_ldn = sum(ld[i] * n[i] for i in range(n_snp)) / n_snp
        tot_agg = m * (mean_y - agg_int) / mean_ldn
        all_idxs = list(range(n_snp))
        initial_w = self._weights(tot_agg, agg_int, ld, w_ld, n, m, all_idxs)
        x_scaled = [n[i] * ld[i] / self.nbar for i in range(n_snp)]

        if intercept is not None:
            yp = [yi - intercept for yi in y]
            design = [[xi] for xi in x_scaled]
            est, delete_values = self._irwls(
                design, yp, initial_w, ld, w_ld, n, m, all_idxs, intercept,
                separators(n_snp, n_blocks))
        elif twostep is None:
            design = [[xi, 1.0] for xi in x_scaled]
            est, delete_values = self._irwls(
                design, y, initial_w, ld, w_ld, n, m, all_idxs, None,
                separators(n_snp, n_blocks))
        else:
            step1_idxs = [i for i in range(n_snp) if step1_filter[i] < twostep]
            n1 = len(step1_idxs)
            design1 = [[x_scaled[i], 1.0] for i in step1_idxs]
            y1 = [y[i] for i in step1_idxs]
            w1 = [initial_w[i] for i in step1_idxs]
            sep1 = separators(n1, n_blocks)
            est1, delete1 = self._irwls(design1, y1, w1, ld, w_ld, n, m,
                                        step1_idxs, None, sep1)
            step1_int = est1[1]
            yp = [yi - step1_int for yi in y]
            design2 = [[xi] for xi in x_scaled]
            sep2 = [0] + [step1_idxs[sep1[b]] for b in range(1, n_blocks)] + [n_snp]
            est2, delete2 = self._irwls(design2, yp, initial_w, ld, w_ld, n, m,
                                        all_idxs, step1_int, sep2)
            c_num = sum(initial_w[i] * x_scaled[i] for i in range(n_snp))
            c_den = sum(initial_w[i] * x_scaled[i] ** 2 for i in range(n_snp))
            c = c_num / c_den
            est = [est2[0], step1_int]
            delete_values = [[delete2[b][0] - c * (delete1[b][1] - step1_int),
                              delete1[b][1]] for b in range(n_blocks)]

        _, cov, se, _ = jknife_from_delete(est, delete_values)
        self.est = est
        self.coef = est[0] / self.nbar
        self.tot = m * self.coef
        self.tot_se = (m * m * cov[0][0] / self.nbar ** 2) ** 0.5
        if intercept is not None:
            self.intercept = intercept
            self.intercept_se = None
        else:
            self.intercept = est[1]
            self.intercept_se = se[1]
        self.tot_delete_values = [d[0] * m / self.nbar for d in delete_values]
        self.n_blocks = len(delete_values)

    def _weights(self, param, intercept, ld, w_ld, n, m, idxs):
        if self.kind == 'hsq':
            return hsq_weights(ld, w_ld, n, m, param, intercept, idxs)
        h1, h2, int_h1, int_h2, n1, n2 = self.gencov_args
        return gencov_weights(ld, w_ld, n1, n2, m, h1, h2, param, intercept,
                              int_h1, int_h2, idxs)

    def _irwls(self, design, y, initial_w, ld, w_ld, n, m, idxs, intercept,
               sep):
        w = [wi ** 0.5 for wi in initial_w]
        null_intercept = 1.0 if self.kind == 'hsq' else 0.0
        for _ in range(2):
            coef = wls(design, y, w)
            param = m * coef[0] / self.nbar
            cur_int = intercept if intercept is not None else coef[1]
            if cur_int is None:
                cur_int = null_intercept
            new_w = self._weights(param, cur_int, ld, w_ld, n, m, idxs)
            w = [wi ** 0.5 for wi in new_w]
        return lstsq_jknife(design, y, w, sep)


def ratio_jknife(est, numer_delete, denom_delete):
    n_blocks = len(numer_delete)
    delete_values = [[numer_delete[b] / denom_delete[b]]
                     for b in range(n_blocks)]
    _, cov, se, _ = jknife_from_delete([est], delete_values)
    return se[0]


# ***** file reading (the same inner joins, in LD Score file order) *****

def read_table(path):
    with open(path) as f:
        header = f.readline().lstrip('#').split()
        rows = [line.split() for line in f if line.strip()]
    return header, rows


def read_ldscore(path):
    header, rows = read_table(path)
    id_col = header.index('SNP') if 'SNP' in header else header.index('ID')
    l2_col = header.index('L2')
    return [r[id_col] for r in rows], [float(r[l2_col]) for r in rows]


def read_sumstats(path):
    header, rows = read_table(path)
    id_col = header.index('SNP') if 'SNP' in header else header.index('ID')
    z_col = header.index('Z')
    n_col = header.index('N')
    return {r[id_col]: (float(r[z_col]), float(r[n_col])) for r in rows}


def merge(ref_ids, ref_l2, w_map, ss1, ss2=None):
    ld, w, z1, n1, z2, n2 = [], [], [], [], [], []
    for i, snp in enumerate(ref_ids):
        if snp not in ss1 or snp not in w_map:
            continue
        if ss2 is not None and snp not in ss2:
            continue
        ld.append(ref_l2[i])
        w.append(w_map[snp])
        z1.append(ss1[snp][0])
        n1.append(ss1[snp][1])
        if ss2 is not None:
            z2.append(ss2[snp][0])
            n2.append(ss2[snp][1])
    return ld, w, z1, n1, z2, n2


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['h2', 'rg'])
    ap.add_argument('--ref-ld', required=True)
    ap.add_argument('--w-ld', required=True)
    ap.add_argument('--sumstats', required=True,
                    help='one file, or two comma-separated for rg')
    ap.add_argument('--M', type=float, required=True)
    ap.add_argument('--n-blocks', type=int, default=200)
    ap.add_argument('--intercept-h2', type=float, default=None)
    ap.add_argument('--two-step', type=float, default=None)
    ap.add_argument('--no-two-step', action='store_true')
    args = ap.parse_args()

    ref_ids, ref_l2 = read_ldscore(args.ref_ld)
    w_ids, w_l2 = read_ldscore(args.w_ld)
    w_map = dict(zip(w_ids, w_l2))
    paths = args.sumstats.split(',')
    ss1 = read_sumstats(paths[0])
    ss2 = read_sumstats(paths[1]) if args.mode == 'rg' else None
    ld, w, z1, n1, z2, n2 = merge(ref_ids, ref_l2, w_map, ss1, ss2)
    n_snp = len(ld)
    n_blocks = min(args.n_blocks, n_snp)
    twostep = args.two_step
    if twostep is None and not args.no_two_step and args.intercept_h2 is None:
        twostep = 30.0
    if args.no_two_step or args.intercept_h2 is not None:
        twostep = None

    chisq1 = [zi * zi for zi in z1]
    hsq1 = Reg('hsq', chisq1, ld, w, n1, args.M, n_blocks,
               intercept=args.intercept_h2, twostep=twostep,
               step1_filter=chisq1)
    if args.mode == 'h2':
        mean_chisq = sum(chisq1) / n_snp
        print('h2 %.12g' % hsq1.tot)
        print('h2_se %.12g' % hsq1.tot_se)
        print('intercept %.12g' % hsq1.intercept)
        if hsq1.intercept_se is not None:
            print('intercept_se %.12g' % hsq1.intercept_se)
            if mean_chisq > 1:
                print('ratio %.12g' % ((hsq1.intercept - 1) / (mean_chisq - 1)))
                print('ratio_se %.12g' % (hsq1.intercept_se / (mean_chisq - 1)))
        print('mean_chisq %.12g' % mean_chisq)
        return

    chisq2 = [zi * zi for zi in z2]
    hsq2 = Reg('hsq', chisq2, ld, w, n2, args.M, n_blocks,
               intercept=args.intercept_h2, twostep=twostep,
               step1_filter=chisq2)
    y = [z1[i] * z2[i] for i in range(n_snp)]
    sqrt_n1n2 = [(n1[i] * n2[i]) ** 0.5 for i in range(n_snp)]
    step1_filter = [max(chisq1[i], chisq2[i]) for i in range(n_snp)]
    gencov = Reg('gencov', y, ld, w, sqrt_n1n2, args.M, n_blocks,
                 intercept=None, twostep=twostep, step1_filter=step1_filter,
                 gencov_args=(hsq1.tot, hsq2.tot, hsq1.intercept,
                              hsq2.intercept, n1, n2))
    rg = gencov.tot / (hsq1.tot * hsq2.tot) ** 0.5
    denom = [(hsq1.tot_delete_values[b] * hsq2.tot_delete_values[b]) ** 0.5
             for b in range(n_blocks)]
    rg_se = ratio_jknife(rg, gencov.tot_delete_values, denom)
    print('h2_p1 %.12g' % hsq1.tot)
    print('h2_obs %.12g' % hsq2.tot)
    print('gcov %.12g' % gencov.tot)
    print('gcov_se %.12g' % gencov.tot_se)
    print('gcov_int %.12g' % gencov.intercept)
    print('rg %.12g' % rg)
    print('se %.12g' % rg_se)
    print('z %.12g' % (rg / rg_se))


if __name__ == '__main__':
    main()
