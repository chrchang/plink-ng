#!/usr/bin/env python3
"""Recompute plink2's WEIGHTED_Z from the study files it was given.

METAL's weighted Z is sum_i sqrt(N_i) z_i / sqrt(sum_i N_i), with each study's
z the signed normal deviate of its p-value.  This derives it from the .assoc
files directly, with an inverse normal accurate to about 1e-15, so it pins
plink2's arithmetic without going through PLINK 1.9's approximation of the
same function.

Deliberately numpy-free: the test runners do not have numpy.
"""
import math
import sys


def normal_isf(p):
    """The x with P(X > x) = p, for a standard normal.

    Acklam's rational approximation plus one Halley step, which takes the
    relative error to around 1e-15.
    """
    if p <= 0.0:
        return math.inf
    if p >= 1.0:
        return -math.inf
    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00)
    p_low = 0.02425
    if p < p_low:
        u = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u
             + c[5]) / ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1.0)
    elif p <= 1.0 - p_low:
        u = p - 0.5
        r = u * u
        x = ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r
              + a[5]) * u
             / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r
                + 1.0))
    else:
        u = math.sqrt(-2.0 * math.log1p(-p))
        x = -(((((c[0] * u + c[1]) * u + c[2]) * u + c[3]) * u + c[4]) * u
              + c[5]) / ((((d[0] * u + d[1]) * u + d[2]) * u + d[3]) * u + 1.0)
    # Halley refinement on the lower tail, then flip to the upper one.
    e = 0.5 * math.erfc(-x / math.sqrt(2.0)) - p
    u = e * math.sqrt(2.0 * math.pi) * math.exp(x * x / 2.0)
    x = x - u / (1.0 + x * u / 2.0)
    return -x


def read_assoc(path):
    with open(path) as f:
        header = f.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        effect = 'OR' if 'OR' in idx else 'BETA'
        null_effect = 1.0 if effect == 'OR' else 0.0
        rows = {}
        for line in f:
            fields = line.split()
            if 'TEST' in idx and fields[idx['TEST']] != 'ADD':
                continue
            try:
                eff = float(fields[idx[effect]])
                p = float(fields[idx['P']])
                n = float(fields[idx['NMISS']])
            except (ValueError, KeyError):
                continue
            rows[fields[idx['SNP']]] = (eff - null_effect, p, n,
                                        fields[idx['A1']])
    return rows


def main():
    meta_path = sys.argv[1]
    studies = [read_assoc(p) for p in sys.argv[2:]]
    with open(meta_path) as f:
        header = f.readline().split()
        i_id = header.index('ID') if 'ID' in header else header.index('SNP')
        i_a1 = header.index('A1')
        i_n = header.index('N')
        i_wz = header.index('WEIGHTED_Z')
        rows = [line.split() for line in f if line.strip()]

    checked = 0
    worst = 0.0
    for fields in rows:
        vid = fields[i_id]
        if fields[i_wz] == 'NA':
            continue
        # Only the variants every study contributed to, and where plink2 says
        # so: anything else depends on which rows it dropped, which the other
        # checks cover.
        present = [s[vid] for s in studies if vid in s]
        if len(present) != len(studies):
            continue
        if int(fields[i_n]) != len(studies):
            continue
        num = 0.0
        den = 0.0
        skip = False
        for signed_effect, p, n, a1 in present:
            if not (0.0 < p <= 1.0) or signed_effect == 0.0:
                skip = True
                break
            z = normal_isf(p / 2.0)
            # plink2 reports the effect of the A1 it settled on, so a study
            # whose A1 is the other allele contributes with the sign flipped.
            if (signed_effect < 0.0) != (a1 != fields[i_a1]):
                z = -z
            num += math.sqrt(n) * z
            den += n
        if skip or den <= 0.0:
            continue
        want = num / math.sqrt(den)
        got = float(fields[i_wz])
        # plink2 prints about six significant digits.
        scale = max(abs(want), 1.0)
        rel = abs(got - want) / scale
        worst = max(worst, rel)
        if rel > 1e-5:
            sys.exit('WEIGHTED_Z on %s: plink2 says %s, recomputation says '
                     '%.9g (relative %.2e)' % (vid, fields[i_wz], want, rel))
        checked += 1
    if checked < 10:
        sys.exit('only %d variants checked; the comparison proves nothing'
                 % checked)
    print('%d WEIGHTED_Z values match the recomputation, worst relative '
          'difference %.2e' % (checked, worst))


if __name__ == '__main__':
    main()
