#!/usr/bin/env python3
"""Checks a --pmerge result against the documented conflict rules.

Rather than freezing a transcript of today's output, this states the rules
themselves and applies them to whatever the fixtures import to:

  identical    a sample that is the same in both inputs comes out unchanged,
               in every mode, whether or not other samples conflict
  one-sided    where only one input has a value, that value wins
  nm-match     where both have values and they differ, the result is missing
  nm-first     where both have values, the first input's wins

A sample's value is its whole call, hardcall and dosage together: plink treats
them as a unit, since the hardcall is derived from the dosage.  Phase is
compared loosely, because a merge may drop it and the rules say nothing about
it.

Usage: check_merge.py <a body> <b body> <merged body> <nm-match|nm-first>
"""
import sys

a_f, b_f, m_f, mode = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]


def load(path):
    rows = {}
    header = None
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            header = f
            continue
        rows[f[2]] = dict(zip(header[9:], f[9:]))
    return rows, header[9:]


def norm(call):
    """Phase-insensitive form, with an absent dosage spelled the same way."""
    parts = call.replace('|', '/').split(':')
    gt = parts[0]
    ds = parts[1] if len(parts) > 1 else '.'
    if ds == '':
        ds = '.'
    return gt, ds


def is_missing(call):
    gt, ds = norm(call)
    return gt in ('./.', '.') and ds == '.'


a, samples = load(a_f)
b, _ = load(b_f)
m, _ = load(m_f)

bad = 0
for vid in sorted(set(a) & set(b)):
    for s in samples:
        av, bv, mv = a[vid][s], b[vid][s], m[vid][s]
        na, nb, nm_ = norm(av), norm(bv), norm(mv)
        if na == nb:
            want, why = na, 'identical in both inputs'
        elif is_missing(av):
            want, why = nb, 'only the second input has a value'
        elif is_missing(bv):
            want, why = na, 'only the first input has a value'
        elif mode == 'nm-match':
            want, why = ('./.', '.'), 'values differ under nm-match'
        else:
            want, why = na, 'first nonmissing value under nm-first'
        if nm_ != want:
            print('%s %s: %s, expected %s, got %s (a=%s b=%s)' % (
                vid, s, why, want, nm_, av, bv))
            bad += 1
if bad:
    sys.exit(1)
print('%s: %d variants x %d samples follow the rules' % (mode, len(set(a) & set(b)), len(samples)))
