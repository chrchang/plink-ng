#!/usr/bin/env python3
"""Verify that --freq counts agrees with --geno-counts on chrX.

Both are computed from the same hardcalls, so for a biallelic variant

  ALT_CTS == HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS + HAP_ALT_CTS
  OBS_CT  == 2 * (HOM_REF_CT + HET_REF_ALT_CTS + TWO_ALT_GENO_CTS)
             + HAP_REF_CT + HAP_ALT_CTS

must hold exactly, with het haploid calls landing in MISSING_CT and thus in
neither side.  ALT_CTS must also be an integer.
"""
import sys


def load(fname):
    rows = {}
    with open(fname) as f:
        header = f.readline().split()
        id_idx = header.index('ID')
        for line in f:
            fields = line.split()
            rows[fields[id_idx]] = dict(zip(header, fields))
    return rows


acount = load(sys.argv[1])
gcount = load(sys.argv[2])
assert len(acount) == len(gcount), 'variant count mismatch'
assert acount, 'no variants'

problem_ct = 0
for variant_id, arow in acount.items():
    grow = gcount[variant_id]
    hom_ref = int(grow['HOM_REF_CT'])
    het = int(grow['HET_REF_ALT_CTS'])
    two_alt = int(grow['TWO_ALT_GENO_CTS'])
    hap_ref = int(grow['HAP_REF_CT'])
    hap_alt = int(grow['HAP_ALT_CTS'])
    expected_alt = het + 2 * two_alt + hap_alt
    expected_obs = 2 * (hom_ref + het + two_alt) + hap_ref + hap_alt
    alt_str = arow['ALT_CTS']
    obs_str = arow['OBS_CT']
    if ('.' in alt_str) or ('.' in obs_str):
        print('%s: non-integer --freq counts output: ALT_CTS %s, OBS_CT %s' %
              (variant_id, alt_str, obs_str))
        problem_ct += 1
        continue
    if (int(alt_str) != expected_alt) or (int(obs_str) != expected_obs):
        print('%s: --freq counts says ALT_CTS %s / OBS_CT %s, '
              '--geno-counts implies %u / %u' %
              (variant_id, alt_str, obs_str, expected_alt, expected_obs))
        problem_ct += 1

if problem_ct:
    sys.exit('%u variant(s) inconsistent' % (problem_ct,))
