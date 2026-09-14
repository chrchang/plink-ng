#!/usr/bin/env python3
"""An independent check of `ldsc --munge`, from the raw file it was given.

Re-derives which variants should have survived and what their Z and N should
be, then compares that against the .sumstats file.  The Z magnitude is checked
by going the other way, p = erfc(|Z| / sqrt(2)), so this does not need an
inverse normal of its own.

Deliberately numpy-free: the test runners do not have numpy.
"""
import math
import sys

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
ALIASES = {
    'SNP': 'SNP', 'MARKERNAME': 'SNP', 'SNPID': 'SNP', 'RS': 'SNP',
    'RSID': 'SNP', 'RS_NUMBER': 'SNP', 'RS_NUMBERS': 'SNP',
    'P': 'P', 'PVALUE': 'P', 'P_VALUE': 'P', 'PVAL': 'P', 'P_VAL': 'P',
    'GC_PVALUE': 'P',
    'A1': 'A1', 'ALLELE1': 'A1', 'ALLELE_1': 'A1', 'EFFECT_ALLELE': 'A1',
    'REFERENCE_ALLELE': 'A1', 'INC_ALLELE': 'A1', 'EA': 'A1',
    'A2': 'A2', 'ALLELE2': 'A2', 'ALLELE_2': 'A2', 'OTHER_ALLELE': 'A2',
    'NON_EFFECT_ALLELE': 'A2', 'DEC_ALLELE': 'A2', 'NEA': 'A2',
    'N': 'N', 'WEIGHT': 'N',
    'NCASE': 'N_CAS', 'CASES_N': 'N_CAS', 'N_CASE': 'N_CAS',
    'N_CASES': 'N_CAS', 'N_CAS': 'N_CAS',
    'N_CONTROLS': 'N_CON', 'N_CON': 'N_CON', 'NCONTROL': 'N_CON',
    'CONTROLS_N': 'N_CON', 'N_CONTROL': 'N_CON',
    'ZSCORE': 'Z', 'Z_SCORE': 'Z', 'GC_ZSCORE': 'Z', 'Z': 'Z',
    'OR': 'OR', 'B': 'BETA', 'BETA': 'BETA', 'EFFECTS': 'BETA',
    'EFFECT': 'BETA', 'LOG_ODDS': 'LOG_ODDS',
    'INFO': 'INFO',
    'EAF': 'FRQ', 'FRQ': 'FRQ', 'MAF': 'FRQ', 'FRQ_U': 'FRQ', 'F_U': 'FRQ',
}
SIGNED_NULLS = {'Z': 0.0, 'OR': 1.0, 'BETA': 0.0, 'LOG_ODDS': 0.0}
MISSING = {'.', '?', 'NA', 'na', ''}


def clean(name):
    return name.upper().replace('-', '_').replace('.', '_')


def valid_snp(a1, a2):
    if a1 not in COMPLEMENT or a2 not in COMPLEMENT or a1 == a2:
        return False
    return COMPLEMENT[a1] != a2


def read_raw(path, daner):
    with open(path) as f:
        header = [clean(h) for h in f.readline().split()]
        rows = [line.split() for line in f if line.strip()]
    fields = {}
    signed = None
    daner_n = None
    for idx, name in enumerate(header):
        if daner and name.startswith('FRQ_A_'):
            daner_n = (float(name[6:]), daner_n[1] if daner_n else None)
            continue
        if daner and name.startswith('FRQ_U_'):
            n_con = float(name[6:])
            daner_n = (daner_n[0] if daner_n else None, n_con)
            fields['FRQ'] = idx
            continue
        internal = ALIASES.get(name)
        if internal is None:
            continue
        if internal in SIGNED_NULLS:
            signed = (idx, internal)
            continue
        if daner and internal in ('N', 'N_CAS', 'N_CON'):
            continue
        fields[internal] = idx
    return header, rows, fields, signed, daner_n


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw', required=True)
    ap.add_argument('--munged', required=True)
    ap.add_argument('--merge-alleles', default=None)
    ap.add_argument('--info-min', type=float, default=0.9)
    ap.add_argument('--maf-min', type=float, default=0.01)
    ap.add_argument('--a1-inc', action='store_true')
    ap.add_argument('--daner', action='store_true')
    ap.add_argument('--keep-maf', action='store_true')
    ap.add_argument('--N', type=float, default=None)
    args = ap.parse_args()

    header, rows, fields, signed, daner_n = read_raw(args.raw, args.daner)
    merge_ids = None
    merge_alleles = {}
    if args.merge_alleles:
        merge_ids = []
        with open(args.merge_alleles) as f:
            mh = [clean(h) for h in f.readline().split()]
            i_snp, i_a1, i_a2 = mh.index('SNP'), mh.index('A1'), mh.index('A2')
            for line in f:
                fl = line.split()
                if fl[i_snp] in merge_alleles:
                    continue
                merge_ids.append(fl[i_snp])
                merge_alleles[fl[i_snp]] = (fl[i_a1].upper(), fl[i_a2].upper())

    kept = []
    seen = set()
    for fl in rows:
        def get(name):
            idx = fields.get(name)
            if idx is None or idx >= len(fl):
                return None
            val = fl[idx]
            return None if val in MISSING else val

        snp = get('SNP')
        p_raw = get('P')
        a1 = get('A1')
        a2 = get('A2')
        sign_raw = None if signed is None else (
            None if fl[signed[0]] in MISSING else fl[signed[0]])
        if snp is None or p_raw is None:
            continue
        if not args.a1_inc and sign_raw is None:
            continue
        if a1 is None or a2 is None:
            continue
        a1, a2 = a1.upper(), a2.upper()
        if merge_ids is not None and snp not in merge_alleles:
            continue
        info = get('INFO')
        if 'INFO' in fields and info is not None:
            if float(info) < args.info_min:
                continue
        frq = get('FRQ')
        if 'FRQ' in fields:
            if frq is None:
                continue
            frq_v = float(frq)
            if frq_v < 0.0 or frq_v > 1.0:
                continue
            if min(frq_v, 1.0 - frq_v) <= args.maf_min:
                continue
        p = float(p_raw)
        if not (0.0 < p <= 1.0):
            continue
        if not valid_snp(a1, a2):
            continue
        if snp in seen:
            continue
        seen.add(snp)
        n_cas = get('N_CAS')
        n_con = get('N_CON')
        n = get('N')
        kept.append({
            'SNP': snp, 'A1': a1, 'A2': a2, 'P': p,
            'sign': None if sign_raw is None else float(sign_raw),
            'N': None if n is None else float(n),
            'N_CAS': None if n_cas is None else float(n_cas),
            'N_CON': None if n_con is None else float(n_con),
            'FRQ': None if frq is None else float(frq),
        })

    # Effective sample size from case and control counts.
    if kept and kept[0]['N_CAS'] is not None and kept[0]['N_CON'] is not None:
        totals = [r['N_CAS'] + r['N_CON'] for r in kept]
        max_n = max(totals)
        fracs = [r['N_CAS'] / (r['N_CAS'] + r['N_CON'])
                 for r in kept if r['N_CAS'] + r['N_CON'] == max_n]
        ref_frac = sum(fracs) / len(fracs)
        for r in kept:
            tot = r['N_CAS'] + r['N_CON']
            r['N'] = tot * (r['N_CAS'] / tot) / ref_frac

    if args.daner:
        for r in kept:
            r['N'] = daner_n[0] + daner_n[1]
    elif args.N is not None:
        for r in kept:
            r['N'] = args.N

    # The sample-size floor: the 90th percentile over 1.5, unless every N is
    # the same (in which case nothing is dropped).
    if kept and kept[0]['N'] is not None and not args.daner and args.N is None:
        ns = sorted(r['N'] for r in kept)
        pos = 0.9 * (len(ns) - 1)
        lo = int(pos)
        frac = pos - lo
        q90 = ns[lo] * (1 - frac) + ns[min(lo + 1, len(ns) - 1)] * frac
        n_min = q90 / 1.5
        kept = [r for r in kept if r['N'] >= n_min]

    with open(args.munged) as f:
        out_header = f.readline().split()
        out_rows = [line.rstrip('\n').split('\t') for line in f if line.strip()]
    got = {}
    for fl in out_rows:
        d = dict(zip(out_header, fl))
        got[d['SNP']] = d

    if merge_ids is not None:
        if [r[0] for r in out_rows] != merge_ids:
            sys.exit('--merge-alleles output must list the merge file\'s '
                     'variants in its order')
        expected_ids = set(r['SNP'] for r in kept)
        for snp in merge_ids:
            present = got[snp]['Z'] != 'NA'
            want = snp in expected_ids
            if present != want:
                sys.exit('%s: %s in the output, expected %s'
                         % (snp, 'present' if present else 'missing',
                            'present' if want else 'missing'))
        kept = [r for r in kept if r['SNP'] in got and got[r['SNP']]['Z'] != 'NA']
    else:
        if len(kept) != len(out_rows):
            sys.exit('expected %d variants, got %d' % (len(kept), len(out_rows)))

    checked = 0
    for r in kept:
        d = got.get(r['SNP'])
        if d is None:
            sys.exit('%s: missing from the output' % r['SNP'])
        if 'A1' in d and (d['A1'] != r['A1'] or d['A2'] != r['A2']):
            sys.exit('%s: alleles %s/%s, expected %s/%s'
                     % (r['SNP'], d['A1'], d['A2'], r['A1'], r['A2']))
        z = float(d['Z'])
        # The p-value has to come back out of |Z|.
        p_back = math.erfc(abs(z) / math.sqrt(2.0))
        if abs(p_back - r['P']) > 2e-3 * max(r['P'], 1e-4) + 2e-3:
            sys.exit('%s: Z = %g implies p = %g, but the input had p = %g'
                     % (r['SNP'], z, p_back, r['P']))
        if r['sign'] is not None and z != 0.0:
            # The sign of Z has to follow the signed statistic.
            null = SIGNED_NULLS[signed[1]]
            if (r['sign'] < null) != (z < 0.0):
                sys.exit('%s: signed statistic %g (null %g) but Z = %g'
                         % (r['SNP'], r['sign'], null, z))
        if abs(float(d['N']) - r['N']) > 1e-3 + 1e-9 * abs(r['N']):
            sys.exit('%s: N = %s, expected %g' % (r['SNP'], d['N'], r['N']))
        if args.keep_maf and 'FRQ' in d:
            if abs(float(d['FRQ']) - r['FRQ']) > 1.1e-3:
                sys.exit('%s: FRQ = %s, expected %g'
                         % (r['SNP'], d['FRQ'], r['FRQ']))
        checked += 1
    if not checked:
        sys.exit('no variants checked; the comparison proves nothing')
    print('%d munged variants verified' % checked)


if __name__ == '__main__':
    main()
