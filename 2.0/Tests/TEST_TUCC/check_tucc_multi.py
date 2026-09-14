#!/usr/bin/env python3
"""Compare a --tucc fileset's multiallelic genotypes against expected pairs."""
import sys


def main():
    vcf_path, expected_path = sys.argv[1], sys.argv[2]
    ids = []
    genotypes = {}
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'):
                ids = fields[9:]
                continue
            vid = fields[2]
            calls = []
            for field in fields[9:]:
                gt = field.split(':', 1)[0]
                if gt == '.':
                    gt = './.'
                alleles = gt.replace('|', '/').split('/')
                if '.' in alleles:
                    calls.append('./.')
                else:
                    lo, hi = sorted(int(a) for a in alleles)
                    calls.append('%d/%d' % (lo, hi))
            genotypes[vid] = calls
    checked = 0
    nonmissing = 0
    het_nonref = 0
    with open(expected_path) as f:
        for line in f:
            fields = line.split()
            vid = fields[0]
            expected = fields[1:]
            got = genotypes.get(vid)
            if got is None:
                sys.exit('%s: missing from the exported VCF' % vid)
            if len(expected) != len(got):
                sys.exit('%s: %d expected values for %d samples' %
                         (vid, len(expected), len(got)))
            for idx, want in enumerate(expected):
                if got[idx] != want:
                    sys.exit('%s / %s: expected %s, got %s' %
                             (vid, ids[idx], want, got[idx]))
                checked += 1
                if want != './.':
                    nonmissing += 1
                    lo, hi = want.split('/')
                    if lo != '0' and hi != '0' and lo != hi:
                        het_nonref += 1
    if not nonmissing:
        sys.exit('every expected genotype is missing; the test proves nothing')
    if not het_nonref:
        sys.exit('no genotype pairs two distinct ALT alleles; the test would '
                 'pass without multiallelic support')
    print('%d multiallelic genotypes verified (%d nonmissing, %d ALT/ALT het)'
          % (checked, nonmissing, het_nonref))


if __name__ == '__main__':
    main()
