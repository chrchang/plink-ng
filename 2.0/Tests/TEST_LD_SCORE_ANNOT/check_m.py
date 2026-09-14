#!/usr/bin/env python3
"""Check the .M and .M_5_50 variant counts against the genotypes."""
import sys

from brute_ldscore import read_geno


def main():
    raw_path, annot_path, m_path, m_5_50_path = sys.argv[1:5]
    ids, cols = read_geno(raw_path)
    with open(annot_path) as f:
        annot_names = f.readline().split()[1:]
        annots = {}
        for line in f:
            fields = line.split()
            annots[fields[0]] = [float(x) for x in fields[1:]]
    n_annot = len(annot_names)
    m = [0.0] * n_annot
    m_common = [0.0] * n_annot
    for j, vid in enumerate(ids):
        present = [v for v in cols[j] if v is not None]
        if not present:
            continue
        # --export A counts one allele, so its mean over 2 is that allele's
        # frequency.
        freq = sum(present) / (2.0 * len(present))
        maf = min(freq, 1.0 - freq)
        for c in range(n_annot):
            m[c] += annots[vid][c]
            if maf > 0.05:
                m_common[c] += annots[vid][c]
    for path, want in ((m_path, m), (m_5_50_path, m_common)):
        got = [float(x) for x in open(path).read().split()]
        if len(got) != n_annot:
            sys.exit('%s has %d values, expected %d'
                     % (path, len(got), n_annot))
        for c in range(n_annot):
            if abs(got[c] - want[c]) > 1e-4 * max(abs(want[c]), 1.0):
                sys.exit('%s column %d: got %g, expected %g'
                         % (path, c, got[c], want[c]))
    print('.M and .M_5_50 match the genotypes (%s)'
          % ', '.join('%s=%g/%g' % (annot_names[c], m[c], m_common[c])
                      for c in range(n_annot)))


if __name__ == '__main__':
    main()
