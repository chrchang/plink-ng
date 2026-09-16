#!/usr/bin/env python3
"""Compare a --tucc fileset's genotypes against the expected ALT counts."""
import sys


def main():
    raw_path, expected_path = sys.argv[1], sys.argv[2]
    with open(raw_path) as f:
        header = f.readline().split()
        rows = [line.split() for line in f if line.strip()]
    # --export A counts the ALT allele when --export-allele names it, and the
    # column header is '<ID>_<allele>'.
    var_ids = [h.rsplit('_', 1)[0] for h in header[6:]]
    ids = [r[1] for r in rows]
    counts = {}
    for row_idx, row in enumerate(rows):
        for col_idx, vid in enumerate(var_ids):
            counts[(vid, row_idx)] = row[6 + col_idx]
    checked = 0
    nonmissing = 0
    with open(expected_path) as f:
        for line in f:
            fields = line.split()
            vid = fields[0]
            expected = fields[1:]
            if len(expected) != len(rows):
                sys.exit('%s: %d expected values for %d samples' %
                         (vid, len(expected), len(rows)))
            for row_idx, want in enumerate(expected):
                got = counts[(vid, row_idx)]
                if got != want:
                    sys.exit('%s / %s: expected %s, got %s' %
                             (vid, ids[row_idx], want, got))
                checked += 1
                if want != 'NA':
                    nonmissing += 1
    if not nonmissing:
        sys.exit('every expected genotype is missing; the test proves nothing')
    print('%d genotypes verified (%d nonmissing)' % (checked, nonmissing))


if __name__ == '__main__':
    main()
