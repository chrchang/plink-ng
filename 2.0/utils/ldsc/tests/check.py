#!/usr/bin/env python3
"""Compare ldsc's output file against the oracle's, field by field."""
import subprocess
import sys


def read_result(path):
    """One-row file -> {column: value}; multi-row -> {column_rowidx: value}."""
    with open(path) as f:
        header = f.readline().split('\t')
        header = [h.strip() for h in header]
        rows = [line.rstrip('\n').split('\t') for line in f if line.strip()]
    if len(rows) == 1:
        return dict(zip(header, rows[0]))
    out = {}
    for row_idx, row in enumerate(rows):
        for col_idx, name in enumerate(header):
            out['%s_%d' % (name, row_idx)] = row[col_idx]
    return out


def main():
    result_path = sys.argv[1]
    oracle_cmd = sys.argv[2:]
    got = {}
    for path in result_path.split(','):
        got.update(read_result(path))
    out = subprocess.run([sys.executable] + oracle_cmd, check=True,
                         capture_output=True, text=True).stdout
    want = {}
    for line in out.split('\n'):
        if not line.strip():
            continue
        key, val = line.split()
        want[key] = float(val)
    checked = 0
    for key, expected in want.items():
        if key not in got:
            sys.exit('%s: oracle reported %s, which %s does not contain' %
                     (key, key, result_path))
        if got[key] == 'NA':
            sys.exit('%s: expected %.12g, got NA' % (key, expected))
        actual = float(got[key])
        # Relative, with a small absolute floor: a few of these quantities are
        # exactly zero up to cancellation (the standard error of a proportion
        # that has to be 1, for instance), and the two implementations round
        # that cancellation differently.
        scale = max(abs(expected), abs(actual))
        if abs(actual - expected) > 1e-6 * scale + 1e-9:
            sys.exit('%s: expected %.12g, got %.12g' % (key, expected, actual))
        checked += 1
    if not checked:
        sys.exit('the oracle reported nothing; the comparison proves nothing')
    print('%d values match the oracle (%s)' % (checked, result_path))


if __name__ == '__main__':
    main()
