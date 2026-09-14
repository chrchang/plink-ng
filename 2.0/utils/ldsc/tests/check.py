#!/usr/bin/env python3
"""Compare ldsc's output file against the oracle's, field by field."""
import subprocess
import sys


def read_result(path):
    with open(path) as f:
        header = f.readline().split()
        values = f.readline().split()
    return dict(zip(header, values))


def main():
    result_path = sys.argv[1]
    oracle_cmd = sys.argv[2:]
    got = read_result(result_path)
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
        actual = float(got[key])
        scale = max(abs(expected), 1e-12)
        if abs(actual - expected) / scale > 1e-6:
            sys.exit('%s: expected %.12g, got %.12g' % (key, expected, actual))
        checked += 1
    if not checked:
        sys.exit('the oracle reported nothing; the comparison proves nothing')
    print('%d values match the oracle (%s)' % (checked, result_path))


if __name__ == '__main__':
    main()
