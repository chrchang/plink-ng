#!/usr/bin/env python3
"""Recover the INFO columns from the source VCF and compare with .vinfo.

This parses the VCF independently of plink2, so the check is against the input
rather than against another plink2 run.
"""
import sys


def parse_vcf(path):
    flags = set()
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('##INFO=<'):
                body = line[len('##INFO=<'):]
                key = body.split('ID=', 1)[1].split(',', 1)[0].split('>', 1)[0]
                if ',Number=0,' in body:
                    flags.add(key)
                continue
            if line.startswith('#'):
                continue
            f8 = line.rstrip('\n').split('\t')
            info = {}
            if f8[7] != '.':
                for sub in f8[7].split(';'):
                    if not sub:
                        continue
                    if '=' in sub:
                        k, v = sub.split('=', 1)
                        info[k] = v
                    else:
                        info[sub] = None
            rows.append((f8[0], f8[1], f8[2], f8[3], f8[4], info))
    return flags, rows


def main():
    vcf, vinfo = sys.argv[1], sys.argv[2]
    flags, rows = parse_vcf(vcf)
    with open(vinfo) as f:
        lines = [ln.rstrip('\n').split('\t') for ln in f]
    hdr = lines[0]
    assert hdr[0] == '#CHROM', hdr
    keys = hdr[5:]
    body = lines[1:]
    assert len(body) == len(rows), 'row count %d vs %d' % (len(body), len(rows))
    for got, (chrom, pos, vid, ref, alt, info) in zip(body, rows):
        assert got[:5] == [chrom, pos, vid, ref, alt], (got[:5], chrom, pos, vid)
        for i, key in enumerate(keys):
            cell = got[5 + i]
            if key in flags:
                want = '1' if key in info else '0'
            elif key not in info:
                want = 'NA'
            elif info[key] is None:
                want = '1'
            else:
                want = info[key]
            assert cell == want, 'variant %s key %s: %r vs %r' % (vid, key, cell, want)
    print('%d rows x %d keys verified against %s' % (len(body), len(keys), vcf))


if __name__ == '__main__':
    main()
