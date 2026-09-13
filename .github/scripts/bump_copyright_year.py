#!/usr/bin/env python3
"""Bump the end year in PLINK's own copyright lines.

Two line forms are rewritten, and only those:

  // This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
  "(C) 2005-2026 S Purcell, C Chang, B Demaille      GNU General Public License v3"

The first identifies the file as PLINK's own, so the bundled third-party
sources (htslib, zstd, libdeflate, SFMT, QD) cannot match: they name their own
authors.  The second is the runtime banner, distinguished from the first by the
licence text on the same line; only the digits are rewritten, so the column
padding is left alone.

Usage: bump_copyright_year.py <year> [root]
Prints the files it changed, one per line.
"""

import os
import re
import sys

SUFFIXES = ('.c', '.cc', '.h', '.hpp', '.cu')

HEADER_RE = re.compile(
    r'(part of (?:the )?PLINK[^,\n]*, copyright \(C\) [0-9]{4})-[0-9]{4}')
BANNER_RE = re.compile(
    r'(\(C\) [0-9]{4})-[0-9]{4}(?=[^\n]*GNU General Public License v3)')


def main():
    if len(sys.argv) < 2:
        sys.exit('usage: bump_copyright_year.py <year> [root]')
    year = sys.argv[1]
    if not re.fullmatch(r'[0-9]{4}', year):
        sys.exit('year must be four digits')
    root = sys.argv[2] if len(sys.argv) > 2 else '.'

    changed = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != '.git']
        for name in filenames:
            if not name.endswith(SUFFIXES):
                continue
            path = os.path.join(dirpath, name)
            try:
                with open(path, encoding='utf-8') as f:
                    old = f.read()
            except (UnicodeDecodeError, OSError):
                continue
            new = BANNER_RE.sub(r'\1-' + year, HEADER_RE.sub(r'\1-' + year, old))
            if new != old:
                with open(path, 'w', encoding='utf-8') as f:
                    f.write(new)
                changed.append(os.path.relpath(path, root))
    for path in sorted(changed):
        print(path)


if __name__ == '__main__':
    main()
