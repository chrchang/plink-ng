#!/usr/bin/env python3
"""Consistency checks between --make-king outputs, --make-king-table, and the
--make-grm-list/--make-grm-sparse pair.

Usage:
  check_king.py king <.kin0> <.king.id> <shape> <file>
      shape is one of triangle, square, square0, bin8tri, bin4tri, bin8sq.
  check_king.py sparse <.grm list file> <.grm.sp file> <threshold>
"""

import math
import struct
import sys


def fail(msg):
    sys.exit("check_king.py: " + msg)


def close(a, b, tol):
    return abs(a - b) <= tol


def load_kin0(fname):
    kin = {}
    with open(fname) as f:
        header = f.readline().rstrip("\n").lstrip("#").split("\t")
        i1 = header.index("IID1")
        i2 = header.index("IID2")
        ik = header.index("KINSHIP")
        for line in f:
            row = line.rstrip("\n").split("\t")
            kin[frozenset((row[i1], row[i2]))] = float(row[ik])
    return kin


def load_ids(fname):
    with open(fname) as f:
        lines = [line.rstrip("\n") for line in f]
    if lines[0].startswith("#"):
        lines = lines[1:]
    return [line.split("\t")[-1] for line in lines]


def check_king(kin0, idfile, shape, fname):
    kin = load_kin0(kin0)
    ids = load_ids(idfile)
    n = len(ids)
    if len(kin) != n * (n - 1) // 2:
        fail("%s has %d pairs, expected %d" % (kin0, len(kin), n * (n - 1) // 2))
    if shape.startswith("bin"):
        width = 8 if shape.startswith("bin8") else 4
        fmt = "<d" if width == 8 else "<f"
        with open(fname, "rb") as f:
            data = f.read()
        vals = [struct.unpack_from(fmt, data, i * width)[0] for i in range(len(data) // width)]
        if shape.endswith("tri"):
            if len(vals) != n * (n - 1) // 2:
                fail("%s has %d values, expected %d" % (fname, len(vals), n * (n - 1) // 2))
            rows = []
            pos = 0
            for i in range(1, n):
                rows.append(vals[pos:pos + i])
                pos += i
            matrix = [[]] + rows
        else:
            if len(vals) != n * n:
                fail("%s has %d values, expected %d" % (fname, len(vals), n * n))
            matrix = [vals[i * n:(i + 1) * n] for i in range(n)]
        # The text formats print ~6 significant digits.
        tol = 1e-5 if width == 8 else 1e-4
    else:
        with open(fname) as f:
            matrix = [[float(x) for x in line.split("\t")] for line in f]
        if shape == "triangle":
            matrix = [[]] + matrix
        tol = 1e-6
    for i in range(n):
        expected_len = i if shape in ("triangle", "bin8tri", "bin4tri") else n
        if len(matrix[i]) != expected_len:
            fail("%s row %d has %d entries, expected %d" % (fname, i, len(matrix[i]), expected_len))
        for j in range(len(matrix[i])):
            v = matrix[i][j]
            if i == j:
                expected = 0.5
            elif j > i:
                if shape == "square0":
                    expected = 0.0
                else:
                    expected = kin[frozenset((ids[i], ids[j]))]
            else:
                expected = kin[frozenset((ids[i], ids[j]))]
            if not close(v, expected, tol):
                fail("%s [%d][%d] = %g, expected %g" % (fname, i, j, v, expected))
    if shape == "square" or shape == "bin8sq":
        for i in range(n):
            for j in range(n):
                if not close(matrix[i][j], matrix[j][i], 1e-12):
                    fail("%s is not symmetric at [%d][%d]" % (fname, i, j))


def check_sparse(grm_list, grm_sp, threshold):
    threshold = float(threshold)
    expected = {}
    with open(grm_list) as f:
        for line in f:
            row = line.split("\t")
            val = float(row[3])
            if val >= threshold:
                expected[(int(row[0]) - 1, int(row[1]) - 1)] = val
    seen = set()
    with open(grm_sp) as f:
        for line in f:
            row = line.split("\t")
            key = (int(row[0]), int(row[1]))
            if key not in expected:
                fail("%s has unexpected entry %s" % (grm_sp, key))
            if key in seen:
                fail("%s has duplicate entry %s" % (grm_sp, key))
            seen.add(key)
            # .grm (list) prints 6 significant digits.
            if not math.isclose(float(row[2]), expected[key], rel_tol=1e-5, abs_tol=1e-6):
                fail("%s %s = %s, expected %g" % (grm_sp, key, row[2].strip(), expected[key]))
    if seen != set(expected):
        fail("%s is missing entries %s" % (grm_sp, sorted(set(expected) - seen)))


def main():
    if len(sys.argv) == 6 and sys.argv[1] == "king":
        check_king(*sys.argv[2:])
    elif len(sys.argv) == 5 and sys.argv[1] == "sparse":
        check_sparse(*sys.argv[2:])
    else:
        sys.exit(__doc__)


main()
