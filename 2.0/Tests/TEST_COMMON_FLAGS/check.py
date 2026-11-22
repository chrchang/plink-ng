#!/usr/bin/env python3
"""Recomputes values from tmp_data.ped, independently of plink."""
import struct
import sys


def load_ped():
    """Returns ALT (C) allele counts per sample, then per variant in map order."""
    with open("tmp_data.ped") as f:
        rows = [line.split() for line in f]
    return [[(r[i] == "C") + (r[i + 1] == "C") for i in range(6, len(r), 2)] for r in rows]


def close(a, b, tol):
    if abs(a - b) > tol:
        sys.exit(f"{a} != {b}")


def check_hudson_fst():
    """Hudson Fst per variant and ratio of averages; odd samples are P1, even ones P2."""
    counts = load_ped()
    pops = [counts[0::2], counts[1::2]]
    with open("plink2_fst.P1.P2.fst.var") as f:
        reported = [line.split() for line in f][1:]
    assert len(reported) == len(counts[0])
    num_sum = 0.0
    denom_sum = 0.0
    for v, row in enumerate(reported):
        n1, n2 = (2 * len(pop) for pop in pops)
        p1, p2 = (sum(s[v] for s in pop) / (2 * len(pop)) for pop in pops)
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        denom = p1 * (1 - p2) + p2 * (1 - p1)
        close(num / denom, float(row[4]), 1e-5)
        num_sum += num
        denom_sum += denom
    with open("plink2_fst.fst.summary") as f:
        close(num_sum / denom_sum, float(f.readlines()[1].split()[2]), 1e-5)


def check_score_list():
    """Two score files: weights {v0: 1, v2: 2} and {v2: 3, v3 on the REF allele A: 1}."""
    counts = load_ped()
    with open("plink2_score.sscore") as f:
        lines = f.readlines()[1:]
    assert len(lines) == len(counts)
    for line, c in zip(lines, counts):
        fields = line.split()
        close((c[0] + 2 * c[2]) / 4, float(fields[2]), 1e-5)
        close((3 * c[2] + (2 - c[3])) / 4, float(fields[3]), 1e-5)


def check_grm():
    """--make-grm-bin matches plink 1.9."""
    for suffix, tol in (("grm.bin", 1e-5), ("grm.N.bin", 0)):
        with open("plink_grm." + suffix, "rb") as f1, open("plink2_grm." + suffix, "rb") as f2:
            b1, b2 = f1.read(), f2.read()
        assert len(b1) == len(b2) == 4 * 120 * 121 // 2
        for x1, x2 in zip(struct.iter_unpack("f", b1), struct.iter_unpack("f", b2)):
            close(x1[0], x2[0], tol)


{"fst": check_hudson_fst, "score": check_score_list, "grm": check_grm}[sys.argv[1]]()
