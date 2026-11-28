#!/usr/bin/env python3
"""Writes a random BGEN v1.2 layout-2 file with multiallelic variants, for
plink2 --bgen import tests.

Usage:
  bgen_multiallelic_gen.py <output .bgen> <seed> <phased | unphased | mixed>
                           <bits> <none | zlib> [sample_ct] [variant_ct]
                           [min_allele_ct]

Each variant gets 2-6 alleles (mostly 3+), and each sample a ploidy of 0, 1
or 2 (mostly 2) and a missingness flag (rare).  Probabilities are drawn from
a mix of exact calls, near-exact calls (to exercise the hardcall and erase
thresholds), and diffuse vectors, then quantized to <bits> bits with
largest-remainder rounding.  'mixed' alternates phased and unphased
variants.  Sample IDs are embedded (s0, s1, ...).
"""

import random
import struct
import sys
import zlib


def quantize(probs, max_val):
    tot = sum(probs)
    scaled = [p * max_val / tot for p in probs]
    numers = [int(x) for x in scaled]
    deficit = max_val - sum(numers)
    order = sorted(range(len(probs)), key=lambda i: -(scaled[i] - numers[i]))
    for i in order[:deficit]:
        numers[i] += 1
    return numers


def rand_vec(rng, ct):
    mode = rng.random()
    if mode < 0.35:
        vec = [0.0] * ct
        vec[rng.randrange(ct)] = 1.0
        return vec
    if mode < 0.65:
        vec = [rng.random() * 0.04 for _ in range(ct)]
        vec[rng.randrange(ct)] += 1.0
        return vec
    if mode < 0.8:
        # two-allele mixture
        vec = [0.0] * ct
        a = rng.randrange(ct)
        b = rng.randrange(ct)
        w = rng.random()
        vec[a] += w
        vec[b] += 1.0 - w
        return vec
    return [rng.random() ** 3 for _ in range(ct)]


def pack_bits(vals, bits):
    acc = 0
    acc_bits = 0
    out = bytearray()
    for v in vals:
        acc |= v << acc_bits
        acc_bits += bits
        while acc_bits >= 8:
            out.append(acc & 255)
            acc >>= 8
            acc_bits -= 8
    if acc_bits:
        out.append(acc & 255)
    return bytes(out)


def main():
    out_fname = sys.argv[1]
    rng = random.Random(int(sys.argv[2]))
    phase_mode = sys.argv[3]
    bits = int(sys.argv[4])
    compressed = (sys.argv[5] == "zlib")
    sample_ct = int(sys.argv[6]) if len(sys.argv) > 6 else 40
    variant_ct = int(sys.argv[7]) if len(sys.argv) > 7 else 30
    min_allele_ct = int(sys.argv[8]) if len(sys.argv) > 8 else 2
    max_val = (1 << bits) - 1
    sample_block = bytearray()
    sample_block += struct.pack("<I", sample_ct)
    for sidx in range(sample_ct):
        sid = ("s%d" % sidx).encode()
        sample_block += struct.pack("<H", len(sid)) + sid
    sample_block = struct.pack("<I", len(sample_block) + 4) + bytes(sample_block)
    flags = (1 if compressed else 0) | (2 << 2) | (1 << 31)
    header = struct.pack("<III", 20, variant_ct, sample_ct) + b"bgen" + struct.pack("<I", flags)
    body = bytearray()
    for vidx in range(variant_ct):
        allele_ct = max(min_allele_ct, rng.choice([2, 3, 3, 3, 4, 4, 5, 6]))
        if phase_mode == "mixed":
            phased = vidx % 2
        else:
            phased = 1 if phase_mode == "phased" else 0
        varid = ("v%d" % vidx).encode()
        chrom = b"1"
        body += struct.pack("<H", len(varid)) + varid
        body += struct.pack("<H", len(varid)) + varid
        body += struct.pack("<H", len(chrom)) + chrom
        body += struct.pack("<IH", vidx + 1, allele_ct)
        for aidx in range(allele_ct):
            allele = b"A" + b"C" * aidx
            body += struct.pack("<I", len(allele)) + allele
        ploidies = []
        vals = []
        for _ in range(sample_ct):
            r = rng.random()
            ploidy = 0 if r < 0.03 else (1 if r < 0.15 else 2)
            missing = rng.random() < 0.05
            ploidies.append(ploidy | (128 if missing else 0))
            if not ploidy:
                continue
            if phased or (ploidy == 1):
                for _ in range(ploidy):
                    numers = quantize(rand_vec(rng, allele_ct), max_val)
                    if missing:
                        numers = [0] * allele_ct
                    vals.extend(numers[:-1])
            else:
                geno_ct = allele_ct * (allele_ct + 1) // 2
                numers = quantize(rand_vec(rng, geno_ct), max_val)
                if missing:
                    numers = [0] * geno_ct
                vals.extend(numers[:-1])
        block = struct.pack("<IHBB", sample_ct, allele_ct, min(p & 63 for p in ploidies), max(p & 63 for p in ploidies))
        block += bytes(ploidies) + bytes([phased, bits]) + pack_bits(vals, bits)
        if compressed:
            comp = zlib.compress(block)
            body += struct.pack("<II", len(comp) + 4, len(block)) + comp
        else:
            body += struct.pack("<I", len(block)) + block
    # first 4 bytes: offset of the first variant block, relative to byte 4
    data = struct.pack("<I", len(header) + len(sample_block)) + header + sample_block + bytes(body)
    open(out_fname, "wb").write(data)


main()
