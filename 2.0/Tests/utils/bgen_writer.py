#!/usr/bin/env python3
"""Minimal BGEN writer for tests (no external dependencies).

plink2's own --export bgen-1.x only emits a few fixed layouts (one bit depth
per run, constant ploidy per variant, embedded sample IDs, zstd/zlib), so
BGEN fixtures built with it can't reach much of the importer.  This writes:
  * layout 1 (v1.1): three uint16 probabilities per sample, scale 32768,
    optional zlib compression;
  * layout 2 (v1.2/v1.3): any bit depth 1..32, phased or unphased,
    per-sample ploidy (variable-ploidy records get min != max) and missing
    flags, no compression / zlib / zstd;
  * with or without the sample identifier block.

Probabilities are passed as integer numerators, already quantized by the
caller, so the caller knows the exact encoded values.

Variant dict fields:
  chrom, pos, varid, rsid, alleles (list of str)
  layout 2: phased (bool), ploidy (list[int]), missing (list[bool]),
            numers (list of per-sample lists of stored numerators, in BGEN
            order; the implied last probability is omitted)
  layout 1: probs (list of per-sample (n0, n1, n2) uint16 triples)
"""

import shutil
import struct
import subprocess
import zlib

try:
    from compression import zstd as _zstd  # Python 3.14+
except ImportError:
    _zstd = None


def zstd_available():
    return (_zstd is not None) or (shutil.which("zstd") is not None)


def _zstd_compress(data):
    if _zstd is not None:
        return _zstd.compress(data)
    return subprocess.run(["zstd", "-q", "-c", "-3"], input=data,
                          stdout=subprocess.PIPE, check=True).stdout


def _str16(s):
    b = s.encode()
    return struct.pack("<H", len(b)) + b


def _str32(s):
    b = s.encode()
    return struct.pack("<I", len(b)) + b


def _pack_bits(vals, bits):
    acc = 0
    nbits = 0
    for i, v in enumerate(vals):
        assert 0 <= v < (1 << bits)
        acc |= v << (i * bits)
        nbits += bits
    return acc.to_bytes((nbits + 7) // 8, "little")


def _layout2_geno(v, sample_ct, bits):
    ploidy = v["ploidy"]
    missing = v["missing"]
    out = struct.pack("<IHBB", sample_ct, len(v["alleles"]), min(ploidy),
                      max(ploidy))
    out += bytes([(p & 63) | (0x80 if m else 0) for p, m in zip(ploidy, missing)])
    out += bytes([1 if v["phased"] else 0, bits])
    vals = []
    for s in v["numers"]:
        vals.extend(s)
    return out + _pack_bits(vals, bits)


def _layout1_geno(v):
    out = b""
    for n0, n1, n2 in v["probs"]:
        out += struct.pack("<HHH", n0, n1, n2)
    return out


def write_bgen(path, sample_ids, variants, layout=2, bits=16,
               compression="zlib", embed_ids=True):
    sample_ct = len(sample_ids)
    comp_code = {"none": 0, "zlib": 1, "zstd": 2}[compression]
    assert layout in (1, 2)
    assert not (layout == 1 and comp_code == 2)
    flags = comp_code | (layout << 2) | ((1 << 31) if embed_ids else 0)
    header = struct.pack("<III", 20, len(variants), sample_ct) + b"bgen"
    header += struct.pack("<I", flags)
    sid_block = b""
    if embed_ids:
        body = struct.pack("<I", sample_ct) + b"".join(_str16(s) for s in sample_ids)
        sid_block = struct.pack("<I", len(body) + 4) + body
    out = struct.pack("<I", len(header) + len(sid_block)) + header + sid_block
    for v in variants:
        rec = b""
        if layout == 1:
            rec += struct.pack("<I", sample_ct)
        rec += _str16(v["varid"]) + _str16(v["rsid"]) + _str16(v["chrom"])
        rec += struct.pack("<I", v["pos"])
        if layout == 1:
            assert len(v["alleles"]) == 2
        else:
            rec += struct.pack("<H", len(v["alleles"]))
        rec += b"".join(_str32(a) for a in v["alleles"])
        if layout == 1:
            geno = _layout1_geno(v)
            if comp_code:
                c = zlib.compress(geno)
                rec += struct.pack("<I", len(c)) + c
            else:
                rec += geno
        else:
            geno = _layout2_geno(v, sample_ct, bits)
            if comp_code:
                c = zlib.compress(geno) if comp_code == 1 else _zstd_compress(geno)
                rec += struct.pack("<II", len(c) + 4, len(geno)) + c
            else:
                rec += struct.pack("<I", len(geno)) + geno
        out += rec
    with open(path, "wb") as f:
        f.write(out)


def write_sample_file(path, sample_ids):
    with open(path, "w") as f:
        f.write("ID_1 ID_2 missing\n0 0 0\n")
        for s in sample_ids:
            f.write("%s %s 0\n" % (s, s))
