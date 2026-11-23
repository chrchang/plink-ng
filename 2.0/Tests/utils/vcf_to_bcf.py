#!/usr/bin/env python3
"""Minimal VCF -> BCF 2.2 converter for tests (no htslib needed).

usage: vcf_to_bcf.py [--legacy-phase] <in.vcf> <out.bcf>

plink2's own --export bcf only emits GT and dosage fields, so BCF fixtures
built with it can never exercise the importer's handling of e.g. GQ/DP.
This writes any combination of FORMAT fields declared in the header as
Type=String (GT only), Type=Integer or Type=Float, with htslib's encoding:
  * GT: (allele + 1) << 1 | phased, missing allele 0, int8 when every value
    fits and int16 otherwise (so a site with >63 ALT alleles exercises the
    int16 path); samples with lower ploidy than the widest one in the record
    are padded with END_OF_VECTOR, which is how mixed ploidy is represented.
    VCFv4.4+ prefix phasing ("|0/1") is supported; see encode_gt().
    --legacy-phase clears every first-allele phase bit, as HTSlib before 1.22
    did (so "0|1" is 0x02 0x05); input must not use prefix phasing then.
  * Integer: smallest of int8/int16/int32 that holds the record's values,
    '.' -> MISSING, short vectors padded with END_OF_VECTOR.
  * Float: '.' -> 0x7F800001, padding 0x7F800002.
INFO must be '.'.  FILTER must be '.' or 'PASS'.  ID and QUAL may be '.'.
"""

import re
import struct
import sys
import zlib

INT_TYPES = (
    # (BCF type code, struct format, MISSING, END_OF_VECTOR, min valid, max)
    (1, "<b", -128, -127, -120, 127),
    (2, "<h", -32768, -32767, -32760, 32767),
    (3, "<i", -2147483648, -2147483647, -2147483640, 2147483647),
)
FLOAT_MISSING = struct.pack("<I", 0x7F800001)
FLOAT_EOV = struct.pack("<I", 0x7F800002)


def typed_int(x):
    for typ, fmt, _, _, lo, hi in INT_TYPES:
        if lo <= x <= hi:
            return bytes([0x10 | typ]) + struct.pack(fmt, x)
    raise ValueError(x)


def typed_header(ct, typ):
    if ct < 15:
        return bytes([(ct << 4) | typ])
    return bytes([0xF0 | typ]) + typed_int(ct)


def typed_str(s):
    b = s.encode()
    return typed_header(len(b), 7) + b


def bgzf(raw):
    out = b""
    for i in range(0, len(raw), 60000):
        chunk = raw[i:i + 60000]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        cdata = c.compress(chunk) + c.flush()
        out += (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
                + struct.pack("<H", len(cdata) + 25) + cdata
                + struct.pack("<II", zlib.crc32(chunk), len(chunk)))
    out += bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    return out


LEGACY_PHASE = False


def encode_gt(gt):
    """Returns (values, slot count), matching htslib 1.22.

    The first allele's phase bit is the VCFv4.4+ leading phasing indicator if
    present (e.g. "|0/1"); otherwise it's implicitly '|' iff no other
    indicator is '/' (for a haploid call: iff the allele isn't '.').  htslib
    reserves an extra vector slot for an explicit leading indicator, so e.g. a
    record where every call is diploid but one is "|0/1" gets width-3 GT
    vectors, padded with END_OF_VECTOR.
    """
    prefix = gt[0] if gt[0] in "/|" else ""
    parts = re.split(r"([/|])", gt[len(prefix):])
    seps = parts[1::2]
    if LEGACY_PHASE:
        assert not prefix, "--legacy-phase input can't use prefix phasing"
        first_phased = False
    elif prefix:
        first_phased = prefix == "|"
    elif seps:
        first_phased = "/" not in seps
    else:
        first_phased = parts[0] != "."
    vals = []
    for i in range(0, len(parts), 2):
        a = parts[i]
        v = 0 if a == "." else (int(a) + 1) << 1
        if (parts[i - 1] == "|") if i else first_phased:
            v |= 1
        vals.append(v)
    return vals, len(vals) + (prefix != "")


def gt_vec(strs):
    encoded = [encode_gt(s) for s in strs]
    vecs = [v for v, _ in encoded]
    width = max(slot_ct for _, slot_ct in encoded)
    mx = max(max(v) for v in vecs)
    typ, fmt, _, eov = INT_TYPES[0][:4] if mx <= 127 else INT_TYPES[1][:4]
    out = typed_header(width, typ)
    for v in vecs:
        for k in range(width):
            out += struct.pack(fmt, v[k] if k < len(v) else eov)
    return out


def int_vec(strs):
    vecs = [[None if x == "." else int(x) for x in s.split(",")] for s in strs]
    width = max(len(v) for v in vecs)
    present = [x for v in vecs for x in v if x is not None]
    lo = min(present, default=0)
    hi = max(present, default=0)
    for typ, fmt, miss, eov, vlo, vhi in INT_TYPES:
        if vlo <= lo and hi <= vhi:
            break
    out = typed_header(width, typ)
    for v in vecs:
        for k in range(width):
            if k >= len(v):
                x = eov
            elif v[k] is None:
                x = miss
            else:
                x = v[k]
            out += struct.pack(fmt, x)
    return out


def float_vec(strs):
    vecs = [s.split(",") for s in strs]
    width = max(len(v) for v in vecs)
    out = typed_header(width, 5)
    for v in vecs:
        for k in range(width):
            if k >= len(v):
                out += FLOAT_EOV
            elif v[k] == ".":
                out += FLOAT_MISSING
            else:
                out += struct.pack("<f", float(v[k]))
    return out


def main():
    global LEGACY_PHASE
    args = sys.argv[1:]
    if args and args[0] == "--legacy-phase":
        LEGACY_PHASE = True
        args = args[1:]
    vcf_fname, bcf_fname = args[:2]
    header_lines = ['##FILTER=<ID=PASS,Description="All filters passed",IDX=0>']
    dict_idx = {"PASS": 0}
    fmt_types = {}
    contigs = []
    records = []
    sample_line = None
    with open(vcf_fname) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                m = re.match(r"##(contig|FORMAT|INFO|FILTER)=<ID=([^,>]+)", line)
                if m:
                    kind, fid = m.groups()
                    if kind == "contig":
                        line = line[:-1] + ",IDX=%d>" % len(contigs)
                        contigs.append(fid)
                    elif fid != "PASS":
                        if fid not in dict_idx:
                            dict_idx[fid] = len(dict_idx)
                        line = line[:-1] + ",IDX=%d>" % dict_idx[fid]
                        if kind == "FORMAT":
                            fmt_types[fid] = re.search(r"Type=(\w+)", line).group(1)
                    else:
                        continue
                header_lines.append(line)
            elif line.startswith("#CHROM"):
                sample_line = line
            elif line:
                records.append(line.split("\t"))
    # ##fileformat must come first.
    ff = [l for l in header_lines if l.startswith("##fileformat")]
    rest = [l for l in header_lines if not l.startswith("##fileformat")]
    text = "\n".join(ff + rest + [sample_line]) + "\n"
    raw = b"BCF\x02\x02" + struct.pack("<I", len(text) + 1) + text.encode() + b"\0"
    for f in records:
        assert f[7] == ".", "INFO must be '.'"
        alleles = [f[3]] + ([] if f[4] == "." else f[4].split(","))
        qual = FLOAT_MISSING if f[5] == "." else struct.pack("<f", float(f[5]))
        fmt_keys = f[8].split(":") if len(f) > 8 else []
        samples = [s.split(":") for s in f[9:]]
        shared = struct.pack("<iii", contigs.index(f[0]), int(f[1]) - 1, len(f[3]))
        shared += qual
        shared += struct.pack("<II", len(alleles) << 16,
                              len(samples) | (len(fmt_keys) << 24))
        shared += typed_str("" if f[2] == "." else f[2])
        for a in alleles:
            shared += typed_str(a)
        if f[6] == ".":
            shared += bytes([0])
        else:
            assert f[6] == "PASS"
            shared += typed_header(1, 1) + bytes([0])
        indiv = b""
        for k, key in enumerate(fmt_keys):
            indiv += typed_int(dict_idx[key])
            vals = [s[k] if k < len(s) else "." for s in samples]
            if key == "GT":
                indiv += gt_vec(vals)
            elif fmt_types[key] == "Integer":
                indiv += int_vec(vals)
            elif fmt_types[key] == "Float":
                indiv += float_vec(vals)
            else:
                raise ValueError("unsupported FORMAT type for " + key)
        raw += struct.pack("<II", len(shared), len(indiv)) + shared + indiv
    with open(bcf_fname, "wb") as f:
        f.write(bgzf(raw))


if __name__ == "__main__":
    main()
