#!/usr/bin/env python3
# Writes the small .bcf files (and matching .vcf references) used by
# run_tests.sh.
import gzip
import struct

GT = [(0, 1), (1, 1), (0, 0), (0, 1)]  # 4 samples
ALLELES = [b"A", b"C"]


def tstr(b):
    if len(b) < 15:
        return bytes([len(b) << 4 | 7]) + b
    return b"\xf7\x13" + struct.pack("<i", len(b)) + b


def record(pos, vid, gt_key, gts):
    n_sample = len(gts)
    n_fmt = 1 if n_sample else 0
    shared = struct.pack("<iiif", 0, pos - 1, 1, 0.0)
    shared += struct.pack("<II", len(ALLELES) << 16, (n_fmt << 24) | n_sample)
    shared += tstr(vid) + b"".join(tstr(a) for a in ALLELES) + b"\x00"
    fmt = b""
    if n_sample:
        # GT: int8, 2 values per sample, (allele + 1) << 1 (unphased)
        fmt = b"\x11" + bytes([gt_key]) + bytes([2 << 4 | 1])
        fmt += b"".join(bytes([(a + 1) << 1, (b + 1) << 1]) for a, b in gts)
    return struct.pack("<II", len(shared), len(fmt)) + shared + fmt


def write(name, meta, gt_key, gts):
    n_sample = len(gts)
    samples = [b"S%d" % (i + 1) for i in range(n_sample)]
    head = b"##fileformat=VCFv4.2\n" + b"".join(m + b"\n" for m in meta)
    cols = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if n_sample:
        cols += b"\tFORMAT\t" + b"\t".join(samples)
    text = head + cols + b"\n\0"
    body = b"".join(record(100 * (k + 1), b"rs%d" % (k + 1), gt_key, gts) for k in range(3))
    with gzip.open(name + ".bcf", "wb") as f:
        f.write(b"BCF\x02\x02" + struct.pack("<I", len(text)) + text + body)
    with open(name + ".vcf", "wb") as f:
        f.write(head + cols + b"\n")
        for k in range(3):
            line = b"1\t%d\trs%d\tA\tC\t0\t.\t." % (100 * (k + 1), k + 1)
            if n_sample:
                line += b"\tGT\t" + b"\t".join(b"%d/%d" % g for g in gts)
            f.write(line + b"\n")


PASS = b'##FILTER=<ID=PASS,Description="ok">'
GTH = b'##FORMAT=<ID=GT,Number=1,Type=String,Description="g">'
CONTIG = b"##contig=<ID=1,length=1000>"
# DP is both an INFO and a FORMAT field, so it has one dictionary entry:
# PASS = 0, DP = 1, GT = 2.
write("dup_id", [PASS, b'##INFO=<ID=DP,Number=1,Type=Integer,Description="d">',
                 b'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">', GTH, CONTIG], 2, GT)
# An ID-only contig line is legal.
write("contig_id_only", [PASS, GTH, b"##contig=<ID=1>"], 1, GT)
# A FILTER line with no comma.
write("filter_no_comma", [b"##FILTER=<ID=PASS>", GTH, CONTIG], 1, GT)
# No samples at all.
write("no_samples", [PASS, CONTIG], 1, [])
