#!/usr/bin/env python3
# Writes the small .bcf files used by run_tests.sh.
import gzip
import struct


def tstr(b):
    if len(b) < 15:
        return bytes([len(b) << 4 | 7]) + b
    return b"\xf7\x13" + struct.pack("<i", len(b)) + b


def write(fname, contigs, gt_type, gts, n_allele_names=(b"A", b"C")):
    text = b'##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="ok">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
    for contig in contigs:
        text += contig + b"\n"
    n_sample = len(gts) // 2
    text += b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    text += b"".join(b"\tS%d" % (i + 1) for i in range(n_sample)) + b"\n\0"
    shared = struct.pack("<iiif", 0, 99, 1, 0.0) + struct.pack("<II", len(n_allele_names) << 16, (1 << 24) | n_sample)
    shared += tstr(b"rs1")
    for allele in n_allele_names:
        shared += tstr(allele)
    shared += b"\x00"
    fmt = b"\x11\x01" + bytes([2 << 4 | gt_type])
    fmt += b"".join(struct.pack({1: "<b", 2: "<h", 3: "<i"}[gt_type], g) for g in gts)
    body = struct.pack("<II", len(shared), len(fmt)) + shared + fmt
    with gzip.open(fname, "wb") as f:
        f.write(b"BCF\x02\x02" + struct.pack("<I", len(text)) + text + body)


ok_gts = [2, 4, 4, 4]
one = [b"##contig=<ID=1,length=1000>"]
write("ok.bcf", one, 1, ok_gts)
write("gt8.bcf", one, 1, [2, 12, 12, 12])
write("gt32.bcf", one, 3, [0x7ffffffe, 0x7ffffffe, 4, 4])
write("idx_ok.bcf", [b"##contig=<ID=1,IDX=0>", b"##contig=<ID=2,IDX=1>"], 1, ok_gts)
write("idx_swapped.bcf", [b"##contig=<ID=1,IDX=1>", b"##contig=<ID=2,IDX=0>"], 1, ok_gts)
# 10 samples: values reach the word-at-a-time path, not just the tail loop.
many = [b"A"] + [b"C%d" % i for i in range(69)]
write("wide_ok.bcf", one, 1, [2, 4] * 10)
write("wide_gt8.bcf", one, 1, [2, 4] * 3 + [2, 8] + [2, 4] * 6)
write("wide_gt16.bcf", one, 2, [2, 4] * 3 + [2, 8] + [2, 4] * 6)
# 70 alleles with int8 GT (htslib picks the narrowest type that fits the data).
write("many_alleles_ok.bcf", one, 1, [2, 0x7e] * 10, many)
