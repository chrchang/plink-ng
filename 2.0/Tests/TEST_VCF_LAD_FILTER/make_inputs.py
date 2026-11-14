#!/usr/bin/env python3
# Writes lad.vcf, lad.bcf and dp.vcf used by run_tests.sh.
import struct
import zlib

SAMPLES = "abcdefgh"
ALTS12 = ["C", "G", "T", "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC"]
# (id, alts, FORMAT, per-sample fields)
RECORDS = [
    ("v1", ["C"], "GT:GQ:LAD",
     ["0/1:30:5,6", "0/0:30:9", "0/0:15:40", "0/0:30:.", "0/1:30:4,.", "1/1:30:0,12", "./.:.:.", "0/0:30"]),
    ("v2", ["C"], "GT:LAD:GQ",
     ["0/1:5,6:30", "0/0:12:15", "1:12:12", "0/0:9:50", "0/0:3,3,4:25", "1:6:30", "0:4:30", "0/0:10"]),
    ("v3", ["C"], "GT:LAD",
     ["0/1:5,6", "0/0:9", "1:5", "0:4", "0/0:.", "0/0", "1/1:0,10", "./."]),
    ("v4", ALTS12, "GT:GQ:LAD",
     ["11:30:0,6", "11/3:30:0,6", "0/0:30:20", "10:30:4", "2/11:30:5,5", "0/0:30:12", "11|11:30:3,3", "0/0:19:30"]),
    ("v5", ["C"], "GT:GQ:LAD",
     ["0/1:30:100,100", "0/0:30:200", "1:30:6", "0/0:30:9", "1/1:30:0,300", "0/1:10:200,1", "0/0:30:150", "0:30:4"]),
]
# Phased record, VCF only.
VCF_ONLY = [
    ("v6", ["C"], "GT:GQ:LAD",
     ["0|1:30:5,6", "1|0:30:9", "0|1:15:20", "1:30:5", "0|0:30:.", "1|1:30:12", "0|1:30", ".|.:.:."]),
]
HEADER = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=PASS,Description="ok">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">',
    '##FORMAT=<ID=LAD,Number=.,Type=Integer,Description="l">',
    "##contig=<ID=1,length=100000>",
]
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES)

with open("lad.vcf", "w") as f:
    f.write("\n".join(HEADER + [COLS]) + "\n")
    for i, (vid, alts, fmt, fields) in enumerate(RECORDS + VCF_ONLY):
        f.write("\t".join(["1", str(100 * (i + 1)), vid, "A", ",".join(alts), ".", "PASS", ".", fmt] + fields) + "\n")

with open("dp.vcf", "w") as f:
    hdr = HEADER[:4] + ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">'] + HEADER[4:]
    f.write("\n".join(hdr + [COLS]) + "\n")
    f.write("\t".join(["1", "100", "d1", "A", "C", ".", "PASS", ".", "GT:DP:LAD",
                       "0/0:9:50", "0/0:10:1", "0/1:30:.", "1/1:.:2", "0/0:4:9", "0/1:12:3", "0/0:20:20", "1:5:5"]) + "\n")


def tstr(b):
    return bytes([len(b) << 4 | 7]) + b


def typed_int(v):
    return b"\x11" + bytes([v])


MISS = {1: -128, 2: -32768}
EOV = {1: -127, 2: -32767}
PACK = {1: "<b", 2: "<h"}

text = ("\n".join(HEADER + [COLS]) + "\n").encode() + b"\0"
body = b""
for i, (vid, alts, fmt, fields) in enumerate(RECORDS):
    keys = fmt.split(":")
    per = [dict(zip(keys, s.split(":"))) for s in fields]
    alleles = [b"A"] + [a.encode() for a in alts]
    shared = struct.pack("<iiif", 0, 100 * (i + 1) - 1, 1, 0.0)
    shared += struct.pack("<II", len(alleles) << 16, ((2 + ("GQ" in keys)) << 24) | len(SAMPLES))
    shared += tstr(vid.encode()) + b"".join(tstr(a) for a in alleles) + b"\x00"
    indiv = b""
    # GT (dictionary index 1)
    gt = b""
    for p in per:
        g = p["GT"]
        phased = "|" in g
        parts = g.replace("|", "/").split("/")
        vals = [0 if a == "." else (int(a) + 1) << 1 for a in parts]
        if len(vals) == 2 and phased:
            vals[1] |= 1
        if len(vals) == 1:
            vals.append(EOV[1])
        gt += b"".join(struct.pack("<b", v) for v in vals)
    indiv += typed_int(1) + bytes([2 << 4 | 1]) + gt
    # GQ (index 2), LAD (index 3)
    if "GQ" in keys:
        indiv += typed_int(2) + bytes([1 << 4 | 1]) + b"".join(
            struct.pack("<b", MISS[1] if p.get("GQ", ".") == "." else int(p["GQ"])) for p in per)
    lads = [[None if x == "." else int(x) for x in p.get("LAD", ".").split(",")] for p in per]
    width = max(len(v) for v in lads)
    t = 2 if any(x is not None and x > 127 for v in lads for x in v) else 1
    lad = b""
    for v in lads:
        v = [MISS[t] if x is None else x for x in v] + [EOV[t]] * (width - len(v))
        lad += b"".join(struct.pack(PACK[t], x) for x in v)
    indiv += typed_int(3) + bytes([width << 4 | t]) + lad
    body += struct.pack("<II", len(shared), len(indiv)) + shared + indiv



def bgzf(data):
    out = b""
    for start in range(0, len(data), 0xff00):
        chunk = data[start:start + 0xff00]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        cdata = c.compress(chunk) + c.flush()
        out += b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", len(cdata) + 25)
        out += cdata + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    return out + bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


with open("lad.bcf", "wb") as f:
    f.write(bgzf(b"BCF\x02\x02" + struct.pack("<I", len(text)) + text + body))
