#!/usr/bin/env python3
# Writes the .vcf/.bcf inputs used by run_tests.sh.
import struct
import zlib

SAMPLES = "abcdefgh"
ALTS12 = ["C", "G", "T", "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC"]
# (id, alts, FORMAT, per-sample fields)
# LAD only.
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
# AD, alone and alongside LAD (which takes precedence, whichever comes first).
AD_RECORDS = [
    ("a1", ["C"], "GT:AD",
     ["0/1:5,6", "0/0:9,0", "1:0,5", "0:4,0", "0/0:.,.", "0/1", "1/1:0,10", "./.:."]),
    ("a2", ["C"], "GT:AD:LAD",
     ["0/1:150,150:5,4", "0/0:3,0:12", "0/1:5,5:.", "0/0:20,0", "1:2,2:6", "0/0:30,0:10", "1/1:0,8:0,11", "0/0:.,.:."]),
    ("a3", ["C"], "GT:LAD:AD",
     ["0/1:5,6:1,1", "0/0:9:50,0", "0/0:10", "1:4:0,9", "0/0:.:20,0", "1/1:0,300:0,1", "0/1:3,3:3,3", "./.:.:.,."]),
    ("a4", ["C"], "GT:GQ:AD",
     ["0/1:30:5,6", "0/0:15:20,0", "0/0:30:9,0", "0/0:30:.,.", "0/0:30", "1:30:0,5", "0/1:.:6,6", "0/0:25:300,0"]),
]
# Header declares AD but not LAD, so the undeclared LAD is ignored.  VCF only.
AD_ONLY_RECORDS = [
    ("z1", ["C"], "GT:LAD:AD",
     ["0/1:20:3,3", "0/0:2:10,0", "0/0:30", "1:1:0,5", "0/0:9:9,1", "0/1:.:.,.", "1/1:0:0,12", "./."]),
]
GQ_LINE = '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">'
LAD_LINE = '##FORMAT=<ID=LAD,Number=.,Type=Integer,Description="l">'
AD_LINE = '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">'
HEADER = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=PASS,Description="ok">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">',
    GQ_LINE,
    LAD_LINE,
    "##contig=<ID=1,length=100000>",
]
AD_HEADER = HEADER[:5] + [AD_LINE] + HEADER[5:]
AD_ONLY_HEADER = HEADER[:4] + [AD_LINE] + HEADER[5:]
NO_DEPTH_HEADER = HEADER[:4] + HEADER[5:]
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES)


def write_vcf(fname, header, records):
    with open(fname, "w") as f:
        f.write("\n".join(header + [COLS]) + "\n")
        for i, (vid, alts, fmt, fields) in enumerate(records):
            f.write("\t".join(["1", str(100 * (i + 1)), vid, "A", ",".join(alts), ".", "PASS", ".", fmt] + fields) + "\n")


write_vcf("lad.vcf", HEADER, RECORDS + VCF_ONLY)
write_vcf("ad.vcf", AD_HEADER, AD_RECORDS)
write_vcf("adonly.vcf", AD_ONLY_HEADER, AD_ONLY_RECORDS)
write_vcf("nodepth.vcf", NO_DEPTH_HEADER, [("n1", ["C"], "GT:GQ", ["0/1:30"] * 8)])

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


def header_dict(header):
    # BCF string dictionary: PASS first, then FILTER/INFO/FORMAT IDs in order
    keys = ["PASS"]
    for line in header:
        if line.startswith(("##FORMAT=<ID=", "##INFO=<ID=", "##FILTER=<ID=")):
            k = line.split("ID=")[1].split(",")[0]
            if k not in keys:
                keys.append(k)
    return {k: i for i, k in enumerate(keys)}


def write_bcf(fname, header, records):
    sidx = header_dict(header)
    text = ("\n".join(header + [COLS]) + "\n").encode() + b"\0"
    body = b""
    for i, (vid, alts, fmt, fields) in enumerate(records):
        keys = fmt.split(":")
        per = [dict(zip(keys, s.split(":"))) for s in fields]
        alleles = [b"A"] + [a.encode() for a in alts]
        shared = struct.pack("<iiif", 0, 100 * (i + 1) - 1, 1, 0.0)
        shared += struct.pack("<II", len(alleles) << 16, (len(keys) << 24) | len(SAMPLES))
        shared += tstr(vid.encode()) + b"".join(tstr(a) for a in alleles) + b"\x00"
        indiv = b""
        # FORMAT vectors in FORMAT-column order
        for key in keys:
            if key == "GT":
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
                indiv += typed_int(sidx["GT"]) + bytes([2 << 4 | 1]) + gt
            elif key == "GQ":
                indiv += typed_int(sidx["GQ"]) + bytes([1 << 4 | 1]) + b"".join(
                    struct.pack("<b", MISS[1] if p.get("GQ", ".") == "." else int(p["GQ"])) for p in per)
            else:
                # integer vector (LAD or AD); a dropped field is a single
                # missing value
                vecs = [[None if x == "." else int(x) for x in p.get(key, ".").split(",")] for p in per]
                width = max(len(v) for v in vecs)
                t = 2 if any(x is not None and x > 127 for v in vecs for x in v) else 1
                vec = b""
                for v in vecs:
                    v = [MISS[t] if x is None else x for x in v] + [EOV[t]] * (width - len(v))
                    vec += b"".join(struct.pack(PACK[t], x) for x in v)
                indiv += typed_int(sidx[key]) + bytes([width << 4 | t]) + vec
        body += struct.pack("<II", len(shared), len(indiv)) + shared + indiv
    with open(fname, "wb") as f:
        f.write(bgzf(b"BCF\x02\x02" + struct.pack("<I", len(text)) + text + body))


def bgzf(data):
    out = b""
    for start in range(0, len(data), 0xff00):
        chunk = data[start:start + 0xff00]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        cdata = c.compress(chunk) + c.flush()
        out += b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", len(cdata) + 25)
        out += cdata + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    return out + bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


write_bcf("lad.bcf", HEADER, RECORDS)
write_bcf("ad.bcf", AD_HEADER, AD_RECORDS)
