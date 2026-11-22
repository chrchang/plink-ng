#!/usr/bin/env python3
# Writes BGEN v1.1 files (layout 1, zlib-compressed) used by run_tests.sh:
#   ok.bgen      two variants, matching ok.gen
#   short.bgen   the same, but the second variant's block decompresses to
#                6 bytes fewer than 6 * sample count
# plus ok.sample.
import struct
import zlib

N = 4
PROBS = [
    [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.5, 0.0)],
    [(0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.5, 0.5)],
]


def lstr16(b):
    return struct.pack("<H", len(b)) + b


def variant(k, probs, short):
    head = struct.pack("<I", N) + lstr16(b"snp%d" % k) + lstr16(b"rs%d" % k) + lstr16(b"01")
    head += struct.pack("<I", 1000 * k) + struct.pack("<I", 1) + b"A" + struct.pack("<I", 1) + b"G"
    raw = b"".join(struct.pack("<HHH", *(round(p * 32768) for p in sp)) for sp in probs)
    if short:
        raw = raw[:-6]
    comp = zlib.compress(raw)
    return head + struct.pack("<I", len(comp)) + comp


def write(name, short_second):
    # header block: length, variant count, sample count, magic, flags
    # (compressed, layout 1)
    header = struct.pack("<III", 20, len(PROBS), N) + b"bgen" + struct.pack("<I", 1 | (1 << 2))
    body = b"".join(variant(k + 1, p, short_second and k == 1) for k, p in enumerate(PROBS))
    with open(name, "wb") as f:
        f.write(struct.pack("<I", len(header)) + header + body)


write("ok.bgen", False)
write("short.bgen", True)
with open("ok.sample", "w") as f:
    f.write("ID_1 ID_2 missing\n0 0 0\n")
    for i in range(N):
        f.write("s%d s%d 0\n" % (i + 1, i + 1))
with open("ok.gen", "w") as f:
    for k, probs in enumerate(PROBS):
        f.write("1 rs%d %d A G " % (k + 1, 1000 * (k + 1)))
        f.write(" ".join("%g %g %g" % sp for sp in probs) + "\n")
