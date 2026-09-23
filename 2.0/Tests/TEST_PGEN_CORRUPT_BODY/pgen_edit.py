#!/usr/bin/env python3
# Rewrites a single-vblock .pgen (storage mode 0x10, 8-bit record types) with
# edits applied, recomputing the record-length index.
#   pgen_edit.py IN OUT EDIT...
# EDIT is one of
#   vrtype:V=X        set variant V's record type byte to X
#   byte:V:I=X        set byte I of variant V's record to X
#   u16:V:I=X         set little-endian uint16 at byte I of variant V's record
#   record:V=HEX      replace variant V's record
import struct, sys

src, dst = sys.argv[1], sys.argv[2]
d = open(src, "rb").read()
assert d[:3] == b"\x6c\x1b\x10", "not a mode 0x10 .pgen"
vct = struct.unpack_from("<I", d, 3)[0]
ctrl = d[11]
assert vct <= 65536 and (ctrl & 15) >= 4 and (ctrl & 0x30) == 0
lb = 1 + (ctrl & 3)
p = 20
vrtypes = bytearray(d[p:p + vct])
p += vct
lens = [int.from_bytes(d[p + i * lb:p + (i + 1) * lb], "little") for i in range(vct)]
p += vct * lb
rest_hdr = d[p:struct.unpack_from("<Q", d, 12)[0]]
recs = []
q = struct.unpack_from("<Q", d, 12)[0]
for ln in lens:
    recs.append(bytearray(d[q:q + ln]))
    q += ln
for e in sys.argv[3:]:
    kind, val = e.split("=")
    x = int(val, 0) if kind.split(":")[0] != "record" else None
    f = kind.split(":")
    v = int(f[1])
    if f[0] == "vrtype":
        vrtypes[v] = x
    elif f[0] == "byte":
        recs[v][int(f[2])] = x
    elif f[0] == "u16":
        struct.pack_into("<H", recs[v], int(f[2]), x)
    elif f[0] == "record":
        recs[v] = bytearray(bytes.fromhex(val))
    else:
        raise SystemExit("unknown edit " + e)
lens = [len(r) for r in recs]
assert max(lens) < (1 << (8 * lb))
hdr = bytearray(d[:12])
first = 12 + 8 + vct + vct * lb + len(rest_hdr)
hdr += struct.pack("<Q", first) + vrtypes + b"".join(ln.to_bytes(lb, "little") for ln in lens) + rest_hdr
open(dst, "wb").write(bytes(hdr) + b"".join(recs))
