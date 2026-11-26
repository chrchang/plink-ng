#!/usr/bin/env python3
# Writes small hand-built .pgen filesets (prefix tmp_) that the spec permits
# but --validate used to reject, plus equivalent "control" filesets in more
# common encodings.  Also writes two encodings (3 and 4 below) which the
# reader tolerates but --validate deliberately rejects; the spec is to be
# updated to prohibit them.  Usage: make_pgens.py

import random
import shutil
import struct


def side(name, n, alts):
    with open(name + '.pvar', 'w') as f:
        f.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for i, a in enumerate(alts):
            f.write('1\t%d\tv%d\tA\t%s\n' % (100 + i, i, a))
    with open(name + '.psam', 'w') as f:
        f.write('#IID\n' + ''.join('s%d\n' % i for i in range(n)))


def write_var(name, n, vrtypes, recs, alts=None, ctrl=0x47):
    """Mode 0x10, single vblock.  ctrl 0x47: 8-bit vrtypes, 4-byte lengths,
    all REF alleles trusted."""
    assert (ctrl & 15) == 7
    body = bytes(vrtypes) + b''.join(struct.pack('<I', len(r)) for r in recs)
    off = 12 + 8 + len(body)
    with open(name + '.pgen', 'wb') as f:
        f.write(b'\x6c\x1b\x10' + struct.pack('<II', len(recs), n) + bytes([ctrl]) + struct.pack('<Q', off) + body + b''.join(recs))
    side(name, n, alts or ['C'] * len(recs))


# 1. LD-compressed record whose base is all-hom-ref (type 5), vs. the same
#    genotypes with a 2-bit base.  1000 samples, so readers take their sparse
#    (difflist) paths.  LD difflist: sample 5 -> het.
N = 1000
ld = bytes([1]) + struct.pack('<H', 5) + bytes([1])
write_var('tmp_ldb5', N, [5, 2], [b'', ld])
write_var('tmp_ldb0', N, [0, 2], [bytes((N + 3) // 4), ld])

# 2. Fixed-width mode 0x02 whose .pvar lists a second ALT nobody carries.
n = 8
with open('tmp_fw_multi.pgen', 'wb') as f:
    f.write(b'\x6c\x1b\x02' + struct.pack('<II', 2, n) + bytes([0x40]) + bytes([0x01, 0x00]) * 2)
side('tmp_fw_multi', n, ['C,G', 'C'])

# 3. Multiallelic hardcall track with mode byte 0xff (both patch sets empty).
#    Must be rejected by --validate.
het0 = bytes([1, 0])
write_var('tmp_mh_ff', n, [0x08], [het0 + bytes([0xff])], ['C,G'])

# 4. Phased-dosage existence bitarray with no bits set.  Must be rejected by
#    --validate.  Dosage list: sample 0, dosage 1.0.
dl = bytes([1, 0]) + struct.pack('<H', 0x4000)
write_var('tmp_dphase_empty', n, [0xa0], [het0 + dl + bytes([0x00])])

# 5. Fused vrtype/length header modes 8-12 and 14 (pgenlib-defined; see
#    PgfiInitPhase1()), single sample, 3 vblocks with an odd variant count,
#    vs. the same records under header mode 7.
VCT = 140001
d16 = lambda v: struct.pack('<H', v)
# mode -> {code: (vrtype, record, ALT)}
T = {
    8: {0: (5, b'', 'C'), 1: (0, bytes([2]), 'C')},
    9: {0: (5, b'', 'C'), 1: (0, bytes([1]), 'C'),
        2: (8, bytes([1, 0xf0, 0x01]), 'C,G'),
        3: (8, bytes([1, 0xf0, 0x01, 0x01]), 'C,G,T')},
    10: {0: (5, b'', 'C'), 1: (0, bytes([2]), 'C'),
         2: (16, bytes([1, 0x02]), 'C'), 3: (16, bytes([1, 0x03, 0x01]), 'C')},
    11: {0: (5, b'', 'C'), 1: (0, bytes([2]), 'C'),
         2: (8, bytes([1, 0xf0, 0x01]), 'C,G'),
         3: (8, bytes([1, 0xf0, 0x01, 0x01]), 'C,G,T'),
         9: (16, bytes([1, 0x02]), 'C'),
         10: (24, bytes([1, 0xf0, 0x01, 0x00]), 'C,G')},
    12: {0: (5, b'', 'C'), 1: (0, bytes([1]), 'C'),
         2: (0x45, d16(0x0800), 'C'), 3: (0x40, bytes([1]) + d16(0x4000), 'C')},
    14: {0: (5, b'', 'C'), 1: (0, bytes([1]), 'C'),
         3: (16, bytes([1, 0x00]), 'C'),
         4: (0x45, d16(0x0800), 'C'), 5: (0x40, bytes([1]) + d16(0x4000), 'C'),
         12: (0xc5, d16(0x2000) + struct.pack('<h', 0), 'C'),
         13: (0xc0, bytes([1]) + d16(0x4000) + struct.pack('<h', 0x1000), 'C')},
}
LOG2_BITS = {8: 0, 9: 1, 10: 1, 11: 2, 12: 1, 14: 2}


def write_multiblock(name, ctrl, vblock_body, recs, alts):
    vblock_ct = (VCT + 65535) // 65536
    bodies = [vblock_body(i * 65536, min(VCT, (i + 1) * 65536)) for i in range(vblock_ct)]
    pos = 12 + 8 * vblock_ct + sum(len(b) for b in bodies)
    offs = []
    for i in range(vblock_ct):
        offs.append(pos)
        pos += sum(len(r) for r in recs[i * 65536:(i + 1) * 65536])
    with open(name + '.pgen', 'wb') as f:
        f.write(b'\x6c\x1b\x10' + struct.pack('<II', VCT, 1) + bytes([ctrl]) + struct.pack('<%dQ' % vblock_ct, *offs) + b''.join(bodies) + b''.join(recs))
    side(name, 1, alts)
    return 12 + 8 * vblock_ct + sum(len(b) for b in bodies[:-1])


random.seed(1)
for mode in sorted(T):
    tab = T[mode]
    codes = [random.choice(sorted(tab)) for _ in range(VCT)]
    vrtypes = [tab[c][0] for c in codes]
    recs = [tab[c][1] for c in codes]
    alts = [tab[c][2] for c in codes]
    lg = LOG2_BITS[mode]
    per_byte = 8 >> lg

    def fused(a, b):
        out = bytearray((b - a + per_byte - 1) // per_byte)
        for i in range(a, b):
            out[(i - a) // per_byte] |= codes[i] << (((i - a) % per_byte) << lg)
        return bytes(out)

    last_vblock_start = write_multiblock('tmp_m%d' % mode, 0x40 | mode, fused, recs, alts)
    write_multiblock('tmp_s%d' % mode, 0x47, lambda a, b: bytes(vrtypes[a:b]) + b''.join(struct.pack('<I', len(r)) for r in recs[a:b]), recs, alts)
    # Same as tmp_m<mode>, but with a nonzero trailing bit in the last fused
    # entry byte.
    b = bytearray(open('tmp_m%d.pgen' % mode, 'rb').read())
    b[last_vblock_start + ((VCT % 65536) >> (3 - lg))] |= 0x80
    with open('tmp_t%d.pgen' % mode, 'wb') as f:
        f.write(bytes(b))
    for ext in ('pvar', 'psam'):
        shutil.copy('tmp_m%d.%s' % (mode, ext), 'tmp_t%d.%s' % (mode, ext))
