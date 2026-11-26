#!/usr/bin/env python3
# Writes hand-built .pgen filesets (prefix tmp_) for TEST_VALIDATE_GAPS:
# tmp_ok_* must pass --validate, tmp_bad_* must fail.  Usage: make_pgens.py

import struct

N = 8  # samples
u16 = lambda v: struct.pack('<H', v)
i16 = lambda v: struct.pack('<h', v)


def vint(x):
    out = bytearray()
    while True:
        b = x & 127
        x >>= 7
        if x:
            out.append(b | 128)
        else:
            out.append(b)
            return bytes(out)


def side(name, alts):
    with open(name + '.pvar', 'w') as f:
        f.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for i, a in enumerate(alts):
            f.write('1\t%d\tv%d\tA\t%s\n' % (100 + i, i, a))
    with open(name + '.psam', 'w') as f:
        f.write('#IID\n' + ''.join('s%d\n' % i for i in range(N)))


def write(name, vrtypes, recs, alts=None, header_exts=(), footer_exts=(),
          gap=b'', trailing=b'', split=False, footer_gap=b'', ext_flags=None):
    """Variable-width .pgen (8-bit vrtypes, 4-byte lengths, all REF trusted).
    header_exts/footer_exts: lists of (type_idx, body).  Extensions select
    mode 0x11/0x21; split selects 0x20/0x21 (header in .pgen.pgi)."""
    vct = len(recs)
    has_ext = bool(header_exts or footer_exts or ext_flags)
    mode = (0x20 if split else 0x10) | (1 if has_ext else 0)
    body = bytes(vrtypes) + b''.join(struct.pack('<I', len(r)) for r in recs)
    hdr_core = struct.pack('<II', vct, N) + bytes([0x47])
    records = b''.join(recs)
    hflags = sum(1 << t for t, _ in header_exts)
    fflags = sum(1 << t for t, _ in footer_exts)
    if ext_flags:
        hflags, fflags = ext_flags
    ext_part_no_fpos = b''
    hsizes = b''.join(vint(len(b)) for _, b in sorted(header_exts))
    hbodies = b''.join(b for _, b in sorted(header_exts))
    footer = b''.join(vint(len(b)) for _, b in sorted(footer_exts)) + b''.join(b for _, b in sorted(footer_exts))
    # header length without the record-start offset
    ext_len = 0
    if has_ext:
        ext_len = len(vint(hflags)) + len(vint(fflags)) + (8 if fflags else 0) + len(hsizes) + len(hbodies)
    if split:
        records_start = 3 + len(gap)
        pgen_prefix = bytes([0x6c, 0x1b, mode]) + gap
    else:
        records_start = 3 + len(hdr_core) + 8 + len(body) + ext_len + len(gap)
    footer_fpos = records_start + len(records) + len(footer_gap)
    ext_part = b''
    if has_ext:
        ext_part = vint(hflags) + vint(fflags) + (struct.pack('<Q', footer_fpos) if fflags else b'') + hsizes + hbodies
    header = hdr_core + struct.pack('<Q', records_start) + body + ext_part
    if split:
        with open(name + '.pgen', 'wb') as f:
            f.write(pgen_prefix + records + footer_gap + footer + trailing)
        with open(name + '.pgen.pgi', 'wb') as f:
            f.write(bytes([0x6c, 0x1b, mode + 0x10]) + header)
    else:
        with open(name + '.pgen', 'wb') as f:
            f.write(bytes([0x6c, 0x1b, mode]) + header + gap + records + footer_gap + footer + trailing)
    side(name, alts or ['C'] * vct)


def write_fixed(name, mode, recs):
    with open(name + '.pgen', 'wb') as f:
        f.write(bytes([0x6c, 0x1b, mode]) + struct.pack('<II', len(recs), N) + bytes([0x40]) + b''.join(recs))
    side(name, ['C'] * len(recs))


def geno(cats):
    out = bytearray((len(cats) + 3) // 4)
    for i, c in enumerate(cats):
        out[i // 4] |= c << (2 * (i % 4))
    return bytes(out)


CATS = [0, 1, 2, 3, 0, 1, 2, 3]
G = geno(CATS)
# Unconditional dosages at the edges of the 0.5 consistency window, plus
# missing dosage with missing hardcall, and a dosage with missing hardcall.
UNCOND_DOS = [8192, 8192, 24576, 65535, 0, 24576, 32768, 12345]
UNCOND = G + b''.join(u16(v) for v in UNCOND_DOS)
# Dosage list: sample 0 (hardcall 0), dosage 0.5.
DLIST = G + bytes([1, 0]) + u16(8192)
# Dosage bitarray: samples 1, 2 (hardcalls 1, 2).
DBITS = G + bytes([0x06]) + u16(16384 + 8192) + u16(16384 + 8192)
# Unconditional phased dosage; missing entries have delta -32768.
PH_DOS = [0, 16384, 32768, 65535, 0, 16384, 32768, 65535]
PH_DELTA = [0, 16384, 0, -32768, 0, -16384, 0, -32768]
UNCOND_PH = G + b''.join(u16(v) for v in PH_DOS) + b''.join(i16(v) for v in PH_DELTA)
# Dosage list + phased dosage: sample 1 (hardcall 1), dosage 1, left=1 right=0.
DLIST_PH = G + bytes([1, 1]) + u16(16384) + bytes([0x01]) + i16(16384)
OK_VRTYPES = [0x40, 0x20, 0x60, 0xc0, 0xa0]
OK_RECS = [UNCOND, DLIST, DBITS, UNCOND_PH, DLIST_PH]

write('tmp_ok_dosage', OK_VRTYPES, OK_RECS)
# Same records with header extensions (writer ID, padding, developer flag
# 0x1) and a footer extension (developer flag 0x1), in both layouts.
EXTS = dict(header_exts=[(0, b'dev'), (1, 'plïnk test'.encode('utf-8')), (2, bytes(37))],
            footer_exts=[(0, b'footer data')])
write('tmp_ok_ext', OK_VRTYPES, OK_RECS, **EXTS)
write('tmp_ok_ext_pgi', OK_VRTYPES, OK_RECS, split=True, **EXTS)
write('tmp_ok_pgi', OK_VRTYPES, OK_RECS, split=True)
write_fixed('tmp_ok_fw3', 0x03, [UNCOND])
write_fixed('tmp_ok_fw4', 0x04, [UNCOND_PH])


def with_dos(dos):
    return G + b''.join(u16(v) for v in dos)


def ph(dos, delta):
    return G + b''.join(u16(v) for v in dos) + b''.join(i16(v) for v in delta)


BAD = {
    # name: (vrtypes, recs, kwargs)
    'incons_uncond': ([0x40], [with_dos([8193] + UNCOND_DOS[1:])], {}),
    'incons_uncond_hi': ([0x40], [with_dos(UNCOND_DOS[:2] + [24575] + UNCOND_DOS[3:])], {}),
    'missing_uncond': ([0x40], [with_dos([65535] + UNCOND_DOS[1:])], {}),
    'incons_list': ([0x20], [G + bytes([1, 0]) + u16(8193)], {}),
    'incons_bits': ([0x60], [G + bytes([0x06]) + u16(8191) + u16(16384)], {}),
    'range_uncond': ([0x40], [with_dos(UNCOND_DOS[:7] + [40000])], {}),
    'dphase_m32768_hi': ([0xc0], [ph(PH_DOS, [0, 16384, -32768] + PH_DELTA[3:])], {}),
    'dphase_m32768_lo': ([0xc0], [ph(PH_DOS, [-32768] + PH_DELTA[1:])], {}),
    'dphase_m32768_list': ([0xa0], [G + bytes([1, 2]) + u16(32768) + bytes([0x01]) + i16(-32768)], {}),
    'mh_no_mode_byte': ([0x08], [G], {'alts': ['C,G']}),
    'gap_header': (OK_VRTYPES, OK_RECS, {'gap': b'\x00'}),
    'gap_header_ext': (OK_VRTYPES, OK_RECS, dict(EXTS, gap=b'\x00')),
    'gap_pgi': (OK_VRTYPES, OK_RECS, {'split': True, 'gap': b'\x00'}),
    'trailing': (OK_VRTYPES, OK_RECS, {'trailing': b'\x00'}),
    'trailing_footer': (OK_VRTYPES, OK_RECS, dict(EXTS, trailing=b'\x00')),
    'footer_gap': (OK_VRTYPES, OK_RECS, dict(EXTS, footer_gap=b'\x00')),
    'padding_nonzero': (OK_VRTYPES, OK_RECS, {'header_exts': [(2, b'\x00\x01')]}),
    'writer_id_utf8': (OK_VRTYPES, OK_RECS, {'header_exts': [(1, b'\xc0\xaf')]}),
    'writer_id_truncated': (OK_VRTYPES, OK_RECS, {'header_exts': [(1, b'ab\xe2\x82')]}),
    'flag_bit_256': (OK_VRTYPES, OK_RECS, {'ext_flags': (1 << 256, 0)}),
}
for name, (vrtypes, recs, kwargs) in BAD.items():
    write('tmp_bad_' + name, vrtypes, recs, **kwargs)

# .pgi with a trailing byte.
write('tmp_bad_pgi_trailing', OK_VRTYPES, OK_RECS, split=True)
with open('tmp_bad_pgi_trailing.pgen.pgi', 'ab') as f:
    f.write(b'\x00')
# Fixed-width dosage modes.
write_fixed('tmp_bad_fw3_incons', 0x03, [with_dos([8193] + UNCOND_DOS[1:])])
write_fixed('tmp_bad_fw3_missing', 0x03, [with_dos([65535] + UNCOND_DOS[1:])])
write_fixed('tmp_bad_fw3_range', 0x03, [with_dos(UNCOND_DOS[:7] + [40000])])
write_fixed('tmp_bad_fw4_dphase', 0x04, [ph(PH_DOS, [-32768] + PH_DELTA[1:])])
