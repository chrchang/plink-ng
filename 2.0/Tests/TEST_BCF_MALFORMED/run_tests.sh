#!/bin/bash

# A key without '=' in a contig/FILTER/INFO/FORMAT line of the BCF text
# header must be rejected as malformed.  BcfHeaderLineIdxCheck() used to
# search for the '=' without stopping at the end of the line, so on the last
# such line it ran past the end of the header block.

set -exo pipefail

$1/plink2 $2 $3 --dummy 8 20 --seed 1 --export bcf --out tmp_data

# Sanity check: the unmodified file imports.
$1/plink2 $2 $3 --bcf tmp_data.bcf --out tmp_ok

# Replace ",IDX=<n>" on the last structured header line with a key that has no
# '=', then write the file back out as BGZF.
python3 -c "
import gzip, struct, zlib
raw = gzip.decompress(open('tmp_data.bcf', 'rb').read())
hl = struct.unpack('<I', raw[5:9])[0]
lines = raw[9:9 + hl].split(b'\n')
k = max(i for i, l in enumerate(lines) if l.startswith(b'##') and (b',IDX=' in l))
lines[k] = lines[k][:lines[k].rfind(b',IDX=')] + b',noequals>'
h = b'\n'.join(lines)
raw = raw[:5] + struct.pack('<I', len(h)) + h + raw[9 + hl:]
out = b''
for i in range(0, len(raw), 60000):
    chunk = raw[i:i + 60000]
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = c.compress(chunk) + c.flush()
    out += b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00' + struct.pack('<H', len(cdata) + 25) + cdata + struct.pack('<II', zlib.crc32(chunk), len(chunk))
out += bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
open('tmp_bad.bcf', 'wb').write(out)
"

if $1/plink2 $2 $3 --bcf tmp_bad.bcf --out tmp_bad 2> tmp_err.txt; then
    echo "expected --bcf to reject the malformed header line"
    exit 1
fi
grep -q "BCF text header block is malformed" tmp_err.txt
