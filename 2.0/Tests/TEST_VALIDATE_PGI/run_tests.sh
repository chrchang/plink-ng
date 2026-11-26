#!/bin/bash

# --validate on .pgen files whose header is in an external .pgi index file
# (storage modes 0x20/0x21).  PgrValidate() used to re-read the header from
# byte 11 of the .pgen itself, so valid files of this type were rejected with
# "vblock-start index is inconsistent with variant record length index".

set -exo pipefail

# 140000 variants = 3 variant blocks, so vblock offsets past the first (which
# only PgrValidate() checks) are exercised.
$1/plink2 $2 $3 --dummy 7 140000 0.05 --seed 1 --make-pgen --out tmp_data

# Rewrite as mode 0x20: .pgen = 3 magic bytes + variant records; .pgen.pgi =
# the original header with third byte 0x30 and vblock offsets rebased.
# Optional third argument: corrupt the last vblock offset.
cat > tmp_split.py <<'PYEOF'
import struct, sys
src, dst = sys.argv[1], sys.argv[2]
b = open(src + '.pgen', 'rb').read()
assert b[2] == 0x10
vct = struct.unpack('<I', b[3:7])[0]
vblock_ct = (vct + 65535) // 65536
offs = struct.unpack('<%dQ' % vblock_ct, b[12:12 + 8 * vblock_ct])
newoffs = [o - offs[0] + 3 for o in offs]
if len(sys.argv) > 3:
    newoffs[-1] += 1
pgi = b'\x6c\x1b\x30' + b[3:12] + struct.pack('<%dQ' % vblock_ct, *newoffs) + b[12 + 8 * vblock_ct:offs[0]]
open(dst + '.pgen', 'wb').write(b'\x6c\x1b\x20' + b[offs[0]:])
open(dst + '.pgen.pgi', 'wb').write(pgi)
PYEOF
python3 tmp_split.py tmp_data tmp_split
cp tmp_data.pvar tmp_split.pvar
cp tmp_data.psam tmp_split.psam

$1/plink2 $2 $3 --pfile tmp_split --validate --out tmp_split_v

# Genotypes are unchanged by the rewrite.
$1/plink2 $2 $3 --pfile tmp_data --export A --out tmp_data_e
$1/plink2 $2 $3 --pfile tmp_split --export A --out tmp_split_e
cmp tmp_data_e.raw tmp_split_e.raw

# Nondefault .pgi name via --pgi.
mkdir -p alt
mv tmp_split.pgen.pgi alt/tmp_other.pgi
$1/plink2 $2 $3 --pfile tmp_split --pgi alt/tmp_other.pgi --validate --out tmp_split_pgi_v
mv alt/tmp_other.pgi tmp_split.pgen.pgi

# A bad vblock offset in the .pgi is still caught.
python3 tmp_split.py tmp_data tmp_bad corrupt
cp tmp_data.pvar tmp_bad.pvar
cp tmp_data.psam tmp_bad.psam
if $1/plink2 $2 $3 --pfile tmp_bad --validate --out tmp_bad_v 2> tmp_err.txt; then
    echo "expected --validate to reject the corrupted .pgi vblock offset"
    exit 1
fi
grep -q "pgi file vblock-start index is inconsistent" tmp_bad_v.log

# --haps import writes an intermediate mode-0x20 .pgen + .pgi.
python3 -c "
import random
random.seed(1)
with open('tmp.haps', 'w') as f:
    for i in range(40):
        f.write('1 v%d %d A C %s\n' % (i, 1000 + i, ' '.join(random.choice('01') for _ in range(8))))
with open('tmp.sample', 'w') as f:
    f.write('ID_1 ID_2 missing\n0 0 0\n' + ''.join('s%d s%d 0\n' % (i, i) for i in range(4)))
"
$1/plink2 $2 $3 --haps tmp.haps --sample tmp.sample --validate --out tmp_haps_v
grep -q "done" tmp_haps_v.log
