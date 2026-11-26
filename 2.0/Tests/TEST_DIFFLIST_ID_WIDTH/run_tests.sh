#!/bin/bash

# Difflist group-start sample IDs are as wide as the number of bytes needed to
# represent the sample count N itself: uint8 iff N < 2^8, uint16 iff
# 2^8 <= N < 2^16, and so on.  (The spec used to say "uint8 iff N <= 2^8",
# which disagreed with pgenlib, and hence with every existing .pgen, at
# N = 2^8, 2^16 and 2^24.)  For N at either side of the 2^8 and 2^16
# boundaries, a difflist record encoded with the specified width must validate
# and decode correctly, and one encoded with the old spec's width must not.

set -exo pipefail

# One variant: record type 4 (difflist vs. hom-REF) with two entries,
# sample 3 = het and sample N-1 = hom-ALT.  id_byte_ct is the group-start
# sample ID width.
cat > tmp_make.py <<'PYEOF'
import struct, sys
name, n, id_byte_ct = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
def vint(x):
    out = bytearray()
    while x > 127:
        out.append((x & 127) | 128)
        x >>= 7
    out.append(x)
    return bytes(out)
rec = vint(2) + (3).to_bytes(id_byte_ct, 'little') + bytes([1 | (2 << 2)]) + vint(n - 1 - 3)
body = bytes([0x04]) + struct.pack('<I', len(rec))
hdr = b'\x6c\x1b\x10' + struct.pack('<II', 1, n) + bytes([0x47]) + struct.pack('<Q', 12 + 8 + len(body)) + body
open(name + '.pgen', 'wb').write(hdr + rec)
open(name + '.pvar', 'w').write('#CHROM\tPOS\tID\tREF\tALT\n1\t100\tv1\tA\tC\n')
with open(name + '.psam', 'w') as f:
    f.write('#IID\n' + ''.join('s%d\n' % i for i in range(n)))
PYEOF

# N, specified width, old spec's width
for case in "255 1 1" "256 2 1" "65535 2 2" "65536 3 2"; do
    set -- "$1" "$2" "$3" $case
    n=$4
    python3 tmp_make.py tmp_ok_$n $n $5
    $1/plink2 $2 $3 --pfile tmp_ok_$n --validate --out tmp_ok_${n}_v
    # Decoded genotypes: exactly sample 3 het and sample N-1 hom-ALT.
    $1/plink2 $2 $3 --pfile tmp_ok_$n --export vcf --out tmp_ok_${n}_e
    grep -v '^#' tmp_ok_${n}_e.vcf | tr '\t' '\n' | tail -n +10 | awk -v n=$n '
        (NR == 4 && $0 != "0/1") || (NR == n && $0 != "1/1") || (NR != 4 && NR != n && $0 != "0/0") {bad = 1}
        END {exit (bad || (NR != n))}'
    if [ $5 != $6 ]; then
        python3 tmp_make.py tmp_old_$n $n $6
        if $1/plink2 $2 $3 --pfile tmp_old_$n --validate --out tmp_old_${n}_v; then
            echo "expected --validate to reject the old-spec width at N=$n"
            exit 1
        fi
    fi
done
