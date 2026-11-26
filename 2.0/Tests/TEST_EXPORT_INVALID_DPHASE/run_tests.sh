#!/bin/bash

# PglDosagesAreInvalid() (called on dosages read for --export vcf
# vcf-dosage=HDS) used to accept a phased-dosage delta of -32768 when the
# dosage was 0 or 32768: both per-haplotype values wrap around in uint16 and
# pass the <= 32768 range check.  Such a record must now be rejected instead
# of being exported as garbage.

set -exo pipefail

# One variant, 8 samples, record type 0xc0 (unconditional phased dosage).
# Hardcalls are all hom-ALT except sample 1 (hom-REF); dosages 32768 except
# sample 1 (0).  The delta for sample 0 (or 1) is the last argument.
cat > tmp_make.py <<'PYEOF'
import struct, sys
name, s0_delta, s1_delta = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
geno = bytes([0b10101010 & ~0b1100, 0b10101010])
dos = [32768, 0] + [32768] * 6
delta = [s0_delta, s1_delta] + [0] * 6
rec = geno + b''.join(struct.pack('<H', v) for v in dos) + b''.join(struct.pack('<h', v) for v in delta)
body = bytes([0xc0]) + struct.pack('<I', len(rec))
hdr = b'\x6c\x1b\x10' + struct.pack('<II', 1, 8) + bytes([0x47]) + struct.pack('<Q', 12 + 8 + len(body)) + body
open(name + '.pgen', 'wb').write(hdr + rec)
open(name + '.pvar', 'w').write('#CHROM\tPOS\tID\tREF\tALT\n1\t100\tv1\tA\tC\n')
open(name + '.psam', 'w').write('#IID\n' + ''.join('s%d\n' % i for i in range(8)))
PYEOF
python3 tmp_make.py tmp_ok 0 0
python3 tmp_make.py tmp_bad_hi -32768 0
python3 tmp_make.py tmp_bad_lo 0 -32768

$1/plink2 $2 $3 --pfile tmp_ok --export vcf vcf-dosage=HDS --out tmp_ok_e
for f in tmp_bad_hi tmp_bad_lo; do
    if $1/plink2 $2 $3 --pfile $f --export vcf vcf-dosage=HDS --out ${f}_e; then
        echo "expected --export vcf vcf-dosage=HDS to reject $f"
        exit 1
    fi
done
