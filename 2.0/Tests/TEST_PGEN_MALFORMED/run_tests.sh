#!/bin/bash

# A .pgen record whose type byte says "1-bit" but whose body is not a valid
# 1-bit record must be rejected, not counted.
#
# The first body byte of a 1-bit record packs two genotype codes, and
# CountparseOnebitSubset() uses both of them to index a four-element
# genocounts[] array.  Before the accompanying fix, a corrupted byte indexed
# that array with a value of up to 66, and --freq still exited 0 with a report
# full of wrong numbers.

set -exo pipefail

$1/plink2 $2 $3 --dummy 400 600 --seed 1 --make-pgen --out tmp_data

# Sanity check: the unmodified fileset is fine.
$1/plink2 $2 $3 --pfile tmp_data --freq --out plink2_ok
test "$(grep -cv '^#' plink2_ok.afreq)" -eq 600

# Byte 92 is inside the record-type array, and setting it to 0x91 marks a
# record as 1-bit whose body then decodes to an out-of-range genotype code.
# (If a .pgen layout change moves that array, this offset is what needs
# updating; the scan that found it just tried every header byte.)
python3 -c "
b = bytearray(open('tmp_data.pgen', 'rb').read())
b[92] = 0x91
open('tmp_bad.pgen', 'wb').write(bytes(b))
"
cp tmp_data.pvar tmp_bad.pvar
cp tmp_data.psam tmp_bad.psam

# --validate has always caught this one.
if $1/plink2 $2 $3 --pfile tmp_bad --validate --out plink2_val 2> tmp_val_err.txt; then
    echo "expected --validate to reject the corrupted .pgen"
    exit 1
fi
grep -q "Invalid 1-bit genotype record" tmp_val_err.txt

# --freq has to reject it too, rather than reporting counts.
if $1/plink2 $2 $3 --pfile tmp_bad --freq --out plink2_bad 2> tmp_err.txt; then
    echo "expected --freq to reject the corrupted .pgen"
    exit 1
fi
grep -q "Failed to unpack" tmp_err.txt
test ! -e plink2_bad.afreq
