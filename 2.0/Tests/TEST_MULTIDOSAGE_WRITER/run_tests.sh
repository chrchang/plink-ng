#!/bin/bash

# Round-trip unit test for multiallelic dosages: the writer
# (SpgwAppendMultiallelicGenovecDphase16()) and the readers (PgrGetMDp(),
# PgrGet1D(), PgrGetMDCounts(), PgrGetRaw(), etc.), and the validator
# (PgrValidate()).  Builds its own driver from the
# pgenlib sources with plain flags, rather than reusing the build directory's
# objects (whose SIMD flags vary).  Set CXX to override the compiler.

set -exo pipefail

INC=../../include
${CXX:-c++} -std=c++17 -O1 -g -I$INC -o tmp_multidosage_writer_test multidosage_writer_test.cc $INC/pgenlib_read.cc $INC/pgenlib_write.cc $INC/pgenlib_misc.cc $INC/plink2_base.cc $INC/plink2_bits.cc $INC/plink2_simd.cc
./tmp_multidosage_writer_test

# plink2 --validate end to end (the validator workspace is sized in plink2.cc):
# the round-trip file must pass, and the implicit-dosage-phase file must fail.
$1/plink2 $2 $3 --pfile tmp_small --validate --out tmp_small_v
if $1/plink2 $2 $3 --pfile tmp_tiny --validate --out tmp_tiny_v; then
    echo "plink2 --validate accepted tmp_tiny.pgen"
    exit 1
fi
grep -q "no explicit dosage-phase value" tmp_tiny_v.log

# plink2 --freq and --missing dosage on multiallelic dosages (PgrGetMDCounts()):
# observation counts must reflect dosage-only samples.
$1/plink2 $2 $3 --pfile tmp_small --freq --missing variant-only vcols=chrom,nmissdosage --out tmp_small_c
awk 'NR > 1 {print $2"\t"$6/2}' tmp_small_c.afreq > tmp_small_freq_nm.txt
diff -q tmp_small_freq_nm.txt tmp_small_expected_nm.txt
awk 'NR > 1 {print $2"\t"(300 - $3)}' tmp_small_c.vmiss > tmp_small_vmiss_nm.txt
diff -q tmp_small_vmiss_nm.txt tmp_small_expected_nm.txt
# Without --freq, plink2 only needs ref vs. nonref counts here, and uses
# PgrGetDCounts() instead.
$1/plink2 $2 $3 --pfile tmp_small --missing variant-only vcols=chrom,nmissdosage --out tmp_small_m
awk 'NR > 1 {print $2"\t"(300 - $3)}' tmp_small_m.vmiss > tmp_small_vmiss2_nm.txt
diff -q tmp_small_vmiss2_nm.txt tmp_small_expected_nm.txt
