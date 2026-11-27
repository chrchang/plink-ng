#!/bin/bash

# Unit test for the multiallelic dosage writer
# (SpgwAppendMultiallelicGenovecDphase16()).  Builds its own driver from the
# pgenlib sources with plain flags, rather than reusing the build directory's
# objects (whose SIMD flags vary).  Set CXX to override the compiler.

set -exo pipefail

INC=../../include
${CXX:-c++} -std=c++17 -O1 -g -I$INC -o tmp_multidosage_writer_test multidosage_writer_test.cc $INC/pgenlib_write.cc $INC/pgenlib_misc.cc $INC/plink2_base.cc $INC/plink2_bits.cc
./tmp_multidosage_writer_test
