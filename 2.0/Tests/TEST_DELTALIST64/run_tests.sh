#!/bin/bash

# Unit test for the 64-bit deltalist encoding (multiallelic dosage tracks #5
# and #9).  Builds its own small driver from the pgenlib sources, with plain
# flags, rather than reusing the build directory's objects (whose SIMD flags
# vary).  Set CXX to override the compiler.

set -exo pipefail

INC=../../include
${CXX:-c++} -std=c++17 -O1 -g -I$INC -o tmp_deltalist64_test deltalist64_test.cc $INC/pgenlib_misc.cc $INC/plink2_base.cc $INC/plink2_bits.cc
./tmp_deltalist64_test
