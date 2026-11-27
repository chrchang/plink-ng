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

# --make-pgen on multiallelic dosages.  Each output must pass --validate and
# match its input semantically (samples matched by IID, variants by ID,
# alleles by string; see CompareFilesets()).
P="$1/plink2 $2 $3"
C="./tmp_multidosage_writer_test compare"
mp_check() {
    $P --pfile $2 --validate --out ${2}_v
    $C $1 $2 $3
}
awk 'NR > 1 && (NR % 3 != 0) {print $1}' tmp_small.psam > tmp_keep.txt
awk '{a[NR] = $1} END {for (i = NR; i; --i) print a[i]}' tmp_keep.txt > tmp_order.txt
# REF := ALT2, and ALT1 := ALT3 where present
awk '!/^#/ && (split($5, a, ",") >= 2) {print $3"\t"a[2]}' tmp_small.pvar > tmp_ref.txt
awk '!/^#/ && (split($5, a, ",") >= 3) {print $3"\t"a[3]}' tmp_small.pvar > tmp_alt1.txt
# default multithreaded writer
$P --pfile tmp_small --make-pgen --out tmp_mp1
mp_check tmp_small tmp_mp1
# sample subset, and subset + reorder
$P --pfile tmp_small --keep tmp_keep.txt --make-pgen --out tmp_mp2
mp_check tmp_small tmp_mp2
$P --pfile tmp_small --keep tmp_keep.txt --indiv-sort f tmp_order.txt --make-pgen --out tmp_mp3
mp_check tmp_small tmp_mp3
# allele permutation
$P --pfile tmp_small --ref-allele force tmp_ref.txt --alt1-allele force tmp_alt1.txt --make-pgen --out tmp_mp4
mp_check tmp_small tmp_mp4
# single-threaded MakePgenRobust() path, with everything at once
$P --pfile tmp_small --sort-vars --keep tmp_keep.txt --indiv-sort f tmp_order.txt --ref-allele force tmp_ref.txt --make-pgen --out tmp_mp5
mp_check tmp_small tmp_mp5
# trim-alts, erase-phase, erase-dosage
$P --pfile tmp_small --make-pgen trim-alts --out tmp_mp6
mp_check tmp_small tmp_mp6
$P --pfile tmp_small --make-pgen erase-phase --out tmp_mp7
mp_check tmp_small tmp_mp7 strip-phase
$P --pfile tmp_small --make-pgen erase-dosage --out tmp_mp8
mp_check tmp_small tmp_mp8 strip-dosage
# --set-mixed-mt-missing: hets (including altx/alty) become missing, and lose
# their dosages
cp tmp_small.pgen tmp_mt.pgen
cp tmp_small.psam tmp_mt.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "MT"; print}' tmp_small.pvar > tmp_mt.pvar
$P --pfile tmp_mt --set-mixed-mt-missing --make-pgen --out tmp_mp9
mp_check tmp_mt tmp_mp9 hets-missing
# large rare-alt lists, multithreaded; and subsetted
$P --pfile tmp_big --threads 2 --make-pgen --out tmp_mpb1
mp_check tmp_big tmp_mpb1
awk 'NR > 1 && (NR % 2) {print $1}' tmp_big.psam > tmp_keep_big.txt
$P --pfile tmp_big --threads 2 --keep tmp_keep_big.txt --make-pgen --out tmp_mpb2
mp_check tmp_big tmp_mpb2
# fixed-width input: ALT1 := ALT2 moves every alt dosage into a rare-alt list
awk '!/^#/ {split($5, a, ","); print $3"\t"a[2]}' tmp_fixed.pvar > tmp_fixed_alt1.txt
$P --pfile tmp_fixed --make-pgen --out tmp_mpf1
mp_check tmp_fixed tmp_mpf1
$P --pfile tmp_fixed --alt1-allele force tmp_fixed_alt1.txt --make-pgen --out tmp_mpf2
mp_check tmp_fixed tmp_mpf2
# not supported yet
if $P --pfile tmp_small --hard-call-threshold 0.1 --make-pgen --out tmp_mpx; then
    echo "--hard-call-threshold unexpectedly accepted multiallelic dosages"
    exit 1
fi
grep -q "don't support multiallelic dosages yet" tmp_mpx.log
