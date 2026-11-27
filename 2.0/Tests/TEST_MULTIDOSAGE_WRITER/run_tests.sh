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
    $C "$@"
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

# --hard-call-threshold, 'fill-missing-from-dosage', and
# --dosage-erase-threshold on multiallelic dosages (thresholds passed to the
# comparison in 16384ths of a copy: 0.1 -> 1638, 0.2 -> 3276).
$P --pfile tmp_small --hard-call-threshold 0.1 --make-pgen --out tmp_hc1
mp_check tmp_small tmp_hc1 hct=1638
$P --pfile tmp_small --make-pgen fill-missing-from-dosage --out tmp_hc2
mp_check tmp_small tmp_hc2 fill
$P --pfile tmp_small --hard-call-threshold 0.1 --make-pgen fill-missing-from-dosage --out tmp_hc3
mp_check tmp_small tmp_hc3 hct=1638 fill
$P --pfile tmp_small --dosage-erase-threshold 0.2 --make-pgen --out tmp_hc4
mp_check tmp_small tmp_hc4 erase=3276
# (--dosage-erase-threshold can't exceed --hard-call-threshold)
$P --pfile tmp_small --keep tmp_keep.txt --hard-call-threshold 0.2 --dosage-erase-threshold 0.1 --make-pgen fill-missing-from-dosage --out tmp_hc5
mp_check tmp_small tmp_hc5 hct=3276 fill erase=1638
$P --pfile tmp_small --hard-call-threshold 0.1 --make-pgen erase-dosage --out tmp_hc6
mp_check tmp_small tmp_hc6 hct=1638 strip-dosage
# (fill is left out here since the comparison's brute-force nearest-hardcall
# search is slow with 255 alleles)
$P --pfile tmp_big --threads 2 --hard-call-threshold 0.1 --make-pgen --out tmp_hcb
mp_check tmp_big tmp_hcb hct=1638

# --flip + --flip-subset on multiallelic dosages.  tmp_flip is tmp_small with
# A/C/G/T/N alleles wherever there are at most 5 alleles (in several orders);
# the 20- and 255-allele variants keep their non-SNP alleles, and are skipped
# by 'permissive'.
cp tmp_small.pgen tmp_flip.pgen
cp tmp_small.psam tmp_flip.psam
awk 'BEGIN {
  FS = "\t"; OFS = "\t";
  split("A,T C,G T,A G,C", s2, " ");
  split("A,T,N N,G,C C,N,G T,A,N", s3, " ");
  split("A,T,C,G C,A,G,T G,C,T,A T,G,A,C", s4, " ");
  split("A,C,G,T,N N,T,G,C,A G,A,N,C,T C,G,A,T,N", s5, " ");
}
/^#/ {print; next}
{
  n = split($5, a, ",") + 1;
  k = (int(vidx / 6) % 4) + 1;
  ++vidx;
  lab = "";
  if (n == 2) lab = s2[k];
  else if (n == 3) lab = s3[k];
  else if (n == 4) lab = s4[k];
  else if (n == 5) lab = s5[k];
  if (lab != "") {
    $4 = substr(lab, 1, 1);
    $5 = substr(lab, 3);
  }
  print;
}' tmp_small.pvar > tmp_flip.pvar
awk '!/^#/ && (NR % 7) {print $3}' tmp_flip.pvar > tmp_flip_vars.txt
awk 'NR == 1 || (NR % 2) {print $1}' tmp_flip.psam > tmp_flip_samples.txt
awk '!/^#/ && (split($5, a, ",") >= 2) {print $3"\t"a[2]}' tmp_flip.pvar > tmp_flip_ref.txt
FL="--flip tmp_flip_vars.txt permissive --flip-subset tmp_flip_samples.txt"
FC="flip-vars=tmp_flip_vars.txt flip-samples=tmp_flip_samples.txt"
$P --pfile tmp_flip $FL --make-pgen --out tmp_fl1
mp_check tmp_flip tmp_fl1 $FC
# sample subset + reorder + allele permutation, multithreaded and
# single-threaded (MakePgenRobust())
$P --pfile tmp_flip $FL --keep tmp_keep.txt --indiv-sort f tmp_order.txt --ref-allele force tmp_flip_ref.txt --make-pgen --out tmp_fl2
mp_check tmp_flip tmp_fl2 $FC
$P --pfile tmp_flip $FL --sort-vars --keep tmp_keep.txt --ref-allele force tmp_flip_ref.txt --make-pgen --out tmp_fl3
mp_check tmp_flip tmp_fl3 $FC
# hardcalls only: --flip-subset can turn 0/1 hets into rare-alt hets, which
# previously overflowed the in-place patch buffers when no sample subsetting
# was needed
$P --pfile tmp_flip --make-pgen erase-dosage --out tmp_flipnd
$P --pfile tmp_flipnd $FL --make-pgen --out tmp_fl0
mp_check tmp_flipnd tmp_fl0 $FC
# after --dosage-erase-threshold / --hard-call-threshold erasures
$P --pfile tmp_flip $FL --hard-call-threshold 0.2 --dosage-erase-threshold 0.1 --make-pgen --out tmp_fl4
mp_check tmp_flip tmp_fl4 hct=3276 erase=1638 $FC
# after --set-mixed-mt-missing erasures
cp tmp_flip.pgen tmp_flipmt.pgen
cp tmp_flip.psam tmp_flipmt.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "MT"; print}' tmp_flip.pvar > tmp_flipmt.pvar
$P --pfile tmp_flipmt $FL --set-mixed-mt-missing --make-pgen --out tmp_fl5
mp_check tmp_flipmt tmp_fl5 hets-missing $FC
# tmp_flipfail: sample s0 has a missing hardcall, and dosages of 1 copy each
# for REF and ALT1.  With alleles C/G/A, flipping s0 is fine (the A has no
# complement, but s0 has no A dosage); with A/C/G, it's an error.
echo "v0" > tmp_flipfail_vars.txt
printf '#IID\ns0\n' > tmp_flipfail_samples.txt
FFL="--flip tmp_flipfail_vars.txt --flip-subset tmp_flipfail_samples.txt"
mv tmp_flipfail.pvar tmp_flipfail_orig.pvar
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$4 = "C"; $5 = "G,A"; print}' tmp_flipfail_orig.pvar > tmp_flipfail.pvar
$P --pfile tmp_flipfail $FFL --make-pgen --out tmp_ff1
mp_check tmp_flipfail tmp_ff1 flip-vars=tmp_flipfail_vars.txt flip-samples=tmp_flipfail_samples.txt
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$4 = "A"; $5 = "C,G"; print}' tmp_flipfail_orig.pvar > tmp_flipfail.pvar
if $P --pfile tmp_flipfail $FFL --make-pgen --out tmp_ff2; then
    echo "--flip-subset accepted a dosage on an allele without a complement"
    exit 1
fi
grep -q "Cannot flip" tmp_ff2.log

# Multiallelic split ('varid-dup' keeps the original IDs, which the
# comparison matches on): each split variant gets its ALT allele's dosage and
# dosage-phase value, with the other ALT alleles lumped with REF.
$P --pfile tmp_small --make-pgen multiallelics=- varid-dup --out tmp_sp1
mp_check tmp_small tmp_sp1 split
$P --pfile tmp_small --keep tmp_keep.txt --indiv-sort f tmp_order.txt --make-pgen multiallelics=- varid-dup --out tmp_sp2
mp_check tmp_small tmp_sp2 split
$P --pfile tmp_small --make-pgen multiallelics=- varid-dup erase-phase --out tmp_sp3
mp_check tmp_small tmp_sp3 split strip-phase
# '-snps': the A/C/G/T/N variants are split, and the others keep their
# multiallelic dosage records
$P --pfile tmp_flip --make-pgen multiallelics=-snps varid-dup --out tmp_sp4
mp_check tmp_flip tmp_sp4 split
# sample indexes past 2^16, and (with --memory 640, 623-variant write blocks)
# 255-allele splits interrupted by write block boundaries
awk 'NR > 65001 {print $1}' tmp_big.psam > tmp_keep_bigsplit.txt
$P --pfile tmp_big --keep tmp_keep_bigsplit.txt --memory 640 --make-pgen multiallelics=- varid-dup --out tmp_spb
mp_check tmp_big tmp_spb split

# Variant-join (multiallelics=+...) with dosages.
# Joining the split filesets must reproduce the originals exactly (the split
# variants share their original's ID, which the joined variant keeps).
$P --pfile tmp_sp1 --make-pgen multiallelics=+any --out tmp_jn1
mp_check tmp_small tmp_jn1
$P --pfile tmp_sp1 --keep tmp_keep.txt --indiv-sort f tmp_order.txt --make-pgen multiallelics=+any --out tmp_jn2
mp_check tmp_small tmp_jn2
# '+snps' rejoins the A/C/G/T/N variants; the others (including 255-allele
# dosage records) are written as single-variant groups
$P --pfile tmp_sp4 --make-pgen multiallelics=+snps --out tmp_jn3
mp_check tmp_flip tmp_jn3
# tmp_joinsrc covers rescaled and missing overflows, inconsistent and
# conflicting hardcalls, implicit and hardcall-derived dosage-phase, etc.;
# compared against CompareJoin()'s reference implementation.
CJ="./tmp_multidosage_writer_test comparejoin"
$P --pfile tmp_joinsrc --validate --out tmp_joinsrc_v
$P --pfile tmp_joinsrc --make-pgen multiallelics=+any --out tmp_jn4
$P --pfile tmp_jn4 --validate --out tmp_jn4_v
$CJ tmp_joinsrc tmp_jn4
awk 'NR > 1 && (NR % 3 != 0) {print $1}' tmp_joinsrc.psam > tmp_keep_join.txt
awk '{a[NR] = $1} END {for (i = NR; i; --i) print a[i]}' tmp_keep_join.txt > tmp_order_join.txt
$P --pfile tmp_joinsrc --keep tmp_keep_join.txt --indiv-sort f tmp_order_join.txt --make-pgen multiallelics=+any --out tmp_jn5
$P --pfile tmp_jn5 --validate --out tmp_jn5_v
$CJ tmp_joinsrc tmp_jn5
$P --pfile tmp_joinsrc --make-pgen multiallelics=+any erase-phase --out tmp_jn6
$P --pfile tmp_jn6 --validate --out tmp_jn6_v
$CJ tmp_joinsrc tmp_jn6 strip-phase
# not supported yet with dosages
for f in "--hard-call-threshold 0.1" "--dosage-erase-threshold 0.1"; do
    if $P --pfile tmp_joinsrc $f --make-pgen multiallelics=+any --out tmp_jn_err; then
        echo "variant-join accepted $f"
        exit 1
    fi
done
