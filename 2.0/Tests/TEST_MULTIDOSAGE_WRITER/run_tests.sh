#!/bin/bash

# Round-trip unit test for multiallelic dosages: the writer
# (SpgwAppendMultiallelicGenovecDphase16()) and the readers (PgrGetMDp(),
# PgrGet1D(), PgrGetMDCounts(), PgrGetRaw(), etc.), and the validator
# (PgrValidate()).  Builds its own driver from the
# pgenlib sources with plain flags, rather than reusing the build directory's
# objects (whose SIMD flags vary).  Set CXX to override the compiler.

set -exo pipefail

INC=../../include
${CXX:-c++} -std=c++17 -O1 -g -I$INC -o tmp_multidosage_writer_test multidosage_writer_test.cc $INC/pgenlib_read.cc $INC/pgenlib_write.cc $INC/pgenlib_misc.cc $INC/plink2_base.cc $INC/plink2_bits.cc $INC/plink2_simd.cc -lz
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

# 'keep-dosage': cleared heterozygous hardcalls keep their dosages, and those
# without one get a dosage equal to the hardcall (one copy of each allele;
# rare-alt entries for ALT2+).  Dosage-phase values of affected samples are
# erased.
$P --pfile tmp_mt --set-mixed-mt-missing keep-dosage --make-pgen --out tmp_kd1
mp_check tmp_mt tmp_kd1 hets-missing keep-dosage
# chrY without sex information: all hets
cp tmp_small.pgen tmp_y.pgen
cp tmp_small.psam tmp_y.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "Y"; print}' tmp_small.pvar > tmp_y.pvar
$P --pfile tmp_y --set-invalid-haploid-missing keep-dosage --make-pgen --out tmp_kd2
mp_check tmp_y tmp_kd2 hets-missing keep-dosage
# chrX: only males are affected (and females keep their phase)
cp tmp_small.pgen tmp_x.pgen
awk 'BEGIN {OFS = "\t"} NR == 1 {print "#IID", "SEX"; next} {print $1, (NR % 2)? 1 : 2}' tmp_small.psam > tmp_x.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "X"; print}' tmp_small.pvar > tmp_x.pvar
awk 'NR > 1 && $2 == 1 {print $1}' tmp_x.psam > tmp_x_males.txt
$P --pfile tmp_x --set-invalid-haploid-missing keep-dosage --make-pgen --out tmp_kd3
mp_check tmp_x tmp_kd3 hets-missing keep-dosage hm-samples=tmp_x_males.txt
$P --pfile tmp_x --set-invalid-haploid-missing --make-pgen --out tmp_kd4
mp_check tmp_x tmp_kd4 hets-missing hm-samples=tmp_x_males.txt
# chrY with females (whose calls are all set to missing, dosages included):
# just --validate
cp tmp_x.psam tmp_ys.psam
cp tmp_small.pgen tmp_ys.pgen
cp tmp_y.pvar tmp_ys.pvar
$P --pfile tmp_ys --set-invalid-haploid-missing keep-dosage --make-pgen --out tmp_kd5
$P --pfile tmp_kd5 --validate --out tmp_kd5_v
$P --pfile tmp_ys --set-invalid-haploid-missing --make-pgen --out tmp_kd6
$P --pfile tmp_kd6 --validate --out tmp_kd6_v
# hardcall-only multiallelic input: new dosages only
cp tmp_flipnd.pgen tmp_ndmt.pgen
cp tmp_flipnd.psam tmp_ndmt.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "MT"; print}' tmp_flipnd.pvar > tmp_ndmt.pvar
$P --pfile tmp_ndmt --set-mixed-mt-missing keep-dosage --make-pgen --out tmp_kd7
mp_check tmp_ndmt tmp_kd7 hets-missing keep-dosage
# multithreaded writer with sample subsetting
$P --pfile tmp_mt --threads 2 --keep tmp_keep.txt --set-mixed-mt-missing keep-dosage --make-pgen --out tmp_kd8
mp_check tmp_mt tmp_kd8 hets-missing keep-dosage

# --export vcf/bcf on multiallelic dosages: DS has one entry per ALT allele,
# and HDS lists the left haplotype's ALT dosages followed by the right's.
# Checked against the source by CheckExport() (BCF output too), on autosomes,
# chrX (males haploid), chrY with females, and MT.
CE="./tmp_multidosage_writer_test checkexport"
for f in tmp_small tmp_x tmp_ys tmp_mt; do
    for m in DS HDS DS-force HDS-force; do
        $P --pfile $f --export vcf bcf vcf-dosage=$m --out tmp_ex_${f}_$m
        $CE $f tmp_ex_${f}_$m.vcf $m
        $CE $f tmp_ex_${f}_$m.bcf $m
    done
done
$P --pfile tmp_small --export vcf bcf vcf-dosage=DS-only --out tmp_ex_ds_only
$CE tmp_small tmp_ex_ds_only.vcf DS-only
$CE tmp_small tmp_ex_ds_only.bcf DS-only
# sample subset (the checker matches samples by ID), and fixed-width input
$P --pfile tmp_x --keep tmp_keep.txt --export vcf bcf vcf-dosage=HDS --out tmp_ex_keep
$CE tmp_x tmp_ex_keep.vcf HDS
$CE tmp_x tmp_ex_keep.bcf HDS
$P --pfile tmp_fixed --export vcf bcf vcf-dosage=HDS-force --out tmp_ex_fixed
$CE tmp_fixed tmp_ex_fixed.vcf HDS-force
$CE tmp_fixed tmp_ex_fixed.bcf HDS-force
# GP isn't supported for multiallelic dosage records
$P --pfile tmp_small --export vcf vcf-dosage=GP --out tmp_ex_gp
grep -q "GP field not exported for 180 multiallelic variants with dosages" tmp_ex_gp.log
# 'A-transpose' with ALT2 counted (PgrGet1D()): must match the second DS
# entry of the DS-force VCF (3 vs. up to 5 decimal places)
awk 'BEGIN {OFS = "\t"} !/^#/ && (split($5, a, ",") >= 2) {print $3, a[2]}' tmp_small.pvar > tmp_ex_alt2.txt
$P --pfile tmp_small --export A-transpose --export-allele tmp_ex_alt2.txt --out tmp_ex_alt2
awk 'BEGIN {FS = "\t"}
FNR == 1 {++file_idx}
file_idx == 1 && !/^#/ && (split($5, a, ",") >= 2) {
  for (i = 10; i <= NF; ++i) {
    split($i, f, ":");
    if (f[2] == ".") {
      ds[$3, i - 9] = "NA";
    } else {
      split(f[2], v, ",");
      ds[$3, i - 9] = v[2];
    }
  }
  next;
}
file_idx == 2 && FNR > 1 && (($2, 1) in ds) {
  for (i = 7; i <= NF; ++i) {
    e = ds[$2, i - 6];
    ++ct;
    if ((e == "NA") != ($i == "NA")) {
      ++bad;
    } else if ((e != "NA") && ((e - $i > 0.0006) || ($i - e > 0.0006))) {
      ++bad;
    }
  }
}
END {
  print ct " A-transpose values, " (bad + 0) " mismatches";
  exit (bad || (ct != 180 * 300));
}' tmp_ex_tmp_small_DS-force.vcf tmp_ex_alt2.traw

# --glm (linear, logistic, multinomial) on multiallelic dosages: each
# non-omitted allele's column holds its dosages.  Checked against plink2's
# biallelic code paths by glm_reference.sh (see its comments), on autosomes,
# chrX under both --xchr-model 1 and 2, and MT (the last three without the
# 255-allele variants, to save time).
GC="awk -f glm_compare.awk"
glm_check() {
    ./glm_reference.sh "$P" $1 $2 $3 "$4"
    $GC ${2}_ref.txt ${2}_multi.Q.glm.linear ${2}_multi.Q2.glm.linear ${2}_multi.B.glm.logistic.hybrid ${2}_multi.M.glm.multinomial.hybrid
    $GC ${2}_ref.txt ${2}_multi1.Q.glm.linear
    $GC ${2}_ref.txt ${2}_multis.B.glm.logistic.hybrid
    $GC ${2}_permref.txt ${2}_dflt.Q.glm.linear ${2}_dflt.Q2.glm.linear ${2}_dflt.B.glm.logistic.hybrid ${2}_dflt.M.glm.multinomial.hybrid
    # alleles with zero dosage everywhere must be reported as constant (with
    # at most 5 alleles; variant-level errors like SAMPLE_CT<=PREDICTOR_CT
    # take precedence for the bigger ones)
    awk 'BEGIN {FS = "\t"} FNR == 1 {for (i = 1; i <= NF; ++i) c[$i] = i; next}
      ($c["TEST"] == "ADD") && ($c["A1_FREQ"] == 0) && (split($c["ALT"], alts, ",") <= 4) {++ct; if ($c["ERRCODE"] != "CONST_ALLELE") ++bad}
      END {print ct " zero-frequency alleles, " (bad + 0) " not CONST_ALLELE"; exit (bad || (!ct))}' ${2}_multi.Q.glm.linear ${2}_multi.B.glm.logistic.hybrid ${2}_multi.M.glm.multinomial.hybrid
}
glm_check tmp_small tmp_gl_small 6
glm_check tmp_x tmp_gl_x2 3 "--max-alleles 20"
glm_check tmp_x tmp_gl_x1 3 "--max-alleles 20 --xchr-model 1"
glm_check tmp_mt tmp_gl_mt 3 "--max-alleles 20"
# --condition on a multiallelic variant (previously hardcalls only): must
# match the same allele dosages given as covariates.
awk 'BEGIN {OFS = "\t"} /^##/ {next} /^#CHROM/ {for (i = 10; i <= NF; ++i) id[i] = $i; ns = NF; next}
$3 == "v73" {
  n = split($5, al, ",");
  hdr = "#IID";
  for (j = 1; j <= n; ++j) hdr = hdr "\tD" j;
  print hdr;
  for (i = 10; i <= ns; ++i) {
    split($i, f, ":");
    row = id[i];
    if (f[2] == ".") {
      for (j = 1; j <= n; ++j) row = row "\tNA";
    } else {
      split(f[2], d, ",");
      for (j = 1; j <= n; ++j) row = row "\t" d[j];
    }
    print row;
  }
}' tmp_gl_small_ds.vcf > tmp_gl_cond_cov.txt
GCOND="--glm allow-no-covars omit-ref hide-covar multinomial cols=+a1freq"
$P --pfile tmp_small --pheno tmp_gl_small_ph.txt --condition v73 multiallelic $GCOND --out tmp_gl_cond
$P --pfile tmp_small --pheno tmp_gl_small_ph.txt --covar tmp_gl_cond_cov.txt $GCOND --out tmp_gl_condref
rm -f tmp_gl_condref.txt
for ph in Q.glm.linear Q2.glm.linear B.glm.logistic.hybrid M.glm.multinomial.hybrid; do
    awk -v ph=$ph 'BEGIN {OFS = "\t"} !/^#/ && ($3 != "v73") {print ph, $3, $0}' tmp_gl_condref.$ph >> tmp_gl_condref.txt
done
$GC tmp_gl_condref.txt tmp_gl_cond.Q.glm.linear tmp_gl_cond.Q2.glm.linear tmp_gl_cond.B.glm.logistic.hybrid tmp_gl_cond.M.glm.multinomial.hybrid

# --make-rel/--make-grm-list/--pca on multiallelic dosages: each allele's row
# holds its centered and scaled dosage.  CheckGrm() computes the relationship
# matrix independently (all four 'cov'/'meanimpute' combinations), and
# CheckPca() checks the PCs against the --make-rel output.
$P --pfile tmp_small --freq --out tmp_rel_freq
for m in "" "meanimpute" "cov" "cov meanimpute"; do
    $P --pfile tmp_small --read-freq tmp_rel_freq.afreq --make-grm-list $m --out tmp_rel_grm
    ./tmp_multidosage_writer_test checkgrm tmp_small tmp_rel_freq.afreq tmp_rel_grm.grm $m
done
# GRM variant loading is parallel (decode threads) except with --threads 1;
# 3 threads give uneven decode segments.  Each decode thread tracks variants
# with missing calls separately (merged afterward), which the default mode's
# per-pair observation counts depend on.
for t in 1 3; do
    for m in "" "cov"; do
        $P --pfile tmp_small --read-freq tmp_rel_freq.afreq --threads $t --make-grm-list $m --out tmp_rel_grm
        ./tmp_multidosage_writer_test checkgrm tmp_small tmp_rel_freq.afreq tmp_rel_grm.grm $m
    done
done
$P --pfile tmp_small --read-freq tmp_rel_freq.afreq --make-rel square --out tmp_rel_sq
$P --pfile tmp_small --read-freq tmp_rel_freq.afreq --make-rel square meanimpute --out tmp_rel_sqmi
# (allele-wts exercises the other variant-loading passes)
$P --pfile tmp_small --read-freq tmp_rel_freq.afreq --pca 4 allele-wts --out tmp_rel_pca
./tmp_multidosage_writer_test checkpca tmp_rel_sq.rel tmp_rel_pca.eigenvec tmp_rel_pca.eigenval
$P --pfile tmp_small --read-freq tmp_rel_freq.afreq --pca 4 approx allele-wts --out tmp_rel_pcaa
./tmp_multidosage_writer_test checkpca tmp_rel_sqmi.rel tmp_rel_pcaa.eigenvec tmp_rel_pcaa.eigenval approx
# --pca variant loading is parallel (decode threads), except with --threads 1.
# With this few variants there's a single compute thread either way, so the
# results must be byte-identical.  3 threads give decode segments that split
# multiallelic variants at different points; the variant subset (every third
# variant excluded) makes batch boundaries land on excluded variants, and its
# last batch (84 rows) has a partial second segment under 3 threads.
awk '!/^#/ && (NR % 3) {print $3}' tmp_small.pvar > tmp_rel_subset.txt
for subset in "" "--extract tmp_rel_subset.txt"; do
    for mode in "approx" ""; do
        for t in 1 3 4; do
            $P --pfile tmp_small $subset --read-freq tmp_rel_freq.afreq --threads $t --seed 1 --pca 4 $mode allele-wts --out tmp_rel_pcat$t
        done
        for t in 3 4; do
            for ext in eigenval eigenvec eigenvec.allele; do
                cmp tmp_rel_pcat1.$ext tmp_rel_pcat$t.$ext
            done
        done
    done
done

# --score with rare ALT alleles counted (PgrGet1D()): ALT2 for every
# multiallelic variant, plus ALT3 where present (two rows for the same
# variant), under the additive/dominant/recessive models, with and without
# mean imputation, and on MT (haploid: dosages halved).  Expected sums are
# computed from the DS-force VCF by score_compare.awk.  (The 255-allele
# variants are left out to keep the VCF small.)
$P --pfile tmp_small --max-alleles 20 --export vcf vcf-dosage=DS-force --out tmp_sc_ds
$P --pfile tmp_small --max-alleles 20 --freq --out tmp_sc_freq
awk 'BEGIN {OFS = "\t"} !/^#/ {n = split($5, a, ","); if (n >= 2) print $3, a[2], ((NR * 37) % 11 - 5) / 7; if (n >= 3) print $3, a[3], ((NR * 53) % 13 - 6) / 9}' tmp_small.pvar > tmp_sc_file.txt
for mod in "" "dominant" "recessive"; do
    for imp in "no-mean-imputation" ""; do
        $P --pfile tmp_small --max-alleles 20 --read-freq tmp_sc_freq.afreq --score tmp_sc_file.txt $mod $imp cols=+scoresums --out tmp_sc
        awk -v mod="$mod" -v imp="$imp" -v haploid=0 -f score_compare.awk tmp_sc_freq.afreq tmp_sc_file.txt tmp_sc_ds.vcf tmp_sc.sscore
    done
done
for imp in "no-mean-imputation" ""; do
    $P --pfile tmp_mt --max-alleles 20 --read-freq tmp_sc_freq.afreq --score tmp_sc_file.txt $imp cols=+scoresums --out tmp_sc_mt
    awk -v mod="" -v imp="$imp" -v haploid=1 -f score_compare.awk tmp_sc_freq.afreq tmp_sc_file.txt tmp_sc_ds.vcf tmp_sc_mt.sscore
done
# PCA projection: projecting the same samples with the .eigenvec.allele
# weights must give back their PCs, up to one constant factor per PC.
# ('meanimpute' PCs, since a projection treats missing calls the same way;
# the default PCs use per-pair nonmissing counts.)
$P --pfile tmp_small --read-freq tmp_rel_freq.afreq --pca 4 meanimpute allele-wts --out tmp_sc_pca
$P --pfile tmp_small --read-freq tmp_rel_freq.afreq --score tmp_sc_pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize cols=+scoresums --score-col-nums 6-9 --out tmp_sc_proj
awk 'function absf(x) { return (x < 0)? -x : x }
     FNR == 1 { next }
     NR == FNR { for (k = 1; k <= 4; ++k) e[FNR, k] = $(k + 1); next }
     { for (k = 1; k <= 4; ++k) s[FNR, k] = $(NF - 4 + k); n = FNR - 1 }
     END {
       for (k = 1; k <= 4; ++k) {
         se = ss = see = sss = ses = 0;
         for (i = 2; i <= n + 1; ++i) {
           x = e[i, k]; y = s[i, k];
           se += x; ss += y; see += x * x; sss += y * y; ses += x * y;
         }
         c = (ses - se * ss / n) / sqrt((see - se * se / n) * (sss - ss * ss / n));
         print "PC" k " projection correlation " c;
         if (absf(c) < 0.99999) bad = 1;
       }
       exit bad;
     }' tmp_sc_pca.eigenvec tmp_sc_proj.sscore
# --variant-score uses total ALT dosages (PgrGetD()'s collapse), with
# missing dosages imputed from the total ALT frequency.
awk 'BEGIN {OFS = "\t"; print "#IID", "W1", "W2"} NR > 1 {print $1, ((NR * 29) % 17 - 8) / 5, ((NR * 41) % 7) / 3}' tmp_small.psam > tmp_sc_vw.txt
for sp in "" "single-prec"; do
    $P --pfile tmp_small --max-alleles 20 --read-freq tmp_sc_freq.afreq --variant-score tmp_sc_vw.txt $sp cols=+altfreq,+nmiss --out tmp_sc_vs
    awk -f vscore_compare.awk tmp_sc_vw.txt tmp_sc_ds.vcf tmp_sc_vs.vscore
done

# --r2-[un]phased and --ld on multiallelic dosages: the major allele's dosage
# and dosage-phase vs. the rest, compared against the corresponding split
# variants (see ld_reference.sh), on autosomes, chrX, and MT.
mat_compare() {
    # all-pairs r^2 matrices; nan must match nan
    awk 'function absf(x) { return (x < 0)? -x : x }
         FNR == NR { for (i = 1; i <= NF; ++i) v[FNR, i] = $i; rows = FNR; next }
         { for (i = 1; i <= NF; ++i) { ++n; x = v[FNR, i]; y = $i;
             if (((x == "nan") != (y == "nan")) || ((x != "nan") && (absf(x - y) > 1e-6 + 1e-4 * absf(x)))) ++bad } }
         END { print n " r^2 entries, " (bad + 0) " mismatches"; exit (bad || (n != rows * rows) || (!n)) }' "$1" "$2"
}
for f in tmp_small tmp_x tmp_mt; do
    ./ld_reference.sh "$P" $f tmp_ld_$f
    for m in unphased phased; do
        $P --pfile tmp_ld_${f}_perm --r2-$m square --out tmp_ld_${f}_m
        $P --pfile tmp_ld_${f}_split --extract tmp_ld_${f}_split_ext.txt --r2-$m square --out tmp_ld_${f}_s
        mat_compare tmp_ld_${f}_m.$m.vcor2 tmp_ld_${f}_s.$m.vcor2
    done
done
# Unpermuted (93 multiallelic variants have REF as the major allele): same
# matrices.
for m in unphased phased; do
    $P --pfile tmp_small --max-alleles 5 --r2-$m square --out tmp_ld_orig
    mat_compare tmp_ld_orig.$m.vcor2 tmp_ld_tmp_small_m.$m.vcor2
done
# --ld on pairs of multiallelic variants: r^2, |D'|, sample count, and the
# haplotype frequencies (as a set, since a split variant's major allele can be
# "the rest" when frequencies are close to 0.5).
ld_summary() {
    sed -n '/^--ld /,$p' $1 | awk '/valid samples/ {s = $1} /r\^2 =/ {r = $3 " " $6}
      /^ *MAJOR  [0-9]/ || /^ *MINOR  [0-9]/ {f[++n] = $2; f[++n] = $3}
      END {for (i = 1; i <= n; ++i) for (j = i + 1; j <= n; ++j) if (f[j] < f[i]) {t = f[i]; f[i] = f[j]; f[j] = t}
           printf("%s %s", s, r); for (i = 1; i <= n; ++i) printf(" %s", f[i]); printf("\n")}'
}
for pair in "v72 v73" "v84 v85" "v96 v97" "v73 v97" "v102 v103" "v1 v2"; do
    set -- $pair
    id1=$(awk -v v=$1 'FNR == NR {if (!/^#/) {++k; if ($3 == v) idx = k}; next} FNR == idx' tmp_ld_tmp_small_perm.pvar tmp_ld_tmp_small_split_ext.txt)
    id2=$(awk -v v=$2 'FNR == NR {if (!/^#/) {++k; if ($3 == v) idx = k}; next} FNR == idx' tmp_ld_tmp_small_perm.pvar tmp_ld_tmp_small_split_ext.txt)
    $P --pfile tmp_ld_tmp_small_perm --ld $1 $2 --out tmp_ld_pair_m
    $P --pfile tmp_ld_tmp_small_split --ld $id1 $id2 --out tmp_ld_pair_s
    test "$(ld_summary tmp_ld_pair_m.log)" == "$(ld_summary tmp_ld_pair_s.log)"
done
# --indep-pairwise uses hardcalls (dosages only affect allele frequencies):
# with shared frequencies, dosages must make no difference.
$P --pfile tmp_small --freq --out tmp_ld_ipf
$P --pfile tmp_small --make-pgen erase-dosage --out tmp_ld_nd
$P --pfile tmp_small --read-freq tmp_ld_ipf.afreq --indep-pairwise 50 5 0.01 --out tmp_ld_ip
$P --pfile tmp_ld_nd --read-freq tmp_ld_ipf.afreq --indep-pairwise 50 5 0.01 --out tmp_ld_ipnd
cmp tmp_ld_ip.prune.in tmp_ld_ipnd.prune.in
cmp tmp_ld_ip.prune.out tmp_ld_ipnd.prune.out

# --rm-dup on duplicate-ID multiallelic dosage records (see WriteDups()):
# identical pairs are duplicates, while the pairs differing only in one
# rare-alt dosage or dosage-phase value are mismatches.
./tmp_multidosage_writer_test writedups tmp_small tmp_rd
$P --pfile tmp_rd --validate --out tmp_rd_v
$P --pfile tmp_rd --rm-dup retain-mismatch list --make-just-pvar --out tmp_rd_rm
diff -q tmp_rd_rm.rmdup.mismatch tmp_rd_mismatch.txt
test "$(grep -v '^#' tmp_rd_rm.pvar | cut -f 3 | tr '\n' ' ')" == "v72 v73 v74 v75 v78 v78 v90 v90 "
$P --pfile tmp_rd --rm-dup exclude-mismatch --make-just-pvar --out tmp_rd_ex
test "$(grep -v '^#' tmp_rd_ex.pvar | cut -f 3 | tr '\n' ' ')" == "v72 v73 v74 v75 "

# --sample-diff 'dosage[=tol]' on multiallelic dosage records: a pair differs
# when some allele's dosage differs by more than the tolerance; the report's
# geno column has one DS value per ALT allele.  Checked against the DS-force
# VCF by sdiff_compare.awk (main report and counts-only), on autosomes and
# chrX (female pairs diploid, male pairs haploid).
SD="--min-alleles 3 --max-alleles 20"
$P --pfile tmp_small $SD --export vcf vcf-dosage=DS-force --out tmp_sd_ds
for tol in 0 0.1; do
    for im in 0 1; do
        mods="dosage=$tol"
        if [ $im = 1 ]; then mods="$mods include-missing"; fi
        $P --pfile tmp_small $SD --sample-diff $mods counts-cols=+nobs,+ibs0,+ibs1 base=s0 s1 s2 s3 --out tmp_sd
        awk -v tol=$tol -v include_missing=$im -f sdiff_compare.awk tmp_sd_ds.vcf tmp_sd.sdiff.summary tmp_sd.s0.sdiff
        $P --pfile tmp_small $SD --sample-diff $mods counts-only counts-cols=+nobs,+ibs0,+ibs1 base=s0 s1 s2 s3 --out tmp_sd_c
        awk -v tol=$tol -v include_missing=$im -f sdiff_compare.awk tmp_sd_ds.vcf tmp_sd_c.sdiff.summary
    done
done
cp tmp_x.pgen tmp_sd_xauto.pgen
cp tmp_x.psam tmp_sd_xauto.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "1"; print}' tmp_x.pvar > tmp_sd_xauto.pvar
$P --pfile tmp_sd_xauto $SD --export vcf vcf-dosage=DS-force --out tmp_sd_xds
for tol in 0 0.1; do
    $P --pfile tmp_x $SD --sample-diff dosage=$tol include-missing counts-cols=+nobs,+ibs0,+ibs1 base=s0 s2 s4 s6 --out tmp_sd_xf
    awk -v tol=$tol -v include_missing=1 -f sdiff_compare.awk tmp_sd_xds.vcf tmp_sd_xf.sdiff.summary tmp_sd_xf.s0.sdiff
    $P --pfile tmp_x $SD --sample-diff dosage=$tol include-missing counts-cols=+nobs,+ibs0,+ibs1 base=s1 s3 s5 s7 --out tmp_sd_xm
    awk -v tol=$tol -v include_missing=1 -v haploid=1 -f sdiff_compare.awk tmp_sd_xds.vcf tmp_sd_xm.sdiff.summary tmp_sd_xm.s1.sdiff
    $P --pfile tmp_x $SD --sample-diff dosage=$tol include-missing counts-only counts-cols=+nobs,+ibs0,+ibs1 base=s1 s3 s5 s7 --out tmp_sd_xmc
    awk -v tol=$tol -v include_missing=1 -v haploid=1 -f sdiff_compare.awk tmp_sd_xds.vcf tmp_sd_xmc.sdiff.summary
done

# --pgen-diff 'dosage[=tol]' with more than two merged alleles: a sample
# differs when some allele's dosage differs by more than the tolerance.
# Checked by pdiff_compare.awk against the second fileset being a
# dosage-erased copy, an ALT-permuted copy (no differences), both, and a
# trim-alts + dosage-erased copy (mostly biallelic); and on chrX (males
# haploid).
PD="--min-alleles 3 --max-alleles 5"
cp tmp_small.pgen tmp_pd_a.pgen
cp tmp_x.psam tmp_pd_a.psam
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "1"; print}' tmp_small.pvar > tmp_pd_a.pvar
awk '!/^#/ && (split($5, a, ",") >= 2) && (split($5, a, ",") <= 19) {print $3 "\t" a[2]}' tmp_pd_a.pvar > tmp_pd_alt2.txt
$P --pfile tmp_pd_a $PD --make-pgen erase-dosage --out tmp_pd_b_erase
$P --pfile tmp_pd_a $PD --alt1-allele force tmp_pd_alt2.txt --make-pgen --out tmp_pd_b_perm
$P --pfile tmp_pd_b_perm --make-pgen erase-dosage --out tmp_pd_b_permerase
$P --pfile tmp_pd_a $PD --make-pgen trim-alts erase-dosage --out tmp_pd_b_trim
$P --pfile tmp_pd_a $PD --export vcf vcf-dosage=DS-force --out tmp_pd_va
# (all tolerance/include-missing combinations for the first one, to save
# time)
for b in b_erase b_perm b_permerase b_trim; do
    $P --pfile tmp_pd_$b --export vcf vcf-dosage=DS-force --out tmp_pd_v$b
    for tol in 0 0.1; do
        for im in 0 1; do
            if [ $b != b_erase ] && { [ $tol = 0 ] || [ $im = 0 ]; }; then continue; fi
            mods="dosage=$tol"
            if [ $im = 1 ]; then mods="$mods include-missing"; fi
            $P --pfile tmp_pd_a $PD --pgen-diff tmp_pd_$b $mods --out tmp_pd
            awk -v tol=$tol -v include_missing=$im -f pdiff_compare.awk tmp_pd_va.vcf tmp_pd_v$b.vcf tmp_pd.pdiff
        done
    done
done
test "$(grep -v '^#' tmp_pd.pdiff | wc -l)" -gt 0
# trim-alts copy that keeps its dosages: no differences.  (Variant order is
# reversed first, so that its biallelic records with dosages follow
# multiallelic ones with rare-alt dosages.)
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$2 = 1000 - $2; print}' tmp_pd_a.pvar > tmp_pd_rev0.pvar
cp tmp_pd_a.pgen tmp_pd_rev0.pgen
cp tmp_pd_a.psam tmp_pd_rev0.psam
$P --pfile tmp_pd_rev0 $PD --sort-vars --make-pgen --out tmp_pd_rev
$P --pfile tmp_pd_rev --make-pgen trim-alts --out tmp_pd_b_trimonly
$P --pfile tmp_pd_rev --pgen-diff tmp_pd_b_trimonly dosage include-missing --out tmp_pd_trimonly
grep -q "108 variants compared, 0 differences" tmp_pd_trimonly.log
test "$(grep -v '^#' tmp_pd_b_trimonly.pvar | awk -F '\t' 'split($5, a, ",") == 1' | wc -l)" -gt 0
# ALT-permuted copy with up to 20 alleles: still no differences
awk '!/^#/ && (split($5, a, ",") >= 2) && (split($5, a, ",") <= 19) {print $3 "\t" a[2]}' tmp_pd_a.pvar > tmp_pd_alt2.txt
$P --pfile tmp_pd_a --min-alleles 3 --max-alleles 20 --alt1-allele force tmp_pd_alt2.txt --make-pgen --out tmp_pd_b_perm20
$P --pfile tmp_pd_a --min-alleles 3 --max-alleles 20 --pgen-diff tmp_pd_b_perm20 dosage include-missing --out tmp_pd_perm20
test "$(grep -v '^#' tmp_pd_perm20.pdiff | wc -l)" -eq 0
grep -q "144 variants compared, 0 differences" tmp_pd_perm20.log
$P --pfile tmp_x $PD --make-pgen erase-dosage --out tmp_pd_x_erase
awk 'NR > 1 && $2 == 1 {print $1}' tmp_x.psam > tmp_pd_x_males.txt
for tol in 0 0.1; do
    $P --pfile tmp_x $PD --pgen-diff tmp_pd_x_erase dosage=$tol include-missing --out tmp_pd_x
    awk -v tol=$tol -v include_missing=1 -v x_males=tmp_pd_x_males.txt -f pdiff_compare.awk tmp_pd_va.vcf tmp_pd_vb_erase.vcf tmp_pd_x.pdiff
done

# --pmerge/--pmerge-list with multiallelic dosages.  comparemerge recomputes
# each output genotype from the sources (matching variants by ID and alleles
# by string, so allele permutations are covered) under the given
# --merge-mode, and compares hardcalls, phase, per-allele dosages, and
# per-allele dosage-phase values in the merged allele space.
# Sources: tmp_pm_a = tmp_small; tmp_pm_bperm = ALT-permuted copy;
# tmp_pm_bnd = ALT-permuted copy without dosages (dosage conflicts);
# tmp_pm_ahm = copy with heterozygous hardcalls (and their dosages) set
# missing; tmp_pm_alo/tmp_pm_bhi/tmp_pm_bndhi = overlapping sample subsets;
# tmp_pm_brev/tmp_pm_bndrev = sample-order-reversed copies; tmp_pm_bswap =
# copy with ALT2 and ALT3 relabeled (where there are at least 3 ALTs), so that
# rare-alt dosages move between alleles while REF/ALT1 hardcalls and alt-sums
# stay the same.
cp tmp_small.pgen tmp_pm_a.pgen
cp tmp_small.pvar tmp_pm_a.pvar
cp tmp_small.psam tmp_pm_a.psam
awk '!/^#/ && (split($5, a, ",") >= 2) {print $3 "\t" a[2]}' tmp_pm_a.pvar > tmp_pm_alt2.txt
$P --pfile tmp_pm_a --alt1-allele force tmp_pm_alt2.txt --make-pgen --out tmp_pm_bperm
$P --pfile tmp_pm_bperm --make-pgen erase-dosage --out tmp_pm_bnd
# (--set-mixed-mt-missing applies to chrM only)
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "MT"; print}' tmp_pm_a.pvar > tmp_pm_amt.pvar
cp tmp_pm_a.pgen tmp_pm_amt.pgen
cp tmp_pm_a.psam tmp_pm_amt.psam
$P --pfile tmp_pm_amt --set-mixed-mt-missing --make-pgen --out tmp_pm_ahm0
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$1 = "1"; print}' tmp_pm_ahm0.pvar > tmp_pm_ahm.pvar
cp tmp_pm_ahm0.pgen tmp_pm_ahm.pgen
cp tmp_pm_ahm0.psam tmp_pm_ahm.psam
awk 'NR > 1 && NR <= 201 {print $1}' tmp_pm_a.psam > tmp_pm_lo.txt
awk 'NR > 101 {print $1}' tmp_pm_a.psam > tmp_pm_hi.txt
awk 'NR > 1 {ids[NR] = $1} END {for (i = NR; i > 1; --i) print ids[i]}' tmp_pm_a.psam > tmp_pm_rev.txt
$P --pfile tmp_pm_a --keep tmp_pm_lo.txt --make-pgen --out tmp_pm_alo
$P --pfile tmp_pm_bperm --keep tmp_pm_hi.txt --make-pgen --out tmp_pm_bhi
$P --pfile tmp_pm_bnd --keep tmp_pm_hi.txt --make-pgen --out tmp_pm_bndhi
$P --pfile tmp_pm_bhi --indiv-sort f tmp_pm_rev.txt --make-pgen --out tmp_pm_brev
$P --pfile tmp_pm_bndhi --indiv-sort f tmp_pm_rev.txt --make-pgen --out tmp_pm_bndrev
test "$(awk 'NR == 2 {print $1}' tmp_pm_brev.psam)" = "s299"
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {n = split($5, a, ","); if (n >= 3) {t = a[2]; a[2] = a[3]; a[3] = t; s = a[1]; for (i = 2; i <= n; ++i) s = s "," a[i]; $5 = s} print}' tmp_pm_a.pvar > tmp_pm_bswap.pvar
cp tmp_pm_a.pgen tmp_pm_bswap.pgen
cp tmp_pm_a.psam tmp_pm_bswap.psam
# Each scenario: sources, in merge order.  Two sources use --pmerge, three use
# --pmerge-list (and are repeated with 2-fileset passes, so that a temporary
# fileset gets merged).
pm_scenarios="tmp_pm_a,tmp_pm_bperm tmp_pm_a,tmp_pm_bnd tmp_pm_bnd,tmp_pm_a tmp_pm_ahm,tmp_pm_bperm tmp_pm_alo,tmp_pm_bndhi tmp_pm_alo,tmp_pm_bndrev tmp_pm_brev,tmp_pm_alo tmp_pm_ahm,tmp_pm_brev tmp_pm_a,tmp_pm_bswap tmp_pm_ahm,tmp_pm_bndhi,tmp_pm_a"
pm_rare_dphase_seen=0
pm_conflicts_seen=0
for mode in first nm-first nm-match; do
    for scenario in $pm_scenarios; do
        srcs=$(echo $scenario | tr ',' ' ')
        set -- $srcs
        if [ $# = 2 ]; then
            pass_opts="single"
        else
            pass_opts="single multi"
        fi
        for pass_opt in $pass_opts; do
            if [ $# = 2 ]; then
                $P --pfile $1 --pmerge $2 --merge-mode $mode --make-pgen --out tmp_pm_out
            else
                rm -f tmp_pm_list.txt
                for src in $srcs; do echo $src >> tmp_pm_list.txt; done
                extra=""
                if [ $pass_opt = multi ]; then extra="--debug pmerge-pass-size=2"; fi
                $P --pmerge-list tmp_pm_list.txt pfile --merge-mode $mode $extra --make-pgen --out tmp_pm_out
            fi
            $P --pfile tmp_pm_out --validate --out tmp_pm_out_v
            ./tmp_multidosage_writer_test comparemerge $mode tmp_pm_out $srcs | tee tmp_pm_cmp.txt
            if ! grep -q " 0 with rare-alt dosage-phase" tmp_pm_cmp.txt; then pm_rare_dphase_seen=1; fi
            if ! grep -q " 0 conflicts" tmp_pm_cmp.txt; then pm_conflicts_seen=1; fi
        done
    done
done
test $pm_rare_dphase_seen = 1
test $pm_conflicts_seen = 1
# Record pairs that differ only in a rare-alt dosage (or only in a rare-alt
# dosage-phase value) of one sample, from writedups above: the two copies are
# split into separate filesets, and nm-match must report exactly one conflict
# for each pair.
awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {$3 = $3 "_" (++ct[$3]); print}' tmp_rd.pvar > tmp_pm_rdtag.pvar
cp tmp_rd.pgen tmp_pm_rdtag.pgen
cp tmp_rd.psam tmp_pm_rdtag.psam
for copy in 1 2; do
    awk -v c=$copy '!/^#/ && ($3 ~ ("_" c "$")) {print $3}' tmp_pm_rdtag.pvar > tmp_pm_rd_ids.txt
    $P --pfile tmp_pm_rdtag --extract tmp_pm_rd_ids.txt --make-pgen --out tmp_pm_rd${copy}t
    awk 'BEGIN {OFS = "\t"} /^#/ {print; next} {sub(/_[12]$/, "", $3); print}' tmp_pm_rd${copy}t.pvar > tmp_pm_rd${copy}.pvar
    cp tmp_pm_rd${copy}t.pgen tmp_pm_rd${copy}.pgen
    cp tmp_pm_rd${copy}t.psam tmp_pm_rd${copy}.psam
done
$P --pfile tmp_pm_rd1 --pmerge tmp_pm_rd2 --make-pgen --out tmp_pm_rdm
$P --pfile tmp_pm_rdm --validate --out tmp_pm_rdm_v
./tmp_multidosage_writer_test comparemerge nm-match tmp_pm_rdm tmp_pm_rd1 tmp_pm_rd2 | tee tmp_pm_cmp.txt
grep -q " 2 conflicts" tmp_pm_cmp.txt
# Variant-disjoint sources (concatenation path).
$P --pfile tmp_pm_a --from v0 --to v107 --make-pgen --out tmp_pm_v1
$P --pfile tmp_pm_bperm --from v108 --to v215 --make-pgen --out tmp_pm_v2
$P --pfile tmp_pm_v1 --pmerge tmp_pm_v2 --make-pgen --out tmp_pm_cat
$P --pfile tmp_pm_cat --validate --out tmp_pm_cat_v
./tmp_multidosage_writer_test comparemerge nm-match tmp_pm_cat tmp_pm_v1 tmp_pm_v2
# --merge-ignore-dosage/--merge-ignore-phase: same as merging copies without
# dosages/phase.
$P --pfile tmp_pm_a --make-pgen erase-dosage --out tmp_pm_a_nd
$P --pfile tmp_pm_a --make-pgen erase-phase --out tmp_pm_a_np
$P --pfile tmp_pm_bperm --make-pgen erase-phase --out tmp_pm_b_np
for ig in dosage phase; do
    if [ $ig = dosage ]; then b_erased=tmp_pm_bnd; else b_erased=tmp_pm_b_np; fi
    $P --pfile tmp_pm_a --pmerge tmp_pm_bperm --merge-ignore-$ig --make-pgen --out tmp_pm_ig
    $P --pfile tmp_pm_a_n${ig:0:1} --pmerge $b_erased --make-pgen --out tmp_pm_ig_ref
    $P --pfile tmp_pm_ig --validate --out tmp_pm_ig_v
    ./tmp_multidosage_writer_test comparemerge nm-match tmp_pm_ig tmp_pm_ig_ref
done
