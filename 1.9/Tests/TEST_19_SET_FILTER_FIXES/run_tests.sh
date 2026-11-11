#!/bin/bash

# Regression tests for PLINK 1.9 set, annotation and filter fixes.  Each check
# failed before its fix.  plink2 is the reference where it has an equivalent
# command.
#
# Usage: ./run_tests.sh {plink 1.9 dir} {plink2 build dir} {up to 2 plink2 args}
# plink (1.9) is taken from {plink 1.9 dir}, or from the system PATH if that
# is empty or omitted.  {plink2 build dir} defaults to
# ../../../2.0/build_dynamic.

set -uo pipefail

# plink 1.9 from {plink 1.9 dir}, or from the PATH if that is empty or omitted.
plink19="${1:+$1/}plink"
plink() {
    command "$plink19" "$@"
}
plink2="${2:-../../../2.0/build_dynamic}/plink2 ${3:-} ${4:-}"

failed=0
check() {
    if [[ "$2" -eq 0 ]]; then
        echo "ok: $1"
    else
        echo "FAILED: $1"
        failed=1
    fi
}

# --- Sets and annotation ---

# 600 variants, 1 kb apart.
plink --dummy 4 600 --seed 1 --make-bed --out tmp_sets > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1000 * NR; print }' tmp_sets.bim > tmp_sets_pos.bim
mv tmp_sets_pos.bim tmp_sets.bim

# A set covering most variants is stored inverted.  save_set_bitfield()
# checked whether its last range ended exactly at the end of the 128-variant
# window ('==' instead of '>='), so a range ending past it got an overlapping
# extra range, and --write-set listed those variants twice.  Two such sets
# (variants 200-599, and everything except 250-259) must match plink2.
{
    echo SUFFIX; for i in $(seq 200 599); do echo snp$i; done; echo END
    echo GAP; for i in $(seq 0 249) $(seq 260 599); do echo snp$i; done; echo END
} > tmp_inv.set
plink --bfile tmp_sets --set tmp_inv.set --write-set --out tmp_inv19 > /dev/null
$plink2 --bfile tmp_sets --set tmp_inv.set --write-set --out tmp_inv2 > /dev/null
test -s tmp_inv2.set && cmp -s tmp_inv19.set tmp_inv2.set
check "--write-set of sets stored inverted matches plink2" $?

# --gene was applied when sets came from --make-set, but not from a --set
# file.  (plink2's --gene is a variant filter that keeps every set, so the
# reference is PLINK 1.9's own --set-names.)
{
    echo S1; for i in $(seq 300 2 398); do echo snp$i; done; echo END
    echo S2; for i in $(seq 0 9); do echo snp$i; done; echo END
    echo S3; for i in $(seq 590 599); do echo snp$i; done; echo END
} > tmp_gene.set
plink --bfile tmp_sets --set tmp_gene.set --gene S3 --write-set --out tmp_gene1 > /dev/null
plink --bfile tmp_sets --set tmp_gene.set --set-names S3 --write-set --out tmp_gene2 > /dev/null
test -s tmp_gene2.set && cmp -s tmp_gene1.set tmp_gene2.set
check "--gene with a --set file keeps only the named set" $?

# --annotate with --annotate-snp-field longer than 3 characters never matched
# the report's header, since header tokens were only compared up to 3
# characters.  The result must match the same report with a 'SNP' header.
printf ' CHR        SNP     BP        P\n   1      rs1    100    0.01\n   1      rs2    200     0.5\n   1      rs3    300   0.002\n' > tmp_rep.txt
sed 's/        SNP/ VARIANT_ID/' tmp_rep.txt > tmp_rep_long.txt
printf 'rs1 coding\nrs3 intron utr\n' > tmp_attr.txt
plink --annotate tmp_rep.txt attrib=tmp_attr.txt --out tmp_an1 > /dev/null
plink --annotate tmp_rep_long.txt attrib=tmp_attr.txt --annotate-snp-field VARIANT_ID --out tmp_an2 > /dev/null
test -s tmp_an1.annot && test "$(sed 's/ VARIANT_ID/        SNP/' tmp_an2.annot 2> /dev/null)" = "$(cat tmp_an1.annot)"
check "--annotate-snp-field with a long field name" $?

# --- Filters and options ---

# --oblig-missing counted duplicate (variant, block) entries twice, so the
# obligatory-missing counts went negative (wrapping around in .lmiss).
# Samples 0-9 form block c1, and their calls at v0-v9 are missing; 4 entries
# of the variant file are repeated.  --geno/--mind must match plink2, and
# --missing must match a run without the repeats.
awk 'BEGIN {
  s = 5
  for (j = 0; j < 30; ++j) print "1 v" j " 0 " (1000 * (j + 1)) > "tmp_om.map"
  for (i = 0; i < 20; ++i) {
    line = "f" i " i" i " 0 0 1 1"
    for (j = 0; j < 30; ++j) {
      s = (s * 69069 + 1) % 4294967296; u = s / 4294967296
      s = (s * 69069 + 1) % 4294967296; r = int(s / 4294967296 * 3)
      if ((i < 10 && j < 10) || u < 0.08) g = "0 0"
      else g = (r == 0)? "A A" : ((r == 1)? "A C" : "C C")
      line = line " " g
    }
    print line > "tmp_om.ped"
    if (i < 10) print "f" i, "i" i, "c1" > "tmp_om_samples.txt"
  }
  for (j = 0; j < 10; ++j) print "v" j, "c1" > "tmp_om_vars.txt"
}'
{ cat tmp_om_vars.txt; head -n 4 tmp_om_vars.txt; } > tmp_om_vars_dup.txt
plink --file tmp_om --make-bed --out tmp_om > /dev/null
om_ok=0
plink --bfile tmp_om --oblig-missing tmp_om_vars_dup.txt tmp_om_samples.txt --geno 0.1 --mind 0.2 --make-just-bim --out tmp_om_geno19 > /dev/null
$plink2 --bfile tmp_om --oblig-missing tmp_om_vars_dup.txt tmp_om_samples.txt --geno 0.1 --mind 0.2 --make-just-bim --out tmp_om_geno2 > /dev/null
test -s tmp_om_geno2.bim && cmp -s tmp_om_geno19.bim tmp_om_geno2.bim || om_ok=1
plink --bfile tmp_om --oblig-missing tmp_om_vars_dup.txt tmp_om_samples.txt --missing --out tmp_om_dup > /dev/null
plink --bfile tmp_om --oblig-missing tmp_om_vars.txt tmp_om_samples.txt --missing --out tmp_om_uniq > /dev/null
cmp -s tmp_om_dup.lmiss tmp_om_uniq.lmiss && cmp -s tmp_om_dup.imiss tmp_om_uniq.imiss || om_ok=1
check "--oblig-missing with repeated entries" $om_ok

# --bp-space compared last_pos + min_bp_space in int32, which overflowed for
# positions near 2^31, so nearby variants were kept.  Positions run from
# 1.90e9 to 2.12e9, 20 Mb apart; the variants kept must match plink2's.
# (Only CHR/ID/POS are compared, since PLINK 1.9 makes the minor allele A1
# while plink2 keeps the input allele order.)
plink --dummy 5 12 --seed 1 --make-bed --out tmp_bp > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 1900000000 + 20000000 * (NR - 1); print }' tmp_bp.bim > tmp_bp_pos.bim
mv tmp_bp_pos.bim tmp_bp.bim
plink --bfile tmp_bp --bp-space 150000000 --make-just-bim --out tmp_bp19 > /dev/null
$plink2 --bfile tmp_bp --bp-space 150000000 --make-just-bim --out tmp_bp2 > /dev/null
test -s tmp_bp2.bim && test "$(cut -f1,2,4 tmp_bp19.bim)" = "$(cut -f1,2,4 tmp_bp2.bim)"
check "--bp-space near the int32 position limit matches plink2" $?

# --cm-map <file> <chr> treated an '@' in the file name as the chromosome
# placeholder anyway, so it looked for map11.txt.  (plink2 rejects an '@'
# filename together with a chromosome, so the reference is the same map
# under a name without '@'.)
plink --dummy 4 10 --seed 1 --make-bed --out tmp_cm > /dev/null
awk 'BEGIN { OFS = "\t" } { $4 = 100000 * NR; print }' tmp_cm.bim > tmp_cm_pos.bim
mv tmp_cm_pos.bim tmp_cm.bim
printf 'position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n50000 1.0 0.05\n400000 2.0 0.4\n1200000 0.5 2.0\n' > 'tmp_map@1.txt'
cp 'tmp_map@1.txt' tmp_map_plain.txt
plink --bfile tmp_cm --cm-map 'tmp_map@1.txt' 1 --make-just-bim --out tmp_cm_at > /dev/null 2>&1
plink --bfile tmp_cm --cm-map tmp_map_plain.txt 1 --make-just-bim --out tmp_cm_plain > /dev/null
test -s tmp_cm_plain.bim && cmp -s tmp_cm_at.bim tmp_cm_plain.bim && \
    test "$(cut -f3 tmp_cm_plain.bim | head -n 1)" != "0"
check "--cm-map with a chromosome and '@' in the file name" $?

# With an '@' pattern, a missing per-chromosome map file is skipped with a
# warning, but it was also logged as an error ("Failed to open ...") while
# the run exited 0; and when no file existed at all, the run still exited 0
# with no cM positions set.  Now a partial match logs only the warning, and
# no match at all is an error.
awk 'BEGIN { OFS = "\t" } { $1 = (NR <= 5)? 1 : 2; print }' tmp_cm.bim > tmp_cm2.bim
cp tmp_cm.bed tmp_cm2.bed
cp tmp_cm.fam tmp_cm2.fam
mkdir -p tmp_cm_part tmp_cm_none
cp tmp_map_plain.txt tmp_cm_part/m1.txt
cmp_ok=0
plink --bfile tmp_cm2 --cm-map 'tmp_cm_part/m@.txt' --make-just-bim --out tmp_cm_part > /dev/null 2>&1 || cmp_ok=1
grep -q '^Warning: --cm-map failed to open tmp_cm_part/m2.txt' tmp_cm_part.log || cmp_ok=1
grep -q '^Error' tmp_cm_part.log && cmp_ok=1
test "$(cut -f3 tmp_cm_part.bim | head -n 5)" = "$(cut -f3 tmp_cm_plain.bim | head -n 5)" || cmp_ok=1
plink --bfile tmp_cm2 --cm-map 'tmp_cm_none/m@.txt' --make-just-bim --out tmp_cm_none > /dev/null 2>&1
rc=$?
test "$rc" -gt 0 && test "$rc" -lt 128 || cmp_ok=1
grep -q 'failed to open any file matching' tmp_cm_none.log || cmp_ok=1
check "--cm-map '@' pattern with missing map files" $cmp_ok
# plink2 already logged only the warning, but also exited 0 when no file
# existed.
cmp2_ok=0
$plink2 --bfile tmp_cm2 --cm-map 'tmp_cm_part/m@.txt' --make-just-bim --out tmp_cm2_part > /dev/null 2>&1 || cmp2_ok=1
grep -q '^Error' tmp_cm2_part.log && cmp2_ok=1
test "$(cut -f3 tmp_cm2_part.bim | head -n 5)" = "$(cut -f3 tmp_cm_part.bim | head -n 5)" || cmp2_ok=1
$plink2 --bfile tmp_cm2 --cm-map 'tmp_cm_none/m@.txt' --make-just-bim --out tmp_cm2_none > /dev/null 2>&1
rc=$?
test "$rc" -gt 0 && test "$rc" -lt 128 || cmp2_ok=1
grep -q 'failed to open any file matching' tmp_cm2_none.log || cmp2_ok=1
check "plink2 --cm-map '@' pattern with missing map files" $cmp2_ok

exit $failed
