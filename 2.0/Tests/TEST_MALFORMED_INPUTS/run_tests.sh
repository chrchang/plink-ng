#!/bin/bash

# Malformed or unusual inputs that used to crash plink2 (segfault, abort, or
# a read past the end of an array) instead of either working or failing with
# an error message.  Each case must now exit normally.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Runs a command that must fail with an error, not a signal.
fails_cleanly() {
    local status=0
    "$@" > /dev/null 2>&1 || status=$?
    if [ $status -eq 0 ] || [ $status -ge 128 ]; then
        echo "expected a clean failure, got exit status $status: $*"
        exit 1
    fi
}

$plink2 --dummy 20 10 0.1 --seed 1 --make-pgen --out tmp_data > /dev/null
id0=$(awk '!/^#/ {print $3; exit}' tmp_data.pvar)

# 1. --update-map with a negative position removes the variant, but the
#    variant count was never written back, so every later step iterated one
#    variant too far.
printf '%s -1\n' $id0 > tmp_um.txt
$plink2 --pfile tmp_data --update-map tmp_um.txt --make-pgen --sort-vars --out tmp_um > /dev/null
test "$(grep -vc '^#' tmp_um.pvar)" -eq 9
$plink2 --pfile tmp_um --freq --out tmp_um > /dev/null

# 2. A .kin0 header without a KINSHIP column: the column scan walked past the
#    end of the header line.
printf '#FID1\tIID1\tFID2\tIID2\tNSNP\n' > tmp_nokin.kin0
fails_cleanly $plink2 --pfile tmp_data --king-cutoff-table tmp_nokin.kin0 0.1 --out tmp_bad
grep -q 'No kinship-coefficient column' tmp_bad.log
fails_cleanly $plink2 --pfile tmp_data --make-king-table --king-table-subset tmp_nokin.kin0 0.01 --out tmp_bad
grep -q 'No kinship-coefficient column' tmp_bad.log

# 3. A --pheno/--covar header with nothing after the sample IDs dereferenced
#    a null token pointer.
awk 'BEGIN {OFS = "\t"} NR == 1 {print "#FID", "IID", "SEX"; next} {print $1, $1, $2}' tmp_data.psam > tmp_fid.psam
awk 'NR == 1 {print "FID IID"; next} {print $1, $1}' tmp_data.psam > tmp_nocols.txt
$plink2 --pgen tmp_data.pgen --pvar tmp_data.pvar --psam tmp_fid.psam --pheno tmp_nocols.txt --freq --out tmp_nocols > /dev/null
$plink2 --pgen tmp_data.pgen --pvar tmp_data.pvar --psam tmp_fid.psam --covar tmp_nocols.txt --freq --out tmp_nocols > /dev/null

# 4. A malformed ##INFO header line: the error message's %s had no argument.
printf '##INFO=<Foo=AF>\n#CHROM\tPOS\tID\tREF\tALT\tINFO\n1\t10\tr1\tG\tA\tAF=0.1\n' > tmp_info.pvar
fails_cleanly $plink2 --pvar tmp_info.pvar --make-just-pvar --out tmp_bad
grep -q 'Header line 1 of tmp_info.pvar is malformed' tmp_bad.log

# 5. A NUL byte after the chromosome code was taken as a delimiter by the
#    main parser but as the end of the line by the INFO reload paths.
printf '##INFO=<ID=NS,Number=1,Type=Integer,Description="x">\n#CHROM\tPOS\tID\tREF\tALT\tINFO\n1\0\t10\tr1\tG\tA\tNS=3\n' > tmp_nul.pvar
fails_cleanly $plink2 --pvar tmp_nul.pvar --make-just-pvar --out tmp_bad
grep -q 'Invalid character after the chromosome code' tmp_bad.log

# 6. --adjust-file test= without a TEST column read past the end of a stack
#    array.
printf 'CHROM\tID\tP\n1\tr1\t0.1\n1\tr2\t0.3\n' > tmp_notest.txt
fails_cleanly $plink2 --adjust-file tmp_notest.txt test=ADD --out tmp_bad
grep -q 'tmp_notest.txt has no TEST' tmp_bad.log

# 7. --meta-analysis with an SE so small that 1/se^2 overflows produced NaN
#    and aborted in ZscoreToLnP(); such a row is now skipped like an invalid
#    one, so rsA is left in one file only.
printf 'CHR SNP BP A1 A2 OR SE\n1 rsA 100 A G 1.2 1e-300\n1 rsB 200 A G 1.5 0.2\n' > tmp_m1.txt
printf 'CHR SNP BP A1 A2 OR SE\n1 rsA 100 A G 1.1 0.1\n1 rsB 200 A G 1.4 0.2\n' > tmp_m2.txt
$plink2 --meta-analysis tmp_m1.txt tmp_m2.txt --out tmp_meta > /dev/null
if grep -qw rsA tmp_meta.meta; then
    echo "rsA should be dropped"
    exit 1
fi
grep -qw rsB tmp_meta.meta

# 8. A repeated --check-sex modifier: the error message's %s had no argument.
fails_cleanly $plink2 --pfile tmp_data --check-sex min-male-xf=0.5 min-male-xf=0.6 --out tmp_bad
grep -q 'Multiple --check-sex min-male-xf= modifiers' tmp_bad.log

# 9. A frequency-file line whose second allele belongs to no allele of the
#    variant: the allele matcher never counted its matches down, so it read
#    one allele past the variant (the next variant's REF) and accepted it.
$plink2 --dummy 60 10 0.1 acgt --seed 1 --make-pgen --out tmp_fr > /dev/null
$plink2 --pfile tmp_fr --freq --out tmp_fr > /dev/null
# Pick a variant whose REF/ALT pair does not include the next variant's REF,
# and replace its ALT with that REF.
read -r bad_id next_ref <<< "$(awk '!/^#/ {if (prev_id != "" && $4 != prev_ref && $4 != prev_alt) {print prev_id, $4; exit} prev_id = $3; prev_ref = $4; prev_alt = $5}' tmp_fr.pvar)"
test -n "$next_ref"
awk -v id=$bad_id -v a=$next_ref 'BEGIN {OFS = "\t"} $2 == id {$4 = a} {print}' tmp_fr.afreq > tmp_fr_bad.afreq
$plink2 --pfile tmp_fr --flip-scan --flip-scan-ref-freq tmp_fr_bad.afreq --out tmp_fs > /dev/null
grep -q '1 entry skipped' tmp_fs.log
