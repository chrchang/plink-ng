#!/bin/bash

# --adjust-file on input without a PROVISIONAL_REF? column.
#
# The default column set includes 'maybeprovref', which only reports the
# column when provisional REF alleles are present.  --adjust-file used to
# treat it as mandatory, so with default columns it failed on any input
# without a PROVISIONAL_REF? column, e.g. PLINK 1.x association reports.
# Forcing the column with 'provref' must still require it.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# A PLINK 1.x-style report: whitespace-delimited, CHR/SNP/BP/A1 headers.
cat > tmp_legacy.assoc <<'EOF'
 CHR  SNP   BP  A1  NMISS  BETA  SE  R2  T  P
   1   v1  100   A    50  0.1  0.1  0.01  1  0.3
   1   v2  200   C    50  0.2  0.1  0.02  2  0.004
   1   v3  300   G    50  0.3  0.1  0.03  3  0.02
   2   v4  400   T    50  0.4  0.1  0.04  4  0.7
EOF

# Default columns: must succeed, and match an explicit column set that leaves
# maybeprovref out.
$plink2 --adjust-file tmp_legacy.assoc --out tmp_default > /dev/null
$plink2 --adjust-file tmp_legacy.assoc cols=chrom,a1,unadj,gc,bonf,holm,sidakss,sidaksd,fdrbh,fdrby --out tmp_explicit > /dev/null
cmp tmp_default.adjusted tmp_explicit.adjusted
test "$(tail -n +2 tmp_default.adjusted | wc -l)" -eq 4

# Forced provref: still an error.
if $plink2 --adjust-file tmp_legacy.assoc cols=+provref --out tmp_forced > /dev/null; then
    echo "cols=+provref succeeded without a PROVISIONAL_REF? column"
    exit 1
fi
grep -q 'No PROVISIONAL_REF? column' tmp_forced.log

# With REF and PROVISIONAL_REF? columns and one provisional REF allele,
# maybeprovref reports the column.
printf '#CHROM\tPOS\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tP\n1\t100\tv1\tG\tA\tN\tA\t0.3\n1\t200\tv2\tA\tC\tY\tC\t0.004\n' > tmp_provref.txt
$plink2 --adjust-file tmp_provref.txt cols=+ref --out tmp_provref > /dev/null
head -n 1 tmp_provref.adjusted | grep -q 'PROVISIONAL_REF?'
awk -F '\t' '$2 == "v2" { if ($4 != "Y") { print "v2 PROVISIONAL_REF? is " $4; exit 1 } }' tmp_provref.adjusted
