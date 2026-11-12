#!/bin/bash

# --adjust-file input field-name flags, none of which had test coverage:
# --adjust-chr-field, --adjust-pos-field, --adjust-id-field,
# --adjust-ref-field, --adjust-alt-field, --adjust-provref-field,
# --adjust-a1-field, --adjust-test-field, and --adjust-p-field.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Two tests per variant, so TEST must be honored.
body='1\t100\tv1\tA\tG\tN\tG\tADD\t0.01\n1\t100\tv1\tA\tG\tN\tG\tDOM\t0.5\n1\t200\tv2\tC\tT,A\tY\tT\tADD\t0.0001\n1\t200\tv2\tC\tT,A\tY\tT\tDOM\t0.2\n2\t300\tv3\tG\tC\tN\tC\tADD\t0.3\n2\t300\tv3\tG\tC\tN\tC\tDOM\t0.003\n'
printf "#CHROM\tPOS\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tTEST\tP\n$body" > tmp_std.txt
printf "CHROMOSOME\tCOORD\tMARKER\tR\tALTS\tPR\tEFF\tMODEL\tPVAL\n$body" > tmp_ren.txt
fields="--adjust-chr-field CHROMOSOME --adjust-pos-field COORD --adjust-id-field MARKER --adjust-ref-field R --adjust-alt-field ALTS --adjust-provref-field PR --adjust-a1-field EFF --adjust-test-field MODEL --adjust-p-field PVAL"

# Expected ADD results, computed by hand for p = 0.0001, 0.01, 0.3 (m = 3).
$plink2 --adjust-file tmp_std.txt test=ADD cols=chrom,pos,ref,alt,provref,a1,unadj,bonf,holm,sidakss,sidaksd,fdrbh,fdrby --out tmp_std > /dev/null
diff tmp_std.adjusted - <<'EOF'
#CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	UNADJ	BONF	HOLM	SIDAK_SS	SIDAK_SD	FDR_BH	FDR_BY
1	200	v2	C	T,A	Y	T	0.0001	0.0003	0.0003	0.00029997	0.00029997	0.0003	0.00055
1	100	v1	A	G	N	G	0.01	0.03	0.02	0.029701	0.0199	0.015	0.0275
2	300	v3	G	C	N	C	0.3	0.9	0.3	0.657	0.3	0.3	0.55
EOF

# Renamed columns: identical output with the field flags, an error without.
for cols in "cols=chrom,pos,ref,alt,provref,a1,unadj,gc,bonf,holm,sidakss,sidaksd,fdrbh,fdrby" "cols=chrom,pos,ref,alt1,provref,a1,unadj,fdrbh"; do
    for test in ADD DOM; do
        $plink2 --adjust-file tmp_std.txt test=$test $cols --out tmp_std_$test > /dev/null
        $plink2 --adjust-file tmp_ren.txt test=$test $cols $fields --out tmp_ren_$test > /dev/null
        cmp tmp_std_$test.adjusted tmp_ren_$test.adjusted
    done
done
if $plink2 --adjust-file tmp_ren.txt test=ADD --out tmp_nofields > /dev/null 2>&1; then
    exit 1
fi
# Without --adjust-test-field, MODEL isn't recognized, so test= is an error.
if $plink2 --adjust-file tmp_ren.txt test=ADD ${fields/--adjust-test-field MODEL/} --out tmp_notest > /dev/null 2>&1; then
    exit 1
fi
# Each field flag matters on its own: requesting a column whose field flag is
# missing is an error.
for pair in "ref:--adjust-ref-field R" "alt:--adjust-alt-field ALTS" "provref:--adjust-provref-field PR" "chrom:--adjust-chr-field CHROMOSOME" "pos:--adjust-pos-field COORD"; do
    col=${pair%%:*}
    flag=${pair#*:}
    if $plink2 --adjust-file tmp_ren.txt test=ADD cols=+$col ${fields/$flag/} --out tmp_missing_$col > /dev/null 2>&1; then
        exit 1
    fi
done

# A1 is silently omitted when it can't be found.
$plink2 --adjust-file tmp_ren.txt test=ADD ${fields/--adjust-a1-field EFF/} --out tmp_noa1 > /dev/null
if head -n 1 tmp_noa1.adjusted | grep -q "$(printf '\tA1\t')"; then
    exit 1
fi
head -n 1 tmp_ren_ADD.adjusted | grep -q "$(printf '\tA1\t')"

# With several names, earlier names take precedence; a missing name falls
# through to the next one.
printf "#CHROM\tID\tP\tPVAL\tA1\n1\tv1\t0.04\t0.001\tG\n1\tv2\t0.5\t0.2\tT\n" > tmp_two_p.txt
$plink2 --adjust-file tmp_two_p.txt cols=unadj --out tmp_default_p > /dev/null
test "$(tail -n +2 tmp_default_p.adjusted | cut -f 2 | tr '\n' ' ')" = "0.04 0.5 "
$plink2 --adjust-file tmp_two_p.txt cols=unadj --adjust-p-field PVAL P --out tmp_pval_first > /dev/null
test "$(tail -n +2 tmp_pval_first.adjusted | cut -f 2 | tr '\n' ' ')" = "0.001 0.2 "
$plink2 --adjust-file tmp_two_p.txt cols=unadj --adjust-p-field P PVAL --out tmp_p_first > /dev/null
cmp tmp_default_p.adjusted tmp_p_first.adjusted
$plink2 --adjust-file tmp_two_p.txt cols=unadj --adjust-p-field NOPE PVAL --out tmp_fallthrough > /dev/null
cmp tmp_pval_first.adjusted tmp_fallthrough.adjusted
# ID: SNP2 before ID.
printf "#CHROM\tID\tSNP2\tP\n1\tv1\trs1\t0.04\n1\tv2\trs2\t0.5\n" > tmp_two_id.txt
$plink2 --adjust-file tmp_two_id.txt cols=unadj --adjust-id-field SNP2 ID --out tmp_id > /dev/null
test "$(tail -n +2 tmp_id.adjusted | cut -f 1 | tr '\n' ' ')" = "rs1 rs2 "

# input-log10 with a renamed column.
printf "MARKER\tLOG10P\n" > tmp_log10.txt
printf "v1\t2\nv2\t4\nv3\t0.5\n" >> tmp_log10.txt
$plink2 --adjust-file tmp_log10.txt input-log10 cols=unadj,bonf --adjust-id-field MARKER --adjust-p-field LOG10P --out tmp_log10 > /dev/null
diff tmp_log10.adjusted - <<'EOF'
#ID	UNADJ	BONF
v2	0.0001	0.0003
v1	0.01	0.03
v3	0.316228	0.948683
EOF
