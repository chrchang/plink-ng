#!/bin/bash

# --test-missing: differential missingness between cases and controls, tested
# against PLINK 1.9 on autosomes, chrX and chrY.
#
# The chrX and chrY parts matter because PLINK 1.x counts a male heterozygous
# call on a haploid chromosome as missing, and restricts chrY to males.

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.05 --simulate-ncases 400 --simulate-ncontrols 400 --out tmp_data

# Relabel the three 100-variant blocks as chr1, chrX and chrY, and give half
# the samples male sex so that the haploid special cases are exercised.
awk 'BEGIN{OFS="\t"} {if (NR <= 100) $1 = 1; else if (NR <= 200) $1 = 23; else $1 = 24; print}' tmp_data.bim > tmp_relabeled.bim
mv tmp_relabeled.bim tmp_data.bim
awk 'BEGIN{OFS=" "} {$5 = (NR % 2) + 1; print}' tmp_data.fam > tmp_relabeled.fam
mv tmp_relabeled.fam tmp_data.fam

compare() {
    # $1: PLINK 1.9 .missing report, $2: plink2 .missing report.
    awk '
        function abs(x) { return (x < 0)? -x : x }
        function close_enough(a, b) { return abs(a - b) <= 1e-6 + 1e-3 * abs(a) }
        FNR == NR {
            if (FNR > 1) {
                fa[$2] = $3; fu[$2] = $4; p[$2] = $5; ++n1;
            }
            next
        }
        /^#/ { next }
        {
            ++n2;
            if (!($3 in fa)) { print "extra row in plink2 output: " $3; failed = 1; exit 1 }
            if (!close_enough(fa[$3], $4) || !close_enough(fu[$3], $5) || !close_enough(p[$3], $6)) {
                print "mismatch on " $3 ": " fa[$3] "/" fu[$3] "/" p[$3] " vs " $4 "/" $5 "/" $6;
                failed = 1; exit 1
            }
            delete fa[$3];
        }
        END {
            if (failed) { exit 1 }
            if (n1 != n2) { print "row count mismatch: " n1 " vs " n2; exit 1 }
            if (n1 == 0) { print "no rows compared"; exit 1 }
            print n1 " rows matched"
        }
    ' "$1" "$2"
}

# 1. Fisher's exact test.
plink --bfile tmp_data --test-missing --out plink19
$1/plink2 $2 $3 --bfile tmp_data --test-missing --out plink2
compare plink19.missing plink2.missing

# 2. mid-p adjustment.
plink --bfile tmp_data --test-missing midp --out plink19_midp
$1/plink2 $2 $3 --bfile tmp_data --test-missing midp --out plink2_midp
compare plink19_midp.missing plink2_midp.missing

# 3. Same results after a round trip through .pgen, which stores the male
#    heterozygous chrX/chrY calls the same way.
$1/plink2 $2 $3 --bfile tmp_data --make-pgen --out tmp_pgen
$1/plink2 $2 $3 --pfile tmp_pgen --test-missing --out plink2_pgen
diff -q plink2.missing plink2_pgen.missing

# 4. 'zs' output is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --test-missing zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.missing.zst > plink2_zs.missing
diff -q plink2.missing plink2_zs.missing

# 5. Column sets.  The counts and denominators have to be consistent with the
#    rates the default report gives, and the column set has to control both
#    the header and the row contents.
$1/plink2 $2 $3 --bfile tmp_data --test-missing cols=chrom,pos,ref,alt,nmissa,nobsa,fmissa,nmissu,nobsu,fmissu,p --out plink2_cols
head -n 1 plink2_cols.missing | grep -q '^#CHROM	POS	ID	REF	ALT	MISSING_CT_A	OBS_CT_A	F_MISS_A	MISSING_CT_U	OBS_CT_U	F_MISS_U	P$'
awk -F '\t' '
    function abs(x) { return (x < 0)? -x : x }
    NR > 1 {
        if (abs($6 / $7 - $8) > 1e-9 || abs($9 / $10 - $11) > 1e-9) {
            print "count/rate mismatch on " $3; exit 1
        }
        ++n
    }
    END { if (n == 0) { print "no rows checked"; exit 1 } }
' plink2_cols.missing
# The rates themselves have to match the default report.
diff -q <(cut -f 1,2,3,8,11,12 plink2_cols.missing) plink2.missing

$1/plink2 $2 $3 --bfile tmp_data --test-missing cols=p --out plink2_p
head -n 1 plink2_p.missing | grep -q '^#ID	P$'
# The leading '#' moves onto the ID column here, so compare the bodies.
diff -q <(tail -n +2 plink2.missing | cut -f 3,6) <(tail -n +2 plink2_p.missing)

# 6. Dosage missingness.  The fixture above is hardcall-only, so build one
#    where the two notions differ: 'dosage-freq=1' gives every call a dosage,
#    and a strict --hard-call-threshold then erases most of the hardcalls
#    without touching the dosages.
$1/plink2 $2 $3 --dummy 400 200 0.05 dosage-freq=1 --seed 7 --hard-call-threshold 0.1 --make-pgen --out tmp_dose
awk 'NR == 1 {print "#IID\tCC"; next} {print $1 "\t" (NR % 2? 2 : 1)}' tmp_dose.psam > tmp_dose_pheno.txt
$1/plink2 $2 $3 --pfile tmp_dose --pheno tmp_dose_pheno.txt --test-missing cols=chrom,pos,nmissa,nobsa,nmissu,nobsu --out plink2_hc
$1/plink2 $2 $3 --pfile tmp_dose --pheno tmp_dose_pheno.txt --test-missing dosage cols=chrom,pos,nmissa,nobsa,nmissu,nobsu --out plink2_dos
# Dosage missingness is a subset of hardcall missingness, and here a strict
# proper subset.
awk -F '\t' 'FNR == NR {if (FNR > 1) {a[$3] = $4; u[$3] = $6}; next}
    FNR > 1 {
        if ($4 > a[$3] || $6 > u[$3]) { print "dosage missingness exceeds hardcall on " $3; exit 1 }
        if ($4 < a[$3] || $6 < u[$3]) { ++strictly_fewer }
        ++n
    }
    END {
        if (n == 0) { print "no rows checked"; exit 1 }
        if (strictly_fewer == 0) { print "fixture has no dosage-only calls"; exit 1 }
    }' plink2_hc.missing plink2_dos.missing

# Both must agree with the totals --missing reports for the same file.
$1/plink2 $2 $3 --pfile tmp_dose --missing variant-only vcols=chrom,pos,nmiss,nmissdosage --out plink2_vmiss
check_totals() {
    # $1: --test-missing report, $2: the --missing column to match.
    awk -F '\t' -v col="$2" '
        FNR == NR {
            if (FNR == 1) {for (i = 1; i <= NF; ++i) {if ($i == col) {idx = i}}; next}
            total[$3] = $idx; next
        }
        FNR > 1 {
            if ($4 + $6 != total[$3]) {
                print "total mismatch on " $3 ": " $4 " + " $6 " != " total[$3]; exit 1
            }
            ++n
        }
        END { if (n == 0) { print "no rows checked"; exit 1 } }
    ' plink2_vmiss.vmiss "$1"
}
check_totals plink2_hc.missing MISSING_CT
check_totals plink2_dos.missing MISSING_DOSAGE_CT

# 7. An ALTx/ALTy call on a haploid chromosome is a heterozygous call, and so
#    a het haploid, even though the raw genotype vector makes it look
#    homozygous.  --missing counts it, so this has to as well.
cat > tmp_hh.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=Y>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	m1	m2	m3	m4
Y	1	yv_alt12	A	G,T	.	.	.	GT	1/2	0/0	0/0	0/0
Y	2	yv_refalt	A	G	.	.	.	GT	0/1	0/0	0/0	0/0
EOF
$1/plink2 $2 $3 --vcf tmp_hh.vcf --make-pgen --out tmp_hh
printf '#IID\tSEX\tCC\nm1\t1\t2\nm2\t1\t2\nm3\t1\t1\nm4\t1\t1\n' > tmp_hh.psam
$1/plink2 $2 $3 --pfile tmp_hh --test-missing cols=chrom,nmissa,nobsa,nmissu,nobsu --out plink2_hh
# Both variants have exactly one het haploid, in the same sample, so the two
# rows have to be identical apart from the ID.
test "$(grep -c '^Y' plink2_hh.missing)" -eq 2
diff -q <(awk 'NR > 1 {print $1, $3, $4, $5, $6}' plink2_hh.missing | sort -u) <(awk 'NR > 1 {print $1, $3, $4, $5, $6}' plink2_hh.missing | sort -u | head -n 1)
# ...and --missing agrees that both are het haploids.
$1/plink2 $2 $3 --pfile tmp_hh --missing variant-only vcols=chrom,hethap --out plink2_hh_m
awk 'NR > 1 && $3 != 1 {print "hethap count on " $2 " is " $3; exit 1}' plink2_hh_m.vmiss
