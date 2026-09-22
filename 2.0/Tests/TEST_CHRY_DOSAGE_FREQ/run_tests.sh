#!/bin/bash

# chrY allele counts with dosages, when unknown-sex samples are present.
#
# By default, unknown-sex samples count toward chrY allele frequencies but not
# toward chrY missingness, so males and unknown-sex samples are counted on a
# path of their own.  That path used to add every sample's dosage, females'
# included, and paired each dosage with the first sample's hardcall instead of
# its own.  --y-nosex-missing-stats counts the same samples through another
# path, so the two must agree, and both must match the counts worked out by
# hand below.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Three males, one unknown-sex sample, two females.  Haploid dosages far from
# 0 and 1 (0.4, 0.2, 0.8) import with a missing hardcall.
#   y1: nonfemale ALT dosages 1, 0.4, 0, 0: 1.4 over 4 observations.
#   y2: nonfemale ALT dosages 0.2, 1, missing, 0.8: 2 over 3 observations.
#   y3: only a female has a fractional dosage; nonfemale 1, 0, 0, 0: 1 over 4.
cat > tmp_data.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=Y,length=60000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Alternate allele dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M1	M2	M3	N1	F1	F2
Y	5000000	y1	A	G	.	.	.	GT:DS	1:1	0:0.4	0:0	0:0	1:0.9	0:0.3
Y	5000100	y2	C	T	.	.	.	GT:DS	0:0.2	1:1	.:.	1:0.8	0:0	1:1
Y	5000200	y3	G	A	.	.	.	GT:DS	1:1	0:0	0:0	0:0	0:0.7	1:1
EOF
printf '#IID\tSEX\nM1\t1\nM2\t1\nM3\t1\nN1\tNA\nF1\t2\nF2\t2\n' > tmp_data.psam
$plink2 --vcf tmp_data.vcf dosage=DS --psam tmp_data.psam --make-pgen --out tmp_data > /dev/null

$plink2 --pfile tmp_data --freq counts cols=+nobs --out tmp_default > /dev/null
$plink2 --pfile tmp_data --y-nosex-missing-stats --freq counts cols=+nobs --out tmp_nosex_stats > /dev/null

awk 'BEGIN { FS = "\t"; want["y1"] = "1.4\t4"; want["y2"] = "2\t3"; want["y3"] = "1\t4" }
     FNR == 1 { next }
     {
       got = $5 "\t" $6
       if (got != want[$2]) { print FILENAME ": " $2 " has ALT_CTS/OBS_CT " got ", expected " want[$2]; exit 1 }
       ++n
     }
     END { if (n != 6) { print "expected 6 rows, got " n + 0; exit 1 } }' tmp_default.acount tmp_nosex_stats.acount

# The same with the first sample removed, so the counting runs on a sample
# subset and the first included sample is not sample 0.  y1 then has ALT
# dosages 0.4, 0, 0 (0.4 over 3), y2 has 1, missing, 0.8 (1.8 over 2), and y3
# has 0, 0, 0 (0 over 3).
printf 'M1\n' > tmp_remove.txt
$plink2 --pfile tmp_data --remove tmp_remove.txt --freq counts cols=+nobs --out tmp_subset > /dev/null
awk 'BEGIN { FS = "\t"; want["y1"] = "0.4\t3"; want["y2"] = "1.8\t2"; want["y3"] = "0\t3" }
     FNR == 1 { next }
     {
       got = $5 "\t" $6
       if (got != want[$2]) { print FILENAME ": " $2 " has ALT_CTS/OBS_CT " got ", expected " want[$2]; exit 1 }
       ++n
     }
     END { if (n != 3) { print "expected 3 rows, got " n + 0; exit 1 } }' tmp_subset.acount
