#!/bin/bash

# --make-king-table HET1_HOM2 / HET2_HOM1 columns.
#
# HET1_HOM2 counts variants where the IID1 sample is heterozygous and the
# IID2 sample homozygous; HET2_HOM1 the reverse.  The full-table path used to
# report the two exchanged, while --king-table-subset reported them the right
# way around.  KINSHIP and the IBS column are symmetric in the two and were
# not affected.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# Three samples.  v1-v4 are ordinary variants; v5 and v6 are singletons (one
# sample differs from everyone else), which the full-table path counts
# separately.
cat > tmp_data.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
1	100	v1	A	G	.	.	.	GT	0/1	0/0	0/1
1	200	v2	A	G	.	.	.	GT	0/1	1/1	1/1
1	300	v3	A	G	.	.	.	GT	0/0	0/1	1/1
1	400	v4	A	G	.	.	.	GT	0/1	0/1	0/0
1	500	v5	A	G	.	.	.	GT	0/0	0/0	0/1
1	600	v6	A	G	.	.	.	GT	0/1	0/0	0/0
EOF
$plink2 --vcf tmp_data.vcf --make-pgen --out tmp_data > /dev/null

# Expected (IID1 het and IID2 hom, IID2 het and IID1 hom) for each pair:
#   A,B: A het at v1 v2 v6 with B hom, B het at v3 with A hom.
#   A,C: A het at v2 v4 v6 with C hom, C het at v5 with A hom.
#   B,C: B het at v3 v4 with C hom, C het at v1 v5 with B hom.
awk 'BEGIN { OFS = "\t"; print "A", "B", 3, 1; print "A", "C", 3, 1; print "B", "C", 2, 2 }' > tmp_expected.txt

# $1: .kin0 file.
check_kin0() {
    awk 'BEGIN { FS = OFS = "\t" }
         NR == FNR { want[$1 "," $2] = $3 OFS $4; want[$2 "," $1] = $4 OFS $3; next }
         FNR == 1 { sub(/^#/, ""); for (i = 1; i <= NF; ++i) { col[$i] = i }; next }
         {
           got = $col["HET1_HOM2"] OFS $col["HET2_HOM1"]
           k = $col["IID1"] "," $col["IID2"]
           if (got != want[k]) { print FILENAME ": " k " has " got ", expected " want[k]; failed = 1; exit 1 }
           ++n
         }
         END { if (failed) { exit 1 }; if (n != 3) { print "expected 3 pairs, got " n + 0; exit 1 } }' tmp_expected.txt $1
}

$plink2 --pfile tmp_data --make-king-table counts cols=id,hethet,ibs0,ibs1,kinship --out tmp_full > /dev/null
check_kin0 tmp_full.kin0

# The subset path, fed the pairs in both orders.
printf '#IID1\tIID2\nA\tB\nA\tC\nB\tC\n' > tmp_pairs.kin0
$plink2 --pfile tmp_data --make-king-table counts cols=id,hethet,ibs0,ibs1,kinship --king-table-subset tmp_pairs.kin0 --out tmp_subset > /dev/null
check_kin0 tmp_subset.kin0
