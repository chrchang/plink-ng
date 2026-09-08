#!/bin/bash

# "--make-gene-masks": collapse each set into one dosage pseudo-variant.
#
# The masks are checked against the genotypes they came from, recomputed here
# from VCF exports of both filesets.  Going through VCF rather than --export A
# is deliberate: A counts A1, which is not the same allele on both sides, so
# comparing raw A columns would compare two different orientations.  A VCF DS
# field is always the alternate dosage.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

printf '400 null 0.05 0.95 1 1\n' > tmp_sim.txt
plink --simulate tmp_sim.txt --make-bed --out tmp_base --seed 33 > /dev/null

# Sets of 8 consecutive variants, plus a set naming a variant that does not
# exist and one naming nothing that passes the frequency ceiling.
awk '{ i = int((NR - 1) / 8);
       ids[i] = (ids[i] ? ids[i] "," : "") $2; chr[i] = $1;
       if (!(i in pos)) pos[i] = $4 }
     END { for (j = 0; j <= i; ++j) print "S" j, chr[j], pos[j], ids[j] }' tmp_base.bim > tmp_sets.txt
echo "GHOST 1 1 no_such_variant" >> tmp_sets.txt

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --make-gene-masks --set-list tmp_sets.txt --mask-max-af 0.5 --out tmp_max
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --make-gene-masks --build-mask sum --set-list tmp_sets.txt --mask-max-af 0.5 --out tmp_sum

# A set with no usable member produces no mask.
if grep -q '^1[[:space:]].*GHOST' tmp_max.pvar; then
    echo "A set with no usable variant still produced a mask."
    exit 1
fi

# The fileset has to be readable, and to have kept the sample IDs: a mask that
# cannot be joined back to the phenotypes is useless.
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_max --export vcf vcf-dosage=DS --out tmp_maxv
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_sum --export vcf vcf-dosage=DS --out tmp_sumv
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --export vcf --out tmp_basev
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --freq --out tmp_freq

if ! head -2 tmp_max.psam | grep -q 'FID'; then
    echo "Mask .psam lost the FID column."
    exit 1
fi

# Recompute both modes from the source genotypes.
check_mode() {
    local maskvcf=$1
    local mode=$2
    awk -v mode="$mode" -v lbl="$mode" '
        # alternate allele frequency per variant
        FILENAME == ARGV[1] && !/^#/ { af[$2] = $5; next }
        # set definitions
        FILENAME == ARGV[2] { members[$1] = $4; next }
        # source genotypes
        FILENAME == ARGV[3] {
            if ($0 ~ /^##/) next
            if ($0 ~ /^#CHROM/) { for (i = 10; i <= NF; ++i) bcol[i - 9] = $i; bn = NF - 9; next }
            n = split($9, f, ":")
            for (i = 10; i <= NF; ++i) {
                g = $i; split(g, gf, ":"); gt = gf[1]
                v = (gt ~ /\./) ? 0 : (substr(gt, 1, 1) + substr(gt, 3, 1))
                geno[$3 SUBSEP (i - 9)] = v
            }
            next
        }
        # masks
        {
            if ($0 ~ /^##/) next
            if ($0 ~ /^#CHROM/) { for (i = 10; i <= NF; ++i) mcol[i - 9] = $i; mn = NF - 9; next }
            split($9, f, ":"); dsidx = 0
            for (i = 1; i <= length(f); ++i) if (f[i] == "DS") dsidx = i
            if (!dsidx) { print "no DS field"; exit 1 }
            nm = split(members[$3], mem, ",")
            for (i = 10; i <= NF; ++i) {
                split($i, gf, ":")
                got = (gf[dsidx] == ".") ? 0 : gf[dsidx] + 0
                want = 0
                for (k = 1; k <= nm; ++k) {
                    if (!(mem[k] in af)) continue
                    if (af[mem[k]] + 0 > 0.5) continue
                    val = geno[mem[k] SUBSEP (i - 9)]
                    if (mode == "max") { if (val > want) want = val }
                    else { want += val }
                }
                if (mode != "max" && want > 2) want = 2
                d = got - want; if (d < 0) d = -d
                if (d > 1e-3) { print lbl " mismatch on " $3 " sample " i - 9 ": got " got ", want " want; exit 1 }
                ++cells
            }
        }
        END { if (cells < 1000) { print "too few cells compared: " cells; exit 1 }
              printf "%s: %d cells checked\n", lbl, cells }
    ' tmp_freq.afreq tmp_sets.txt tmp_basev.vcf "$maskvcf"
}
check_mode tmp_maxv.vcf max
check_mode tmp_sumv.vcf sum

# The whole point of writing a fileset: --glm applies to it unchanged.
{ echo "#FID IID QT"; awk '{print $1, $2, (NR % 71) / 71.0}' tmp_base.fam; } > tmp_qt.pheno
$BUILD/plink2 $EXTRA1 $EXTRA2 --pfile tmp_max --pheno tmp_qt.pheno --pheno-name QT --glm allow-no-covars --out tmp_g
if [[ $(grep -cv '^#' tmp_g.QT.glm.linear) != $(grep -cv '^#' tmp_max.pvar) ]]; then
    echo "--glm did not test every mask."
    exit 1
fi

# And --acat-file combines the per-mask results.
awk '!/^#/ { print "ALL 1 1 " $3 }' tmp_max.pvar | head -1 > /dev/null
awk 'BEGIN { s = "" } !/^#/ { s = (s ? s "," : "") $3 } END { print "ALLMASKS 1 1 " s }' tmp_max.pvar > tmp_maskset.txt
$BUILD/plink2 $EXTRA1 $EXTRA2 --acat-file tmp_g.QT.glm.linear test=ADD --set-list tmp_maskset.txt --out tmp_ac
if [[ $(grep -cv '^#' tmp_ac.acat) != 1 ]]; then
    echo "--acat-file did not combine the mask results."
    exit 1
fi

# Error cases.
fails() {
    if "$@" > /dev/null 2>&1; then
        echo "expected failure: $*"
        exit 1
    fi
}
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --make-gene-masks --out tmp_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --build-mask sum --out tmp_bad
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --make-gene-masks --set-list tmp_sets.txt --build-mask nonsense --out tmp_bad
# Nothing passes a ceiling of 1e-9 in this dataset.
fails $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_base --make-gene-masks --set-list tmp_sets.txt --mask-max-af 0.000000001 --out tmp_bad

echo "--make-gene-masks tests passed."
