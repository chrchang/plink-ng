#!/bin/bash

# --clump-verbose, --clump-best and --clump-annotate.
#
# The verbose rows are compared against the block PLINK 1.9 embeds in its
# .clumped file.  1.9's verbose mode also reports every index candidate as its
# own block rather than only the clumps it formed -- its table grows from 14
# rows to 52 on this fixture -- so the comparison covers the clumps plink2
# reports and ignores 1.9's extra blocks.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

awk -f make_vcf.awk > tmp_h.vcf
$BUILD/plink2 $EXTRA1 $EXTRA2 --vcf tmp_h.vcf --double-id --make-bed --out tmp_data
awk 'BEGIN{OFS=" "} {print $1, $2, 0, 0, (NR % 2) + 1, (NR % 2) + 1}' tmp_data.fam > tmp_ped.fam
mv tmp_ped.fam tmp_data.fam

plink --bfile tmp_data --assoc --out tmp_a1 > /dev/null
head -n 200 tmp_data.fam | cut -d' ' -f1,2 > tmp_half.txt
plink --bfile tmp_data --keep tmp_half.txt --assoc --out tmp_a2 > /dev/null

# Two files so that the F column is present on both sides, and thresholds
# loose enough for clumps to form: the phenotype is unrelated to the
# genotypes.
FILES="tmp_a1.assoc,tmp_a2.assoc"
# PLINK 1.9 rejects --clump-best with two files unless --clump-index-first is
# also given, so that is on throughout.
P1="--clump-p1 0.5 --clump-p2 0.8 --clump-index-first"

# 1. Verbose rows against PLINK 1.9.
plink --bfile tmp_data --clump $FILES $P1 --clump-verbose --clump-best --out plink19
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump $FILES $P1 --clump-verbose --clump-best --out plink2
awk -f compare_verbose.awk plink19.clumped plink2.clumps.verbose

# 2. The verbose file has to agree with the .clumps table it came from: the
#    index variant plus one row per SP2 member, which is the below-p2 subset
#    of the clump rather than all TOTAL of its members.
awk '
    FNR == NR {
        if (!/^#/ && NF) {
            sp2_ct = ($NF == ".")? 0 : (split($NF, a, ",") + 0);
            expected[$3] = sp2_ct + 1;
        }
        next
    }
    /^#/ { next }
    { ++seen[$1] }
    END {
        for (id in seen) {
            if (!(id in expected)) { print "verbose clump " id " is not in the .clumps table"; exit 1 }
            if (seen[id] != expected[id]) {
                print id ": " seen[id] " verbose rows, expected " expected[id]; exit 1
            }
            ++n;
        }
        if (!n) { print "no clumps checked"; exit 1 }
        print n " clumps agree with the .clumps table"
    }
' plink2.clumps plink2.clumps.verbose

# 3. --clump-best names one member per clump, and it has to be a member with
#    the highest r^2 in that clump.  The index variant's own row is the one
#    whose ALLELES field names a single allele rather than an in-phase pair;
#    with --clump-index-first the other file's entry for the same variant is a
#    genuine member, so it is not excluded.
awk '
    FNR == NR {
        if (!/^#/) {
            if ((index($5, "/") != 0) && ($4 > best[$1])) { best[$1] = $4 }
        }
        next
    }
    /^#/ { next }
    {
        ++n;
        if ($1 in picked) { print "more than one best row for " $1; exit 1 }
        picked[$1] = 1;
        if ($4 < best[$1] - 1e-9) { print $1 ": best names r^2 " $4 ", but " best[$1] " occurs"; exit 1 }
    }
    END {
        if (!n) { print "no best rows checked"; exit 1 }
        print n " best rows are maximal"
    }
' plink2.clumps.verbose plink2.clumps.best

# 4. --clump-annotate copies the named input columns into both reports.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump $FILES $P1 --clump-verbose --clump-best --clump-annotate CHISQ,OR --out plink2_ann
head -n 1 plink2_ann.clumps.verbose | grep -qx '#INDEX_ID	ID	KB	RSQ	ALLELES	F	P	CHISQ	OR'
head -n 1 plink2_ann.clumps.best | grep -qx '#INDEX_ID	ID	KB	RSQ	ALLELES	F	P	CHISQ	OR'
test "$(wc -l < plink2_ann.clumps.verbose)" -eq "$(wc -l < plink2.clumps.verbose)"

# 5. An annotation field that no input file has is tolerated, as in PLINK 1.9:
#    the column appears with every value missing rather than erroring out.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump $FILES $P1 --clump-verbose --clump-annotate NOSUCHFIELD --out plink2_unknown
head -n 1 plink2_unknown.clumps.verbose | grep -qx '#INDEX_ID	ID	KB	RSQ	ALLELES	F	P	NOSUCHFIELD'
test "$(awk '!/^#/ && $8 != "NA"' plink2_unknown.clumps.verbose | wc -l)" -eq 0

# 6. 'zs' compresses the extra reports too.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --clump zs $FILES $P1 --clump-verbose --clump-best --out plink2_zs
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.clumps.verbose.zst > plink2_zs.clumps.verbose
$BUILD/plink2 $EXTRA1 $EXTRA2 --zst-decompress plink2_zs.clumps.best.zst > plink2_zs.clumps.best
diff -q plink2.clumps.verbose plink2_zs.clumps.verbose
diff -q plink2.clumps.best plink2_zs.clumps.best
