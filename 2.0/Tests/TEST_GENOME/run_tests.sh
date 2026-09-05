#!/bin/bash

# --genome, checked against PLINK 1.9's IBD report.
#
# PLINK 1.9 prints four decimals, so the comparison in compare.awk uses a 1e-4
# tolerance on the floating-point columns; RT and PHE have to match exactly.

set -exo pipefail

plink --simulate simulate.txt --simulate-missing 0.02 --simulate-ncases 100 --simulate-ncontrols 100 --out tmp_data > /dev/null

# Give the samples a pedigree, so that the RT column has something to report:
# 50 families of (father, mother, two children).  The phenotype column is kept.
awk 'BEGIN{OFS=" "} {
    fam = int((NR - 1) / 4);
    slot = (NR - 1) % 4;
    $1 = "F" fam;
    if (slot == 0)      { $2 = "pa"; $3 = 0;    $4 = 0;    $5 = 1 }
    else if (slot == 1) { $2 = "ma"; $3 = 0;    $4 = 0;    $5 = 2 }
    else if (slot == 2) { $2 = "c1"; $3 = "pa"; $4 = "ma"; $5 = 1 }
    else                { $2 = "c2"; $3 = "pa"; $4 = "ma"; $5 = 2 }
    print
}' tmp_data.fam > tmp_ped.fam
mv tmp_ped.fam tmp_data.fam

compare() {
    awk -f compare.awk "$1" "$2"
}

# 1. Defaults.
plink --bfile tmp_data --genome --out plink19
$1/plink2 $2 $3 --bfile tmp_data --genome --out plink2
compare plink19.genome plink2.genome

# 2. 'unbounded': (Z0, Z1, Z2) is left outside the simplex.
plink --bfile tmp_data --genome unbounded --out plink19_unb
$1/plink2 $2 $3 --bfile tmp_data --genome unbounded --out plink2_unb
compare plink19_unb.genome plink2_unb.genome

# 3. 'nudge'.
plink --bfile tmp_data --genome nudge --out plink19_nudge
$1/plink2 $2 $3 --bfile tmp_data --genome nudge --out plink2_nudge
compare plink19_nudge.genome plink2_nudge.genome

# 4. --nonfounders, which widens the allele-frequency sample set.
plink --bfile tmp_data --genome --nonfounders --out plink19_nf
$1/plink2 $2 $3 --bfile tmp_data --genome --nonfounders --out plink2_nf
compare plink19_nf.genome plink2_nf.genome

# 5. PI_HAT range filters.
plink --bfile tmp_data --genome --min 0.2 --out plink19_min
$1/plink2 $2 $3 --bfile tmp_data --genome --genome-min-pi-hat 0.2 --out plink2_min
compare plink19_min.genome plink2_min.genome
plink --bfile tmp_data --genome --max 0.15 --out plink19_max
$1/plink2 $2 $3 --bfile tmp_data --genome --genome-max-pi-hat 0.15 --out plink2_max
compare plink19_max.genome plink2_max.genome

# 6. IBS counts, from PLINK 1.9's "--genome full".
plink --bfile tmp_data --genome full --out plink19_full
$1/plink2 $2 $3 --bfile tmp_data --genome cols=+nsnp,+ibs,+homhom,+hethet --out plink2_full
awk '
    function pairkey(f1, i1, f2, i2,   a, b) {
        a = f1 "\t" i1; b = f2 "\t" i2;
        return (a < b)? (a "|" b) : (b "|" a);
    }
    FNR == NR {
        if (FNR > 1) { k = pairkey($1, $2, $3, $4); ibs0[k] = $15; ibs1[k] = $16; ibs2[k] = $17; ++n1 }
        next
    }
    /^#/ { next }
    {
        ++n2;
        k = pairkey($1, $2, $3, $4);
        # #FID1 IID1 FID2 IID2 RT Z0 Z1 Z2 PI_HAT PHE DST NSNP IBS0 IBS1 IBS2 HOMHOM HETHET
        if (ibs0[k] != $13 || ibs1[k] != $14 || ibs2[k] != $15) {
            print "IBS mismatch on " k; failed = 1; exit 1
        }
        if ($12 != $13 + $14 + $15) { print "NSNP is not IBS0+IBS1+IBS2 on " k; failed = 1; exit 1 }
        if ($12 != $16 + $17 + $14) { print "NSNP is not HOMHOM+HETHET+IBS1 on " k; failed = 1; exit 1 }
    }
    END {
        if (failed) { exit 1 }
        if (n1 != n2) { print "pair count mismatch: " n1 " vs " n2; exit 1 }
        print n1 " IBS count triples matched"
    }
' plink19_full.genome plink2_full.genome

# 7. --threads 1 must reproduce the multithreaded result.
$1/plink2 $2 $3 --bfile tmp_data --genome --threads 1 --out plink2_st
diff -q plink2.genome plink2_st.genome

# 8. --parallel: the chunks concatenate to the whole report.
for i in 1 2 3
do
    $1/plink2 $2 $3 --bfile tmp_data --genome --parallel $i 3 --out plink2_par$i
done
cat plink2_par1.genome.1 plink2_par2.genome.2 plink2_par3.genome.3 > plink2_par.genome
diff -q plink2.genome plink2_par.genome

# 9. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --genome zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.genome.zst > plink2_zs.genome
diff -q plink2.genome plink2_zs.genome
