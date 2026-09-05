#!/bin/bash

# --homozyg, checked against PLINK 1.9's .hom, .hom.indiv and .hom.summary.
#
# plink2 requires --homozyg-min-af, which PLINK 1.9 has no equivalent of, so
# the comparison runs pass 0 (no frequency floor, which is what 1.9 does).

set -exo pipefail

awk -f make_vcf.awk > tmp_h.vcf
plink --vcf tmp_h.vcf --double-id --make-bed --out tmp_data > /dev/null

# PLINK 1.9 writes fixed-precision columns (3 decimals), so the floating-point
# columns are compared with a tolerance rather than exactly.  Run boundaries
# and counts have to match exactly.
compare_hom() {
    awk '
        function abs(x) { return (x < 0)? -x : x }
        function close_enough(a, b) { return abs(a - b) <= 2e-3 + 2e-3 * abs(a) }
        FNR == NR {
            if (FNR > 1) {
                # FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET
                key = $1 "_" $2 "_" $5;
                seen[key] = 1;
                snp2[key] = $6; pos1[key] = $7; pos2[key] = $8;
                kb[key] = $9; nsnp[key] = $10; density[key] = $11;
                phom[key] = $12; phet[key] = $13;
                ++n1;
            }
            next
        }
        /^#/ { next }
        {
            # FID IID CHROM ID1 ID2 POS1 POS2 KB NSNP DENSITY PHOM PHET
            ++n2;
            key = $1 "_" $2 "_" $4;
            if (!(key in seen)) { print "run not found in PLINK 1.9 output: " key; failed = 1; exit 1 }
            if (snp2[key] != $5 || pos1[key] != $6 || pos2[key] != $7 || nsnp[key] != $9) {
                print "boundary mismatch on " key; failed = 1; exit 1
            }
            if (!close_enough(kb[key], $8) || !close_enough(density[key], $10) ||
                !close_enough(phom[key], $11) || !close_enough(phet[key], $12)) {
                print "statistic mismatch on " key; failed = 1; exit 1
            }
        }
        END {
            if (failed) { exit 1 }
            if (n1 != n2) { print "run count mismatch: " n1 " vs " n2; exit 1 }
            if (n1 == 0) { print "no runs found, so nothing was compared"; exit 1 }
            print n1 " runs matched"
        }
    ' "$1" "$2"
}

compare_indiv() {
    awk '
        function abs(x) { return (x < 0)? -x : x }
        function close_enough(a, b) { return abs(a - b) <= 2e-3 + 2e-3 * abs(a) }
        FNR == NR {
            if (FNR > 1) { nseg[$2] = $4; kb[$2] = $5; kbavg[$2] = $6; ++n1 }
            next
        }
        /^#/ { next }
        {
            ++n2;
            if (nseg[$2] != $3 || !close_enough(kb[$2], $4) || !close_enough(kbavg[$2], $5)) {
                print "per-sample mismatch on " $2; failed = 1; exit 1
            }
        }
        END {
            if (failed) { exit 1 }
            if (n1 != n2) { print "sample count mismatch: " n1 " vs " n2; exit 1 }
            print n1 " samples matched"
        }
    ' "$1" "$2"
}

compare_summary() {
    awk '
        FNR == NR {
            if (FNR > 1) { aff[$2] = $4; unaff[$2] = $5; ++n1 }
            next
        }
        /^#/ { next }
        {
            ++n2;
            if (aff[$3] != $4 || unaff[$3] != $5) {
                print "summary mismatch on " $3; failed = 1; exit 1
            }
        }
        END {
            if (failed) { exit 1 }
            if (n1 != n2) { print "variant count mismatch: " n1 " vs " n2; exit 1 }
            print n1 " variants matched"
        }
    ' "$1" "$2"
}

# 1. Defaults.
plink --bfile tmp_data --homozyg --out plink19
$1/plink2 $2 $3 --bfile tmp_data --homozyg --homozyg-min-af 0 --out plink2
compare_hom plink19.hom plink2.hom
compare_indiv plink19.hom.indiv plink2.hom.indiv
compare_summary plink19.hom.summary plink2.hom.summary

# 2. Non-default thresholds, to move the run boundaries around.
plink --bfile tmp_data --homozyg --homozyg-snp 50 --homozyg-kb 500 --homozyg-window-snp 30 --homozyg-window-het 0 --homozyg-gap 200 --out plink19_tight
$1/plink2 $2 $3 --bfile tmp_data --homozyg --homozyg-min-af 0 --homozyg-snp 50 --homozyg-kb 500 --homozyg-window-snp 30 --homozyg-window-het 0 --homozyg-gap 200 --out plink2_tight
compare_hom plink19_tight.hom plink2_tight.hom
compare_indiv plink19_tight.hom.indiv plink2_tight.hom.indiv
compare_summary plink19_tight.hom.summary plink2_tight.hom.summary

# 3. --homozyg-het, which caps the heterozygous calls inside a reported run.
plink --bfile tmp_data --homozyg --homozyg-het 2 --out plink19_het
$1/plink2 $2 $3 --bfile tmp_data --homozyg --homozyg-min-af 0 --homozyg-het 2 --out plink2_het
compare_hom plink19_het.hom plink2_het.hom

# 4. --homozyg-min-af drops rare variants, so it can only shorten runs, and
#    it is required rather than defaulted.
$1/plink2 $2 $3 --bfile tmp_data --homozyg --homozyg-min-af 0.05 --out plink2_maf
test "$(grep -vc '^#' plink2_maf.hom)" -le "$(grep -vc '^#' plink2.hom)"
# --homozyg-maf is the same flag.
$1/plink2 $2 $3 --bfile tmp_data --homozyg --homozyg-maf 0.05 --out plink2_maf2
diff -q plink2_maf.hom plink2_maf2.hom
if $1/plink2 $2 $3 --bfile tmp_data --homozyg --out plink2_nomaf > /dev/null 2>&1; then
    echo "--homozyg ran without --homozyg-min-af"
    exit 1
fi

# 5. Multithreading must not change the result.
$1/plink2 $2 $3 --bfile tmp_data --homozyg --homozyg-min-af 0 --threads 1 --out plink2_st
diff -q plink2.hom plink2_st.hom
diff -q plink2.hom.indiv plink2_st.hom.indiv
diff -q plink2.hom.summary plink2_st.hom.summary

# 6. 'zs' compresses the .hom.summary file, and nothing else.
$1/plink2 $2 $3 --bfile tmp_data --homozyg zs --homozyg-min-af 0 --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.hom.summary.zst > plink2_zs.hom.summary
diff -q plink2.hom.summary plink2_zs.hom.summary
diff -q plink2.hom plink2_zs.hom

# 7. cols= drops the columns it doesn't name.
$1/plink2 $2 $3 --bfile tmp_data --homozyg cols=chrom,pos,nsnp --homozyg-min-af 0 --out plink2_cols
head -n 1 plink2_cols.hom | grep -qx '#IID	CHROM	ID1	ID2	POS1	POS2	NSNP'
