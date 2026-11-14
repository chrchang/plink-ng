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
            # FID IID CHROM POS1 POS2 ID1 ID2 KB NSNP DENSITY PHOM PHET
            ++n2;
            key = $1 "_" $2 "_" $6;
            if (!(key in seen)) { print "run not found in PLINK 1.9 output: " key; failed = 1; exit 1 }
            if (snp2[key] != $7 || pos1[key] != $4 || pos2[key] != $5 || nsnp[key] != $9) {
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
            if (FNR > 1) { rohct[$2] = $5; ++n1 }
            next
        }
        /^#/ { next }
        {
            ++n2;
            if (rohct[$3] != $4) {
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

compare_summary_cc() {
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
# --homozyg-min-af is only demanded when the data actually contains
# low-frequency variants.  This fileset's lowest MAF is 0.25, so the flag is
# optional here, and omitting it must match an explicit floor of 0.
$1/plink2 $2 $3 --bfile tmp_data --homozyg --out plink2_nomaf
$1/plink2 $2 $3 --bfile tmp_data --homozyg --homozyg-min-af 0 --out plink2_maf0
diff -q plink2_nomaf.hom plink2_maf0.hom
diff -q plink2_nomaf.hom.indiv plink2_maf0.hom.indiv

# With low-frequency variants present it is demanded again.
{ printf '400 rare 0.001 0.02 1 1\n'; printf '400 common 0.2 0.5 1 1\n'; } > tmp_rare.sim
plink --simulate tmp_rare.sim --simulate-ncases 50 --simulate-ncontrols 50 --out tmp_rare > /dev/null
if $1/plink2 $2 $3 --bfile tmp_rare --homozyg --out plink2_nomaf2 > /dev/null 2>&1; then
    echo "--homozyg ran without --homozyg-min-af on a fileset with rare variants"
    exit 1
fi
$1/plink2 $2 $3 --bfile tmp_rare --homozyg --homozyg-min-af 0.05 --out plink2_rare

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
head -n 1 plink2_cols.hom | grep -qx '#IID	CHROM	POS1	POS2	ID1	ID2	NSNP'
# cols= named no .hom.indiv column, so only IID survives there.
head -n 1 plink2_cols.hom.indiv | grep -qx '#IID'

# 8. FROH: total autosomal run length over the span of the autosomal variants
#    actually scanned.  The fixture is one chromosome, so the denominator is
#    that chromosome's span; it is recomputed here from the .bim rather than
#    taken from plink2.
head -n 1 plink2.hom.indiv | grep -q 'FROH'
awk '
    function abs(x) { return (x < 0)? -x : x }
    FNR == NR {
        if (FNR == 1) { min = $4; max = $4 }
        if ($4 < min) { min = $4 }
        if ($4 > max) { max = $4 }
        next
    }
    FNR == 1 {
        denom = (max + 1 - min) / 1000.0
        if (denom <= 0) { print "empty denominator"; exit 1 }
        # Column positions depend on which optional columns are present, so
        # they are read from the header rather than assumed.
        for (i = 1; i <= NF; i++) {
            name = $i
            sub(/^#/, "", name)
            col[name] = i
        }
        if (!col["KB"] || !col["FROH"] || !col["IID"]) {
            print "missing expected column in .hom.indiv header"; exit 1
        }
        next
    }
    {
        want = $(col["KB"]) / denom
        if (abs($(col["FROH"]) - want) > 1e-5 + 1e-5 * abs(want)) {
            print "FROH mismatch on " $(col["IID"]) ": " $(col["FROH"]) " vs " want
            exit 1
        }
        if ($(col["FROH"]) > 0) { nonzero = 1 }
        ++n
    }
    END {
        if (n != 30) { print "expected 30 samples, saw " n; exit 1 }
        if (!nonzero) { print "every FROH is zero"; exit 1 }
        print n " FROH values verified"
    }
' tmp_data.bim plink2.hom.indiv

# 9. cols= without froh drops the column.
$1/plink2 $2 $3 --bfile tmp_data --homozyg cols=nseg,kbtot --homozyg-min-af 0 --out plink2_nofroh
head -n 1 plink2_nofroh.hom.indiv | grep -qx '#IID	NSEG	KB'

# 10. pheno= stratifies the .hom.summary ROH_CT column by category, with the
#     category names in the header; the counts must add up to the unstratified
#     ROH_CT (every sample has a binary phenotype value here).
awk '{print $1 "\t" $2 "\t" (NR % 2? 2 : 1)}' tmp_data.fam > tmp_pheno.txt
$1/plink2 $2 $3 --bfile tmp_data --pheno tmp_pheno.txt --homozyg pheno=PHENO1 --homozyg-min-af 0 --out plink2_pheno
head -n 1 plink2_pheno.hom.summary | grep -qx '#CHROM	POS	ID	CASE	CONTROL'
awk '
    FNR == NR { if (FNR > 1) { tot[$3] = $4 } next }
    FNR > 1 {
        if ($4 + $5 != tot[$3]) { print "stratified counts do not sum to ROH_CT at " $3; exit 1 }
        if ($4 > 0 || $5 > 0) { nonzero = 1 }
    }
    END { if (!nonzero) { print "every stratified count is zero"; exit 1 } }
' plink2.hom.summary plink2_pheno.hom.summary

# 11. With no runs at all, every .hom.summary count must be zero, even when the
#     workspace was not clean (--randmem), which small inputs otherwise hide.
randmem=--randmem
[[ "$2 $3" == *--randmem* ]] && randmem=
$1/plink2 $2 $3 --dummy 100 200 0.01 acgt --seed 1 --make-pgen --out tmp_noroh
$1/plink2 $2 $3 $randmem --pfile tmp_noroh --homozyg --homozyg-min-af 0 --out plink2_noroh
awk 'FNR > 1 && $4 != 0 { print "nonzero ROH_CT at " $3; exit 1 }' plink2_noroh.hom.summary

# 12. --homozyg-window-missing, --homozyg-window-threshold and
#     --homozyg-density, against PLINK 1.9.  This variant of the fixture has
#     ~3% of calls set missing, and the second half of the chromosome is
#     spaced 3 kb apart, so the three runs there have densities of 2.3-3 kb
#     per variant.  Each setting must also change the runs found, so it's
#     known to have an effect.
awk 'BEGIN { OFS = "\t" }
    /^#/ { print; next }
    {
        v = NR;
        $2 = (v <= 2500)? v * 1000 : 2500000 + (v - 2500) * 3000;
        for (i = 10; i <= NF; i++) {
            if (((NR * 31 + i * 17) % 37) == 0) { $i = "./." }
        }
        print
    }' tmp_h.vcf > tmp_hm.vcf
plink --vcf tmp_hm.vcf --double-id --make-bed --out tmp_hm > /dev/null
$1/plink2 $2 $3 --bfile tmp_hm --homozyg --homozyg-min-af 0 --out plink2_hm
for opt in "--homozyg-window-missing 2" "--homozyg-window-threshold 0.5" "--homozyg-density 2"; do
    plink --bfile tmp_hm --homozyg $opt --out plink19_opt > /dev/null
    $1/plink2 $2 $3 --bfile tmp_hm --homozyg --homozyg-min-af 0 $opt --out plink2_opt
    compare_hom plink19_opt.hom plink2_opt.hom
    compare_indiv plink19_opt.hom.indiv plink2_opt.hom.indiv
    if diff -q plink2_opt.hom plink2_hm.hom > /dev/null; then
        echo "$opt had no effect"
        exit 1
    fi
done
