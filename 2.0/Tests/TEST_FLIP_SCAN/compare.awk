# Compares a PLINK 1.9 .flipscan report against a plink2 one.
#   1.9:     CHR SNP BP A1 A2 F POS R_POS NEG R_NEG NEGSNPS
#   plink2:  #CHROM POS ID REF ALT ALT_FREQ POS_CT R_POS NEG_CT R_NEG NEG_IDS
# PLINK 1.9's F is the A1 (minor allele) frequency and plink2's is the ALT
# frequency, so only the minor-allele form is comparable.
function abs(x) { return (x < 0)? -x : x }
function minor(f) { return (f <= 0.5)? f : (1 - f) }
function same_value(a, b) {
    if (a == "NA" || b == "NA") { return a == b }
    return abs(a - b) <= 1e-3 + 1e-3 * abs(a)
}
FNR == NR {
    if (FNR > 1) {
        f[$2] = $6; posct[$2] = $7; rpos[$2] = $8; negct[$2] = $9; rneg[$2] = $10;
        negids[$2] = ($11 == "")? "." : $11;
        ++n1;
    }
    next
}
/^#/ { next }
{
    ++n2;
    id = $3;
    if (!(id in posct)) { print "variant missing from the PLINK 1.9 report: " id; failed = 1; exit 1 }
    if (posct[id] != $7 || negct[id] != $9) {
        print "neighbor counts differ on " id ": " posct[id] "/" negct[id] " vs " $7 "/" $9; failed = 1; exit 1
    }
    if (negids[id] != $11) { print "flipped-neighbor list differs on " id ": " negids[id] " vs " $11; failed = 1; exit 1 }
    if (!same_value(minor(f[id]), minor($6))) { print "frequency differs on " id ": " f[id] " vs " $6; failed = 1; exit 1 }
    if (!same_value(rpos[id], $8) || !same_value(rneg[id], $10)) {
        print "mean correlation differs on " id ": " rpos[id] "/" rneg[id] " vs " $8 "/" $10; failed = 1; exit 1
    }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "variant count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no variants compared"; exit 1 }
    print n1 " variants matched"
}
