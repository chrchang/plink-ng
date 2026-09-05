# Compares a PLINK 1.9 .tdt report against a plink2 one.
#   1.9:     CHR SNP BP A1 A2 T U OR CHISQ P
#   plink2:  #CHROM POS ID REF ALT PROVISIONAL_REF? A1 T U OR CHISQ P
# The exact tests drop CHISQ, so the column index of P is passed in.
function abs(x) { return (x < 0)? -x : x }
function same_value(a, b) {
    if (a == "NA" || b == "NA") { return a == b }
    return abs(a - b) <= 1e-6 + 1e-3 * abs(a)
}
FNR == NR {
    if (FNR > 1) {
        t[$2] = $6; u[$2] = $7; orat[$2] = $8;
        chisq[$2] = (p19chisq)? $9 : "";
        pval[$2] = $(p19col);
        ++n1;
    }
    next
}
/^#/ { next }
{
    ++n2;
    id = $3;
    if (!(id in t)) { print "variant missing from the PLINK 1.9 report: " id; failed = 1; exit 1 }
    if (t[id] != $8 || u[id] != $9) {
        print "T/U mismatch on " id ": " t[id] "/" u[id] " vs " $8 "/" $9; failed = 1; exit 1
    }
    if (!same_value(orat[id], $10)) { print "OR mismatch on " id ": " orat[id] " vs " $10; failed = 1; exit 1 }
    if (p19chisq && !same_value(chisq[id], $11)) {
        print "CHISQ mismatch on " id ": " chisq[id] " vs " $11; failed = 1; exit 1
    }
    if (!same_value(pval[id], $(p2col))) {
        print "P mismatch on " id ": " pval[id] " vs " $(p2col); failed = 1; exit 1
    }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "variant count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no variants compared"; exit 1 }
    print n1 " variants matched"
}
