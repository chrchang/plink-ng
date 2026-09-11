# Compares the HETERO rows of a PLINK 1.9 .missing.hap report against
# plink2's.  Those are the rows that do not depend on how the two programs
# label the flanking haplotypes, and PLINK 1.9's haplotype rows are not
# comparable in general: its EM step emits numerically degenerate rows (counts
# like 9.99e-16/2.02e-14) where plink2 reports 0/0.
#
#   1.9:     SNP HAPLOTYPE F_0 F_1 M_H1 M_H2 CHISQ P FLANKING
#   plink2:  #ID HAPLOTYPE F_0 F_1 M_H1 M_H2 CHISQ P FLANKING
function abs(x) { return (x < 0)? -x : x }
function printed_tol(s,   mant, expo, dot, digits, epos) {
    epos = index(tolower(s), "e");
    if (epos) { mant = substr(s, 1, epos - 1); expo = substr(s, epos + 1) + 0 }
    else { mant = s; expo = 0 }
    dot = index(mant, ".");
    digits = dot? (length(mant) - dot) : 0;
    return 0.5 * exp((expo - digits) * log(10));
}
function same_value(a, b,   tol) {
    if (a == "NA" || b == "NA") { return a == b }
    tol = printed_tol(a);
    if (1e-3 * abs(a) > tol) { tol = 1e-3 * abs(a) }
    return abs(a - b) <= tol + 1e-12
}
FNR == NR {
    if (FNR > 1 && $2 == "HETERO") {
        f0[$1] = $3; f1[$1] = $4; mh1[$1] = $5; mh2[$1] = $6; chisq[$1] = $7; pval[$1] = $8; flank[$1] = $9;
        ++n1;
    }
    next
}
/^#/ { next }
$2 == "HETERO" {
    ++n2;
    id = $1;
    if (!(id in f0)) { print "locus missing from the PLINK 1.9 report: " id; failed = 1; exit 1 }
    if (mh1[id] != $5 || mh2[id] != $6 || flank[id] != $9) {
        print "counts or flanking variants differ on " id ": " mh1[id] "/" mh2[id] "/" flank[id] \
              " vs " $5 "/" $6 "/" $9;
        failed = 1; exit 1
    }
    if (!same_value(f0[id], $3) || !same_value(f1[id], $4) ||
        !same_value(chisq[id], $7) || !same_value(pval[id], $8)) {
        print "statistic differs on " id ": " f0[id] "/" f1[id] "/" chisq[id] "/" pval[id] \
              " vs " $3 "/" $4 "/" $7 "/" $8;
        failed = 1; exit 1
    }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "HETERO row count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no HETERO rows compared"; exit 1 }
    print n1 " HETERO rows matched"
}
