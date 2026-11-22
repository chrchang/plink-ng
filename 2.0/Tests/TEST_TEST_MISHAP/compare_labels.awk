# Compares every haplotype row of a PLINK 1.9 .missing.hap report against
# plink2's by (variant ID, HAPLOTYPE label), so that a row attached to the
# wrong allele names fails.  Rows where either side is numerically degenerate
# (an EM count like 3.99e-15/5.76e-14 instead of 0/0) are skipped, as are
# loci only one program reports; at least 95% of plink2's rows must still be
# compared.
#
#   1.9:     SNP HAPLOTYPE F_0 F_1 M_H1 M_H2 CHISQ P FLANKING
#   plink2:  #ID HAPLOTYPE F_0 F_1 M_H1 M_H2 CHISQ P FLANKING
function degenerate(s) { return s ~ /e-/ }
FNR == NR {
    if (FNR > 1) { mh1[$1 SUBSEP $2] = $5; mh2[$1 SUBSEP $2] = $6 }
    next
}
FNR > 1 {
    ++total;
    key = $1 SUBSEP $2;
    if (!(key in mh1)) { next }
    if (degenerate(mh1[key]) || degenerate(mh2[key]) || degenerate($5) || degenerate($6)) { next }
    ++compared;
    if (mh1[key] != $5 || mh2[key] != $6) {
        print $1 " " $2 ": PLINK 1.9 has " mh1[key] " " mh2[key] ", plink2 has " $5 " " $6;
        exit 1
    }
}
END {
    if (compared < 0.95 * total) { print "only " compared " of " total " rows compared"; exit 1 }
    print compared " of " total " rows match PLINK 1.9"
}
