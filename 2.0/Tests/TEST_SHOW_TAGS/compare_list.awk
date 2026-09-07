# Compares a PLINK 1.9 .tags.list against a plink2 one.
#   1.9:     SNP CHR BP NTAG LEFT RIGHT KBSPAN TAGS
#   plink2:  #CHROM POS ID NTAG LEFT RIGHT KBSPAN TAGS
function abs(x) { return (x < 0)? -x : x }
FNR == NR {
    if (FNR > 1) { ntag[$1] = $4; left[$1] = $5; right[$1] = $6; kb[$1] = $7; tags[$1] = $8; ++n1 }
    next
}
/^#/ { next }
{
    ++n2;
    id = $3;
    if (!(id in ntag)) { print "variant missing from the PLINK 1.9 report: " id; failed = 1; exit 1 }
    if (ntag[id] != $4 || left[id] != $5 || right[id] != $6 || tags[id] != $8) {
        print "tag set differs on " id ": " ntag[id] "/" left[id] "/" right[id] "/" tags[id] \
              " vs " $4 "/" $5 "/" $6 "/" $8;
        failed = 1; exit 1
    }
    if (abs(kb[id] - $7) > 1e-3) { print "KBSPAN differs on " id ": " kb[id] " vs " $7; failed = 1; exit 1 }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "variant count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no variants compared"; exit 1 }
    print n1 " variants matched"
}
