# Compares a PLINK 1.9 .tags.list against a plink2 one.
#   1.9:     SNP CHR BP NTAG LEFT RIGHT KBSPAN TAGS
#   plink2:  #ID CHROM POS NTAG LEFT RIGHT KBSPAN TAGS
function abs(x) { return (x < 0)? -x : x }
FNR == NR {
    if (FNR > 1) { ntag[$1] = $4; left[$1] = $5; right[$1] = $6; kb[$1] = $7; tags[$1] = $8; ++n1 }
    next
}
/^#/ { next }
{
    ++n2;
    if (!($1 in ntag)) { print "variant missing from the PLINK 1.9 report: " $1; failed = 1; exit 1 }
    if (ntag[$1] != $4 || left[$1] != $5 || right[$1] != $6 || tags[$1] != $8) {
        print "tag set differs on " $1 ": " ntag[$1] "/" left[$1] "/" right[$1] "/" tags[$1] \
              " vs " $4 "/" $5 "/" $6 "/" $8;
        failed = 1; exit 1
    }
    if (abs(kb[$1] - $7) > 1e-3) { print "KBSPAN differs on " $1 ": " kb[$1] " vs " $7; failed = 1; exit 1 }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "variant count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no variants compared"; exit 1 }
    print n1 " variants matched"
}
