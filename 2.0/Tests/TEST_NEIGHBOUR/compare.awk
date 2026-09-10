# Compares a PLINK 1.9 .nearest report against a plink2 one.
#
# PLINK 1.9 prints four significant digits, so IBS and Z are compared with a
# tolerance.  When two candidates share a sample's IBS value exactly, which one
# is named is down to the last bits of two different rescalings, so a row whose
# numbers agree but whose neighbour differs is counted as a tie rather than a
# mismatch; the test caps how many of those it will accept.
function abs(x) { return (x < 0)? -x : x }
function close_enough(a, b) { return abs(a - b) <= 1e-4 + 1e-4 * abs(a) }
FNR == NR {
    if (FNR > 1) {
        # FID IID NN MIN_DST Z FID2 IID2
        k = $1 "\t" $2 "\t" $3;
        dst[k] = $4; z[k] = $5; nbr[k] = $6 "\t" $7;
        ++n1;
    }
    next
}
/^#/ { next }
{
    # FID IID NN IBS Z FID2 IID2
    ++n2;
    k = $1 "\t" $2 "\t" $3;
    if (!(k in dst)) { print "row missing from the 1.9 report: " k; exit 1 }
    if (!close_enough(dst[k], $4)) {
        print "IBS mismatch on " k ": " dst[k] " vs " $4; exit 1
    }
    if (z[k] == "nan" || z[k] == "-nan" || $5 == "NA") {
        # Undefined Z: 1.9 prints a NaN or a zero from 0/0, plink2 prints NA.
        # Both sides have to see it as undefined.
        if (!(($5 == "NA") && (z[k] == "nan" || z[k] == "-nan" || z[k] + 0 == 0))) {
            print "Z definedness mismatch on " k ": " z[k] " vs " $5; exit 1
        }
    } else if (abs(z[k] - $5) > 2e-3) {
        print "Z mismatch on " k ": " z[k] " vs " $5; exit 1
    }
    if (nbr[k] != $6 "\t" $7) { ++ties }
}
END {
    if (n1 != n2) { print "row count mismatch: " n1 " vs " n2; exit 1 }
    if (ties > n2 / 20) { print ties " tied rows out of " n2 ", too many to be ties"; exit 1 }
    print n2 " rows matched (" ties+0 " ties)"
}
