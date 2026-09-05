# Compares a PLINK 1.9 .genome report against a plink2 one.  PLINK 1.9 lists
# each pair in ascending sample order and plink2 in descending, so the pair key
# is normalized; PLINK 1.9 prints four decimals, so the floating-point columns
# are compared with a tolerance.
function abs(x) { return (x < 0)? -x : x }
function close_enough(a, b) { return abs(a - b) <= 1e-4 }
function pairkey(f1, i1, f2, i2,   a, b) {
    a = f1 "\t" i1;
    b = f2 "\t" i2;
    return (a < b)? (a "|" b) : (b "|" a);
}
FNR == NR {
    if (FNR > 1) {
        # FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST ...
        k = pairkey($1, $2, $3, $4);
        rt[k] = $5; z0[k] = $7; z1[k] = $8; z2[k] = $9; pihat[k] = $10;
        phe[k] = $11; dst[k] = $12;
        ++n1;
    }
    next
}
/^#/ { next }
{
    # FID1 IID1 FID2 IID2 RT Z0 Z1 Z2 PI_HAT PHE DST
    ++n2;
    k = pairkey($1, $2, $3, $4);
    if (!(k in rt)) { print "pair missing from the PLINK 1.9 report: " k; failed = 1; exit 1 }
    if (rt[k] != $5) { print "RT mismatch on " k ": " rt[k] " vs " $5; failed = 1; exit 1 }
    if (phe[k] != $10) { print "PHE mismatch on " k ": " phe[k] " vs " $10; failed = 1; exit 1 }
    if (!close_enough(z0[k], $6) || !close_enough(z1[k], $7) || !close_enough(z2[k], $8) ||
        !close_enough(pihat[k], $9) || !close_enough(dst[k], $11)) {
        print "statistic mismatch on " k ": " z0[k] "/" z1[k] "/" z2[k] "/" pihat[k] "/" dst[k] \
              " vs " $6 "/" $7 "/" $8 "/" $9 "/" $11;
        failed = 1; exit 1
    }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "pair count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no pairs compared"; exit 1 }
    print n1 " pairs matched"
}
