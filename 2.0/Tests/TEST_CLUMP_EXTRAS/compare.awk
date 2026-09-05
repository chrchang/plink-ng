# Compares a PLINK 1.9 .clumped report against a plink2 .clumps one.
#   1.9:     CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2
#   plink2:  #CHROM POS ID F P TOTAL NONSIG S0.05 S0.01 S0.001 S0.0001 SP2
# PLINK 1.9 pads its report with blank lines and prints p-values with four
# significant digits.
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
    # Half of PLINK 1.9's last printed digit.
    tol = printed_tol(a);
    return abs(a - b) <= tol + 1e-12
}
FNR == NR {
    if (FNR == 1 || NF == 0) { next }
    id = $3;
    f[id] = $2; p[id] = $5; total[id] = $6; nsig[id] = $7;
    s05[id] = $8; s01[id] = $9; s001[id] = $10; s0001[id] = $11; sp2[id] = $12;
    ++n1;
    next
}
/^#/ { next }
NF == 0 { next }
{
    ++n2;
    id = $3;
    if (!(id in total)) { print "clump missing from the PLINK 1.9 report: " id; failed = 1; exit 1 }
    if (f[id] != $4 || total[id] != $6 || nsig[id] != $7 || s05[id] != $8 || s01[id] != $9 ||
        s001[id] != $10 || s0001[id] != $11) {
        print "clump contents differ on " id; failed = 1; exit 1
    }
    if (sp2[id] != $12) { print "SP2 differs on " id ": " sp2[id] " vs " $12; failed = 1; exit 1 }
    if (!same_value(p[id], $5)) { print "P differs on " id ": " p[id] " vs " $5; failed = 1; exit 1 }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "clump count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no clumps compared"; exit 1 }
    print n1 " clumps matched"
}
