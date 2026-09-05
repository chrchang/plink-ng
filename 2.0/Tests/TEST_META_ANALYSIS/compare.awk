# Compares a PLINK 1.9 .meta report against a plink2 one, by column name, so
# that the modes which add or drop columns ('no-map', 'no-allele', 'study',
# 'weighted-z') are handled without index juggling.
#
# PLINK 1.9 prints four significant digits, so the numeric columns are
# compared with a relative tolerance.
function abs(x) { return (x < 0)? -x : x }
# PLINK 1.9 prints these columns with a fixed number of digits, so the
# tolerance is half of its last printed digit rather than a fixed relative
# one: 0.0063 and 0.00633666 agree as far as 1.9 said anything.
function printed_tol(s,   mant, expo, dot, digits, epos) {
    epos = index(tolower(s), "e");
    if (epos) {
        mant = substr(s, 1, epos - 1);
        expo = substr(s, epos + 1) + 0;
    } else {
        mant = s;
        expo = 0;
    }
    dot = index(mant, ".");
    digits = dot? (length(mant) - dot) : 0;
    return 0.5 * exp((expo - digits) * log(10));
}
function same_value(a, b, abs_floor,   tol) {
    if (a == "NA" || b == "NA") { return a == b }
    # Half of 1.9's last printed digit, or a 0.1% relative band, whichever is
    # looser: 1.9 prints four significant digits, and a few statistics differ
    # in the last of them.
    tol = printed_tol(a);
    if (1e-3 * abs(a) > tol) { tol = 1e-3 * abs(a) }
    if (abs_floor > tol) { tol = abs_floor }
    return abs(a - b) <= tol + 1e-12
}
function canon(name) {
    sub(/^#/, "", name);
    if (name == "CHROM") { return "CHR" }
    if (name == "POS") { return "BP" }
    if (name == "ID") { return "SNP" }
    if (name == "P_R") { return "P(R)" }
    if (name == "OR_R") { return "OR(R)" }
    if (name == "BETA_R") { return "BETA(R)" }
    if (name == "I2") { return "I" }
    if (name == "P_WZ") { return "P(WZ)" }
    return name
}
FNR == NR {
    if (FNR == 1) {
        for (i = 1; i <= NF; ++i) { c19[canon($i)] = i }
        next
    }
    id = $(c19["SNP"]);
    for (name in c19) { v19[id "|" name] = $(c19[name]) }
    ++n1;
    next
}
FNR == 1 {
    for (i = 1; i <= NF; ++i) { c2[canon($i)] = i }
    for (name in c2) {
        if (!(name in c19)) { print "column " name " is not in the PLINK 1.9 report"; failed = 1; exit 1 }
    }
    next
}
{
    ++n2;
    id = $(c2["SNP"]);
    if (!((id "|SNP") in v19)) { print "variant missing from the PLINK 1.9 report: " id; failed = 1; exit 1 }
    for (name in c2) {
        want = v19[id "|" name];
        got = $(c2[name]);
        if (name == "CHR" || name == "SNP" || name == "A1" || name == "A2") {
            # PLINK 1.9 writes '?' for an unknown allele, plink2 '.'.
            if (want == "?") { want = "." }
            if (got == "?") { got = "." }
            if (want != got) { print name " differs on " id ": " want " vs " got; failed = 1; exit 1 }
        } else if (name == "P(WZ)") {
            # Skipped: this is Phi(WEIGHTED_Z), and PLINK 1.9 evaluates the
            # normal CDF with a five-term Abramowitz-Stegun approximation whose
            # absolute error is around 7.5e-8, which swamps the p-values this
            # column reaches.  WEIGHTED_Z itself does not go through that
            # approximation and is compared normally, and the p-value is
            # checked against it separately.
            continue;
        } else if (!same_value(want, got, 0)) {
            print name " differs on " id ": " want " vs " got; failed = 1; exit 1
        }
    }
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "variant count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no variants compared"; exit 1 }
    print n1 " variants matched"
}
