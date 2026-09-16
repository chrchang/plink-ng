# Compares PLINK 1.9's verbose clump block, which is embedded in its .clumped
# file, against plink2's separate .clumps.verbose table.
#
# 1.9 writes "(INDEX) <id> <kb> <rsq> <alleles> <f> <p>" followed by one
# "<id> <kb> <rsq> <alleles> <f> <p>" line per member; plink2 writes
# "#INDEX_ID ID KB RSQ ALLELES F P" rows, the index variant included.  The
# in-phase allele pair is written in the opposite order by the two, so it is
# normalized.
#
# PLINK 1.9's verbose mode also reports every index candidate as its own block
# rather than only the clumps it formed -- its table grows from 14 rows to 52
# on this fixture -- so the comparison covers the clumps plink2 reports and
# ignores 1.9's extra blocks.
function abs(x) { return (x < 0)? -x : x }
function printed_tol(s,   mant, expo, dot, digits, epos) {
    epos = index(tolower(s), "e");
    if (epos) { mant = substr(s, 1, epos - 1); expo = substr(s, epos + 1) + 0 }
    else { mant = s; expo = 0 }
    dot = index(mant, ".");
    digits = dot? (length(mant) - dot) : 0;
    return 0.5 * exp((expo - digits) * log(10));
}
function same_value(a, b) { return abs(a - b) <= printed_tol(a) + 1e-12 }
function norm_alleles(s,   a) {
    if (index(s, "/") == 0) { return s }
    split(s, a, "/");
    return (a[1] <= a[2])? (a[1] "/" a[2]) : (a[2] "/" a[1]);
}
# The same variant ID can appear twice under one index -- as the index row and
# again as the other file's entry for it -- so the source file is part of the
# key rather than a compared field.
FNR == NR {
    if ($1 == "(INDEX)") { cur_index = $2; key = cur_index "|" $2 "|" $6; kb[key] = $3; rsq[key] = $4; al[key] = norm_alleles($5); p[key] = $7; ++n1; next }
    if (cur_index != "" && NF == 6 && $1 != "RANGE:" && $1 != "SPAN:") {
        key = cur_index "|" $1 "|" $5; kb[key] = $2; rsq[key] = $3; al[key] = norm_alleles($4); p[key] = $6; ++n1; next
    }
    if (NF == 0 || $1 == "RANGE:" || $1 == "SPAN:") { next }
    next
}
/^#/ { next }
{
    ++n2;
    key = $1 "|" $2 "|" $6;
    if (!(key in kb)) { print "row missing from the PLINK 1.9 report: " key; failed = 1; exit 1 }
    # The index variant's own row names a single allele, and the two programs
    # name different ones (1.9's A1 against plink2's ALT), so only the member
    # rows' in-phase pair is comparable.
    if (kb[key] != $3 || (index(al[key], "/") && (al[key] != norm_alleles($5)))) {
        print "row differs on " key ": " kb[key] "/" al[key] " vs " $3 "/" norm_alleles($5);
        failed = 1; exit 1
    }
    if (!same_value(rsq[key], $4)) { print "r^2 differs on " key ": " rsq[key] " vs " $4; failed = 1; exit 1 }
    if ((key in p) && !same_value(p[key], $7)) { print "P differs on " key ": " p[key] " vs " $7; failed = 1; exit 1 }
}
END {
    if (failed) { exit 1 }
    if (n2 == 0) { print "no rows compared"; exit 1 }
    print n2 " verbose rows matched"
}
