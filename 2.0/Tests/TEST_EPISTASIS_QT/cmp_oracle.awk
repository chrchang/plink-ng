# Compares an oracle.py table against a plink2 .epi.qt report.  plink2 prints
# six significant digits, so the tolerance is half of the last printed place.
function abs(x) { return (x < 0)? -x : x }
function key(a, b) { return (a < b)? (a "|" b) : (b "|" a) }
function bad(name, k, want, got) {
    if (abs(want - got) <= 5.5e-6 * (1 + abs(want))) { return 0 }
    print name " differs on " k ": " want " vs " got
    return 1
}
FNR == NR { if (FNR > 1) { k = key($1, $2); b[k] = $3 + 0; se[k] = $4 + 0; t[k] = $5 + 0; ++n1 }; next }
/^#/ { next }
{
    ++n2
    k = key($2, $4)
    if (!(k in b)) { print "pair missing from the oracle: " k; exit 1 }
    e += bad("BETA_INT", k, b[k], $5 + 0)
    e += bad("SE", k, se[k], $6 + 0)
    e += bad("T_STAT", k, t[k], $7 + 0)
    if (e > 3) { exit 1 }
}
END {
    if (n1 != n2) { print "pair count mismatch: " n1 " vs " n2; exit 1 }
    if (e) { exit 1 }
    print n1 " pairs matched the oracle"
}
