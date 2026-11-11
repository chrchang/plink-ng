# Compares a PLINK 1.9 --r2 dprime report (first file) with a plink2
# --r2-phased cols=+dprime report (second file).  Pairs are matched by their
# two IDs in either order; PLINK 1.9's self-pairs (from --ld-snp-list) are
# skipped, and a pair it lists in both orders is counted once.  Every pair
# must appear in both reports, with R2 and |D'| within 1e-4, and at least one
# compared pair must involve chrX.
function abs(x) { return (x < 0)? -x : x }
function pairkey(a, b) { return (a < b)? (a " " b) : (b " " a) }
FNR == NR {
  if (FNR > 1 && $3 != $6) { key = pairkey($3, $6); if (!(key in r2)) ++n1; r2[key] = $7; dp[key] = $8 }
  next
}
/^#/ { next }
{
  ++n2
  key = pairkey($3, $6)
  if (!(key in r2)) { print "missing from PLINK 1.9 report: " key; bad = 1; next }
  if (abs(r2[key] - $7) > 1e-4 || abs(abs(dp[key]) - abs($8)) > 1e-4) {
    if (++nbad <= 5) print "differs: " key " " r2[key] " " dp[key] " vs " $7 " " $8
    bad = 1
  }
  if ($1 == "X" || $4 == "X") ++nx
}
END {
  if (n1 != n2) { print "pair count mismatch: " n1 " vs " n2; bad = 1 }
  if (!nx) { print "no chrX pairs compared"; bad = 1 }
  exit bad
}
