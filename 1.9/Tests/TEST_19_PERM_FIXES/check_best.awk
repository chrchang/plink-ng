# Usage: awk [-v lo=1] -f check_best.awk dump best
# Checks that each permutation's value in a --mperm-save .mperm.dump.best
# file equals the most extreme value in that permutation's row of a
# --mperm-save-all dump from a run with the same permutations.  Row 0 (the
# original statistics) is not compared.  Set lo=1 for p-values.
NR == FNR {
  if (FNR > 1) {
    b = "NA"
    for (k = 2; k <= NF; ++k) {
      if ($k != "NA" && (b == "NA" || (lo? ($k + 0 < b + 0) : ($k + 0 > b + 0)))) {
        b = $k
      }
    }
    best[FNR - 1] = b
    ++np
  }
  next
}
FNR > 1 {
  ++n
  d = $2 - best[$1]
  if (d < 0) {
    d = -d
  }
  if (!($1 in best) || d > 1e-5 * $2 + 1e-8) {
    if (++shown <= 5) {
      print "permutation " $1 ": " $2 " vs " best[$1]
    }
    bad = 1
  }
}
END {
  if (n != np) {
    print n " permutations vs " np
    bad = 1
  }
  exit bad
}
