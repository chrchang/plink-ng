# Compares two .mds files dimension by dimension, allowing each dimension's
# sign to differ.  Both files must list the same samples in the same order,
# and values must agree within 1e-4 + 1e-3 * |value|.
function abs(x) { return (x < 0)? -x : x }
FNR == 1 { ++fileno; nd = NF - 3; next }
fileno == 1 { id[FNR] = $1 " " $2; for (d = 1; d <= nd; ++d) a[FNR, d] = $(d + 3); n = FNR; next }
{
  if (id[FNR] != $1 " " $2) { print "sample mismatch on line " FNR; bad = 1 }
  for (d = 1; d <= nd; ++d) b[FNR, d] = $(d + 3)
  m = FNR
}
END {
  if (n < 3 || nd < 1 || m != n) { print "nothing to compare, or row count mismatch"; exit 1 }
  for (d = 1; d <= nd; ++d) {
    # take the sign from the largest-magnitude entry
    best = 2
    for (k = 3; k <= n; ++k) if (abs(a[k, d]) > abs(a[best, d])) best = k
    sgn = (a[best, d] * b[best, d] < 0)? -1 : 1
    for (k = 2; k <= n; ++k) {
      if (abs(a[k, d] - sgn * b[k, d]) > 1e-4 + 1e-3 * abs(a[k, d])) {
        if (++nbad <= 3) print "dimension " d ", line " k ": " a[k, d] " vs " b[k, d]
        bad = 1
      }
    }
  }
  exit bad
}
