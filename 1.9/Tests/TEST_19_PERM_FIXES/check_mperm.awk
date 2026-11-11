# Usage: awk [-v lo=1] -f check_mperm.awk ids dump report
#   ids     the variant IDs of the dump's columns, one per line, in order
#   dump    a --mperm-save-all .mperm.dump.all file
#   report  a .mperm report from a run with the 'perm-count' modifier
# Recomputes EMP1 (the number of permutations at least as extreme as the
# original statistic, with ties counting 1/2) and EMP2 (the number of
# permutations whose most extreme statistic is at least as extreme) from the
# dump, and compares them with the report, matching rows by variant ID.
# Every dump row must have one value per variant.  Set lo=1 when the
# statistic is a p-value (smaller is more extreme).
#
# The dump has 6 significant digits, so values this close to the original
# statistic are "near ties": each may count as more extreme, tied, or less
# extreme, so the reported count must fall in a range.  The report has 4
# significant digits, so the range is rounded to match.
#
# As in PLINK, a variant that was skipped (all permutation values NA) counts
# as tied in every permutation, and a variant whose original statistic is NA
# is either reported as NA or as tied everywhere (EMP1 = half the
# permutations, EMP2 = all of them).  A single NA permutation value, from a
# failed regression, counts as not extreme.
function tol(o) {
  return 1e-5 * ((o < 0)? -o : o) + 1e-8
}
function beyond(v, o) {
  return lo? (v < o - tol(o)) : (v > o + tol(o))
}
function near(v, o) {
  return ((v - o <= tol(o)) && (o - v <= tol(o)))
}
function r4(x) {
  return sprintf("%.4g", x) + 0
}
function in_range(r, a, b) {
  return (r != "NA") && (r + 0 >= r4(a)) && (r + 0 <= r4(b))
}
function fail(msg) {
  if (++shown <= 5) {
    print FILENAME ": " msg
  }
  bad = 1
}
FILENAME == ARGV[1] {
  id[++nv] = $1
  col[$1] = nv
  next
}
FILENAME == ARGV[2] {
  if (NF != nv + 1) {
    fail("line " FNR " has " (NF - 1) " values; expected " nv)
  }
  if (FNR == 1) {
    for (k = 1; k <= nv; ++k) {
      orig[k] = $(k + 1)
    }
    next
  }
  ++np
  rb = "NA"
  for (k = 1; k <= nv; ++k) {
    v = $(k + 1)
    if (v == "NA") {
      ++na_ct[k]
      continue
    }
    if (rb == "NA" || (lo? (v + 0 < rb + 0) : (v + 0 > rb + 0))) {
      rb = v
    }
    if (orig[k] == "NA") {
      continue
    }
    if (near(v, orig[k])) {
      e1_hi[k] += 1
    } else if (beyond(v, orig[k])) {
      e1_lo[k] += 1
      e1_hi[k] += 1
    }
  }
  row_best[np] = rb
  next
}
FNR == 1 {
  next
}
{
  if (!($2 in col)) {
    fail("unexpected variant " $2)
    next
  }
  k = col[$2]
  ++nrep
  if (orig[k] == "NA") {
    if (!(($3 == "NA" && $4 == "NA") || (in_range($3, np / 2, np / 2) && in_range($4, np, np)))) {
      fail($2 " has no original statistic, but EMP1 " $3 " and EMP2 " $4)
    }
    next
  }
  if (na_ct[k] == np) {
    e1_lo[k] = np / 2
    e1_hi[k] = np / 2
  }
  if (!in_range($3, e1_lo[k], e1_hi[k])) {
    fail($2 " EMP1 " $3 "; the dump gives " (e1_lo[k] + 0) "-" (e1_hi[k] + 0))
  }
  e2_lo = 0
  e2_hi = 0
  for (p = 1; p <= np; ++p) {
    if (row_best[p] == "NA") {
      continue
    }
    if (near(row_best[p], orig[k])) {
      ++e2_hi
    } else if (beyond(row_best[p], orig[k])) {
      ++e2_lo
      ++e2_hi
    }
  }
  if (!in_range($4, e2_lo, e2_hi)) {
    fail($2 " EMP2 " $4 "; the dump gives " e2_lo "-" e2_hi)
  }
}
END {
  if (!nrep) {
    print "no usable report rows"
    bad = 1
  }
  exit bad
}
