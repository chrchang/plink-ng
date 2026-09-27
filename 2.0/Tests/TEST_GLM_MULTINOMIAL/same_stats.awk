# Checks that two .glm.multinomial files report the same statistics.
# Usage: awk [-v negate=1] [-v tol=1e-5] -f same_stats.awk <file 1> <file 2>
# * Same variants in the same order, same DF, and CHISQ equal within a
#   relative tolerance of tol.
# * Same ERRCODE up to the allele named after a comma (SEPARATION,REF can
#   become SEPARATION,ALT1 when REF and ALT are swapped).
# * With negate=1, every BETA_ column of file 2 must be the negation of file
#   1's (within tol), as when the counted allele is flipped.

function fail(msg) {
  print msg
  failed = 1
  exit 1
}

function close_enough(a, b) {
  d = a - b
  if (d < 0) { d = -d }
  if (b < 0) { b = -b }
  return d <= tol * b + 1e-9
}

BEGIN {
  FS = "\t"
  if (tol == "") { tol = 1e-5 }
}

FNR == 1 {
  ++file_idx
  for (i = 1; i <= NF; ++i) { col[file_idx, $i] = i }
  next
}

file_idx == 1 {
  line1[FNR] = $0
  ct1 = FNR
  next
}

{
  split(line1[FNR], a, "\t")
  if (a[col[1, "ID"]] != $col[2, "ID"]) { fail("variant order differs at line " FNR) }
  e1 = a[col[1, "ERRCODE"]]
  e2 = $col[2, "ERRCODE"]
  sub(/,.*/, "", e1)
  sub(/,.*/, "", e2)
  if (e1 != e2) { fail($col[2, "ID"] ": ERRCODE " e1 " vs " e2) }
  c1 = a[col[1, "CHISQ"]]
  c2 = $col[2, "CHISQ"]
  if ((c1 == "NA") || (c2 == "NA")) {
    if (c1 != c2) { fail($col[2, "ID"] ": CHISQ " c1 " vs " c2) }
  } else if (!close_enough(c2, c1)) {
    fail($col[2, "ID"] ": CHISQ " c1 " vs " c2)
  }
  if (a[col[1, "DF"]] != $col[2, "DF"]) { fail($col[2, "ID"] ": DF differs") }
  if (negate) {
    for (name in col) {
      split(name, parts, SUBSEP)
      if ((parts[1] != 1) || (parts[2] !~ /^BETA_/)) { continue }
      b1 = a[col[1, parts[2]]]
      b2 = $col[2, parts[2]]
      if ((b1 == "NA") || (b2 == "NA")) {
        if (b1 != b2) { fail($col[2, "ID"] ": " parts[2] " " b1 " vs " b2) }
      } else if (!close_enough(-b2, b1)) {
        fail($col[2, "ID"] ": " parts[2] " " b1 " vs " b2 " (expected negation)")
      }
    }
  }
  ct2 = FNR
}

END {
  if (failed) { exit 1 }
  if (ct1 != ct2) { fail("line counts differ: " ct1 " vs " ct2) }
}
