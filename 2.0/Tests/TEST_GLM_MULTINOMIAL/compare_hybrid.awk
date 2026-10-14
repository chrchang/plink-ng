# Checks a 'firth-fallback' report against 'no-firth' and 'firth' reports of
# the same run.
# Usage: awk [-v min_firth_ct=1] -f compare_hybrid.awk <no-firth> <firth> <hybrid>
# * A FIRTH?=N row must be the no-firth row, apart from the FIRTH? column.
# * A FIRTH?=Y row must be the firth row, apart from that column, and the
#   no-firth row must have had NA statistics.
# * At least min_firth_ct rows (default 1) must be FIRTH?=Y, and every row
#   must be present in all three reports.

function fail(msg) {
  print FILENAME ": " $0
  print "  " msg
  bad = 1
  exit 1
}

BEGIN {
  FS = "\t"
  if (min_firth_ct == "") { min_firth_ct = 1 }
}

FNR == 1 {
  ++file_idx
  for (i = 1; i <= NF; ++i) { col[file_idx, $i] = i }
  next
}

{
  key = $col[file_idx, "ID"] "\t" $col[file_idx, "A1"]
}

file_idx == 1 {
  no_firth[key] = $0
  no_firth_chisq[key] = $col[1, "CHISQ"]
  ++no_firth_ct
  next
}

file_idx == 2 {
  firth[key] = $0
  next
}

{
  if (!(key in no_firth) || !(key in firth)) { fail("row missing from the no-firth or firth report") }
  firth_col = col[3, "FIRTH?"]
  row = ""
  for (i = 1; i <= NF; ++i) {
    if (i != firth_col) { row = row ((row == "")? "" : "\t") $i }
  }
  if ($firth_col == "N") {
    if (row != no_firth[key]) { fail("FIRTH?=N row differs from the no-firth row") }
  } else if ($firth_col == "Y") {
    if (row != firth[key]) { fail("FIRTH?=Y row differs from the firth row") }
    if (no_firth_chisq[key] != "NA") { fail("FIRTH?=Y, but the no-firth fit succeeded") }
    ++firth_ct
  } else {
    fail("FIRTH? is neither Y nor N")
  }
  ++row_ct
}

END {
  if (bad) { exit 1 }
  if (row_ct != no_firth_ct) {
    print "hybrid report has " row_ct " rows, no-firth report " no_firth_ct
    exit 1
  }
  if (firth_ct < min_firth_ct) {
    print "only " firth_ct " FIRTH?=Y rows"
    exit 1
  }
}
