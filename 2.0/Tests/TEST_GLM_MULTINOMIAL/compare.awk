# Compares a .glm.multinomial file against oracle.py output.
# Usage: awk -v stat=<LRT|SCORE|WALD> [-v tol=1e-5] -f compare.awk <ref> <plink2 output>
# * Rows are matched on (ID, A1); every reference row must be reported, with
#   the same OBS_CT and DF.
# * CHISQ, and any BETA_ / SE_ column plink2 reports (comma-separated lists
#   when several columns are tested), must agree within a relative tolerance
#   of tol (plink2 prints 6 significant digits).
# * A reference 'SEP' must come out as SEPARATION, 'CONST' as CONST_ALLELE,
#   'CORR' as CORR_TOO_HIGH, and 'VIF' as VIF_*, with CHISQ and P NA.

function close_enough(a, b) {
  d = a - b
  if (d < 0) { d = -d }
  if (b < 0) { b = -b }
  return d <= tol * b + 1e-9
}

function lists_agree(a, b,    na, nb, pa, pb, i) {
  na = split(a, pa, ",")
  nb = split(b, pb, ",")
  if (na != nb) { return 0 }
  for (i = 1; i <= na; ++i) {
    if ((pa[i] == "NA") || (pb[i] == "NA")) {
      if (pa[i] != pb[i]) { return 0 }
    } else if (!close_enough(pa[i], pb[i])) {
      return 0
    }
  }
  return 1
}

function fail(msg) {
  print FILENAME ": " $0
  print "  " msg
  bad = 1
  exit 1
}

BEGIN {
  FS = "\t"
  if (tol == "") { tol = 1e-5 }
  expected_errs["SEP"] = "^SEPARATION"
  expected_errs["CONST"] = "^CONST_ALLELE$"
  expected_errs["CORR"] = "^CORR_TOO_HIGH$"
  expected_errs["VIF"] = "^VIF_"
}

FNR == 1 {
  file_idx += 1
  for (i = 1; i <= NF; ++i) { col[file_idx, $i] = i }
  if (file_idx == 2) {
    ncoef = 0
    for (i = 1; i <= NF; ++i) {
      if (($i ~ /^BETA_/) || ($i ~ /^SE_/)) { coef[++ncoef] = $i }
    }
  }
  next
}

file_idx == 1 {
  key = $1 "\t" $2
  ref_line[key] = $0
  ++ref_ct
  next
}

{
  key = $col[2, "ID"] "\t" $col[2, "A1"]
  if (!(key in ref_line)) { fail("row missing from reference") }
  split(ref_line[key], r, "\t")
  ++seen_ct
  if (r[col[1, "OBS_CT"]] != $col[2, "OBS_CT"]) { fail("OBS_CT mismatch, expected " r[col[1, "OBS_CT"]]) }
  ref_stat = r[col[1, stat]]
  if (ref_stat in expected_errs) {
    expected_err = expected_errs[ref_stat]
    if (($col[2, "ERRCODE"] !~ expected_err) || ($col[2, "CHISQ"] != "NA") || ($col[2, "P"] != "NA")) {
      fail("expected " ref_stat " with NA statistics")
    }
  } else {
    if ($col[2, "ERRCODE"] != ".") { fail("unexpected ERRCODE") }
    if (r[col[1, "DF"]] != $col[2, "DF"]) { fail("DF mismatch, expected " r[col[1, "DF"]]) }
    if (!close_enough($col[2, "CHISQ"], ref_stat)) { fail("CHISQ mismatch, expected " ref_stat) }
  }
  for (c = 1; c <= ncoef; ++c) {
    name = coef[c]
    if (!lists_agree($col[2, name], r[col[1, name]])) { fail(name " mismatch, expected " r[col[1, name]]) }
  }
}

END {
  if (bad) { exit 1 }
  if (seen_ct != ref_ct) {
    print "reported " seen_ct " rows, reference has " ref_ct
    exit 1
  }
}
