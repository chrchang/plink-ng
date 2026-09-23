# Compares a .glm.multinomial file against oracle.py output.
# Usage: awk -v stat=<LRT|SCORE|WALD> [-v tol=1e-5] -f compare.awk <ref> <plink2 output>
# * Every variant in the reference must be reported, with the same OBS_CT and
#   DF.
# * CHISQ, and any BETA_ / SE_ column plink2 reports, must agree within a
#   relative tolerance of tol (plink2 prints 6 significant digits).
# * A reference 'SEP' must come out as SEPARATION with CHISQ, P and the
#   coefficients all NA.

function close_enough(a, b) {
  d = a - b
  if (d < 0) { d = -d }
  if (b < 0) { b = -b }
  return d <= tol * b + 1e-9
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
  id = $1
  ref_line[id] = $0
  ++ref_ct
  next
}

{
  id = $col[2, "ID"]
  if (!(id in ref_line)) { fail("variant missing from reference") }
  split(ref_line[id], r, "\t")
  ++seen_ct
  if (r[col[1, "OBS_CT"]] != $col[2, "OBS_CT"]) { fail("OBS_CT mismatch, expected " r[col[1, "OBS_CT"]]) }
  ref_stat = r[col[1, stat]]
  if (ref_stat == "SEP") {
    if (($col[2, "ERRCODE"] !~ /^SEPARATION/) || ($col[2, "CHISQ"] != "NA") || ($col[2, "P"] != "NA")) {
      fail("expected SEPARATION with NA statistics")
    }
  } else {
    if ($col[2, "ERRCODE"] != ".") { fail("unexpected ERRCODE") }
    if (r[col[1, "DF"]] != $col[2, "DF"]) { fail("DF mismatch, expected " r[col[1, "DF"]]) }
    if (!close_enough($col[2, "CHISQ"], ref_stat)) { fail("CHISQ mismatch, expected " ref_stat) }
  }
  for (c = 1; c <= ncoef; ++c) {
    name = coef[c]
    ref_val = r[col[1, name]]
    val = $col[2, name]
    if ((ref_val == "NA") || (val == "NA")) {
      if (ref_val != val) { fail(name " mismatch, expected " ref_val) }
    } else if (!close_enough(val, ref_val)) {
      fail(name " mismatch, expected " ref_val)
    }
  }
}

END {
  if (bad) { exit 1 }
  if (seen_ct != ref_ct) {
    print "reported " seen_ct " variants, reference has " ref_ct
    exit 1
  }
}
