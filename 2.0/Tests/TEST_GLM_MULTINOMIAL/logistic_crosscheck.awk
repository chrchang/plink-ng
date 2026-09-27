# With two levels, multinomial logistic regression is ordinary logistic
# regression.  Checks a --glm multinomial=wald run on a 2-level phenotype
# (non-reference level 'hi') against a --glm no-firth hide-covar run on the
# same split as a case/control phenotype, both with cols=+beta (and
# multinomial-ref=lo, so that hi plays the role of the cases).
# Usage: awk -v terms=<comma-separated TEST names> [-v tol=1e-4] -f logistic_crosscheck.awk <.glm.logistic> <.glm.multinomial>
# terms is the logistic TEST name of each tested column: ADD, DOM, REC, or
# HET, or ADD,DOMDEV for 'genotypic' and HOM,HET for 'hethom'.
# * BETA_hi and SE_hi must match the logistic BETA and SE of each tested
#   column; in the additive model, a multiallelic variant's A1 list matches
#   one logistic row per allele (one joint fit on both sides).
# * The Wald statistic must match Z^2 for a single column, and the logistic
#   joint test's statistic times 2 (it is reported as chi-square / df) for
#   two.
# tol is looser than elsewhere because logistic regression stops earlier.
# Rows that either side does not report are skipped; at least min_ct rows
# (default 20) must be compared.

function close_enough(a, b) {
  d = a - b
  if (d < 0) { d = -d }
  if (b < 0) { b = -b }
  return d <= tol * b + 1e-6
}

function fail(msg) {
  print FILENAME ": " $0
  print "  " msg
  bad = 1
  exit 1
}

BEGIN {
  FS = "\t"
  if (tol == "") { tol = 1e-4 }
  if (min_ct == "") { min_ct = 20 }
  term_ct = split(terms, term_list, ",")
}

FNR == 1 {
  ++file_idx
  for (i = 1; i <= NF; ++i) { col[file_idx, $i] = i }
  stat_name = ((file_idx, "Z_STAT") in col)? "Z_STAT" : "Z_OR_F_STAT"
  if (file_idx == 1) { logistic_stat_col = col[1, stat_name] }
  next
}

file_idx == 1 {
  if ($col[1, "ERRCODE"] != ".") { next }
  key = $col[1, "ID"] "\t" $col[1, "A1"] "\t" $col[1, "TEST"]
  beta[key] = $col[1, "BETA"]
  se[key] = $col[1, "SE"]
  stat[key] = $logistic_stat_col
  next
}

$col[2, "ERRCODE"] == "." {
  id = $col[2, "ID"]
  a1_ct = split($col[2, "A1"], a1_list, ",")
  split($col[2, "BETA_hi"], mbeta, ",")
  split($col[2, "SE_hi"], mse, ",")
  entry_ct = 0
  for (a = 1; a <= a1_ct; ++a) {
    for (t = 1; t <= term_ct; ++t) {
      keys[++entry_ct] = id "\t" a1_list[a] "\t" term_list[t]
      if (!(keys[entry_ct] in beta)) { next }
    }
  }
  for (e = 1; e <= entry_ct; ++e) {
    if (!close_enough(mbeta[e], beta[keys[e]])) { fail("BETA " mbeta[e] " vs logistic " beta[keys[e]]) }
    if (!close_enough(mse[e], se[keys[e]])) { fail("SE " mse[e] " vs logistic " se[keys[e]]) }
  }
  chisq = $col[2, "CHISQ"]
  if (entry_ct == 1) {
    z = stat[keys[1]]
    if (!close_enough(chisq, z * z)) { fail("CHISQ vs logistic Z^2 " z * z) }
  } else if (term_ct == 2) {
    joint_key = id "\t" a1_list[1] "\tGENO_2DF"
    if (!(joint_key in stat)) { next }
    if (!close_enough(chisq, 2 * stat[joint_key])) { fail("CHISQ vs logistic joint statistic " 2 * stat[joint_key]) }
  }
  ++compared_ct
}

END {
  if (bad) { exit 1 }
  if (compared_ct < min_ct) {
    print "only " compared_ct " rows compared"
    exit 1
  }
}
