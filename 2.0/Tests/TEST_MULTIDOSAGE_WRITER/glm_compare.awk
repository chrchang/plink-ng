# Compares multiallelic --glm output rows against reference rows.
# Usage: awk -f glm_compare.awk <reference rows> <glm output>...
# Reference rows are "<output suffix>\t<variant ID>\t<glm row>", where the
# glm row has the same columns as the corresponding output, whose suffix is
# everything after the first '.' in its filename (e.g. "Q.glm.linear").  Rows
# are matched on (suffix, ID, A1, occurrence index); every reference row of a
# suffix present among the outputs must be matched.  Numeric columns must
# agree to a relative 1e-4 (odds ratios on the log scale).
function absf(x) { return (x < 0)? -x : x }
function close_enough(x, y) {
  if ((x "") == (y "")) return 1;
  if ((x == "NA") || (y == "NA")) return 0;
  return absf(x - y) <= 1e-4 * (absf(x) + absf(y)) + 1e-6;
}
# log10 of a positive number, computed from its text, since odds ratios can
# exceed double range (e.g. 5e+828 under separation, which awk turns into inf;
# inf - inf is NaN, which awks compare inconsistently).
function text_log10(x,    parts, n) {
  n = split(tolower(x), parts, "e");
  return ((n == 2)? (parts[2] + 0) : 0) + log(parts[1] + 0) / log(10);
}
function or_close_enough(x, y) {
  if ((x "") == (y "")) return 1;
  if ((x == "NA") || (y == "NA") || (x + 0 <= 0) || (y + 0 <= 0)) return close_enough(x, y);
  # same relative tolerance as close_enough(), on the log scale
  return absf(text_log10(x) - text_log10(y)) <= 1e-4;
}
BEGIN {
  FS = "\t";
  name_ct = split("A1_FREQ OBS_CT BETA OR SE LOG(OR)_SE P CHISQ DF MIN_EXPECTED CASE_NON_A1_CT CASE_HET_A1_CT CASE_HOM_A1_CT CTRL_NON_A1_CT CTRL_HET_A1_CT CTRL_HOM_A1_CT", names, " ");
}
FNR == 1 {++file_idx}
file_idx == 1 {
  # (A1 is the 7th glm column)
  key = $1 SUBSEP $2 SUBSEP $9;
  key = key SUBSEP (++ref_occ[key]);
  row = $3;
  for (i = 4; i <= NF; ++i) row = row "\t" $i;
  ref[key] = row;
  ++ref_sfx_ct[$1];
  next;
}
/^#/ {
  sfx = FILENAME;
  sub(/^[^.]*\./, "", sfx);
  seen_sfx[sfx] = 1;
  for (i = 1; i <= NF; ++i) col[sfx, $i] = i;
  next;
}
{
  key = sfx SUBSEP $3 SUBSEP $col[sfx, "A1"];
  key = key SUBSEP (++occ[key]);
  if (!(key in ref)) next;
  ++cmp_ct;
  split(ref[key], r, "\t");
  for (n = 1; n <= name_ct; ++n) {
    c = col[sfx, names[n]];
    if (c && (!((names[n] == "OR")? or_close_enough($c, r[c]) : close_enough($c, r[c])))) {
      if (bad++ < 10) print "mismatch: " sfx " " $3 " " $col[sfx, "A1"] " " names[n] ": " $c " vs. " r[c] > "/dev/stderr";
    }
  }
}
END {
  for (s in seen_sfx) expected_ct += ref_sfx_ct[s];
  print cmp_ct " multiallelic glm rows compared, " (bad + 0) " mismatches";
  exit (bad || (cmp_ct != expected_ct) || (!cmp_ct));
}
