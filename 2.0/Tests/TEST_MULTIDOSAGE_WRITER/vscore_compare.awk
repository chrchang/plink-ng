# Checks a plink2 --vscore ... cols=+altfreq,+nmiss report against variant
# scores computed from total ALT dosages.
# Usage: awk -f vscore_compare.awk <weight file> <DS-force .vcf> <.vscore>
# The weight file has a header line, then <IID> <W1> <W2> lines; each variant
# score is the sum over samples of <weight> * <total ALT dosage>, with missing
# dosages imputed as 2 * <ALT_FREQ column>.  Only the .vcf's variants are
# checked.
function absf(x) { return (x < 0)? -x : x }
BEGIN {FS = "\t"}
FNR == 1 {++file_idx}
file_idx == 1 {
  if (FNR > 1) {
    w1[$1] = $2;
    w2[$1] = $3;
  }
  next;
}
file_idx == 2 {
  if (/^##/) next;
  if (/^#CHROM/) {
    for (i = 10; i <= NF; ++i) id[i] = $i;
    ns = NF;
    next;
  }
  in_vcf[$3] = 1;
  for (i = 10; i <= ns; ++i) {
    split($i, f, ":");
    if (f[2] == ".") {
      miss[$3] = miss[$3] " " id[i];
      continue;
    }
    na = split(f[2], d, ",");
    x = 0;
    for (j = 1; j <= na; ++j) x += d[j];
    s1[$3] += w1[id[i]] * x;
    s2[$3] += w2[id[i]] * x;
  }
  next;
}
FNR == 1 {
  for (i = 1; i <= NF; ++i) {
    h = $i;
    sub(/^#/, "", h);
    c[h] = i;
  }
  next;
}
{
  v = $c["ID"];
  if (!(v in in_vcf)) next;
  af = $c["ALT_FREQ"];
  e1 = s1[v];
  e2 = s2[v];
  nm = split(miss[v], ms, " ");
  for (k = 1; k <= nm; ++k) {
    e1 += w1[ms[k]] * 2 * af;
    e2 += w2[ms[k]] * 2 * af;
  }
  ++n;
  if (($c["MISSING_CT"] != nm) || (absf($c["W1"] - e1) > 2e-3 + 1e-4 * absf(e1)) || (absf($c["W2"] - e2) > 2e-3 + 1e-4 * absf(e2))) {
    if (bad++ < 5) print "mismatch: " v ": " $c["MISSING_CT"] " " $c["W1"] " " $c["W2"] " vs. " nm " " e1 " " e2 > "/dev/stderr";
  }
}
END {
  print n " variant scores, " (bad + 0) " mismatches";
  exit (bad || (!n));
}
