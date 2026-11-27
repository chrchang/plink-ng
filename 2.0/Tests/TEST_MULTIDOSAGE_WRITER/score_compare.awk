# Checks a plink2 --score ... cols=+scoresums .sscore file against score sums
# computed from per-allele dosages.
# Usage: awk -v mod=<"" | dominant | recessive> -v imp=<"" | no-mean-imputation>
#   -v haploid=<0 | 1> -f score_compare.awk <.afreq> <score file> <DS-force .vcf>
#   <.sscore>
# The score file has <variant ID> <allele> <weight> lines; the .vcf (with every
# variant on a diploid chromosome, so DS is on the 0..2 scale) provides each
# sample's allele dosages.  'dominant' uses min(x, 1), 'recessive'
# max(x - 1, 0); without 'no-mean-imputation', missing dosages are imputed as
# 2f, 1 - (1-f)^2, or f^2 respectively.  With haploid=1, (additive) dosages
# and imputed values are halved.
function absf(x) { return (x < 0)? -x : x }
BEGIN {FS = "\t"}
FNR == 1 {++file_idx}
file_idx == 1 {
  if (FNR > 1) {
    nf = split($5, f, ",");
    for (j = 1; j <= nf; ++j) fr[$2, j] = f[j];
  }
  next;
}
file_idx == 2 {
  w[$1, $2] = $3;
  scored[$1] = 1;
  next;
}
file_idx == 3 {
  if (/^##/) next;
  if (/^#CHROM/) {
    for (i = 10; i <= NF; ++i) id[i] = $i;
    ns = NF;
    next;
  }
  if (!($3 in scored)) next;
  na = split($5, alts, ",");
  for (j = 1; j <= na; ++j) {
    if (!(($3, alts[j]) in w)) continue;
    ww = w[$3, alts[j]];
    p = fr[$3, j];
    for (i = 10; i <= ns; ++i) {
      split($i, f, ":");
      if (f[2] == ".") {
        if (imp == "") {
          if (mod == "dominant") x = 1 - (1 - p) ^ 2;
          else if (mod == "recessive") x = p ^ 2;
          else x = haploid? p : 2 * p;
          expected[id[i]] += ww * x;
        }
        continue;
      }
      split(f[2], d, ",");
      x = d[j];
      if (mod == "dominant") {
        if (x > 1) x = 1;
      } else if (mod == "recessive") {
        x = (x > 1)? x - 1 : 0;
      } else if (haploid) {
        x /= 2;
      }
      expected[id[i]] += ww * x;
    }
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
  ++n;
  got = $c["SCORE1_SUM"];
  e = expected[$c["IID"]] + 0;
  if (absf(got - e) > 1e-3 + 1e-4 * absf(e)) {
    if (bad++ < 5) print "mismatch: " $c["IID"] ": " got " vs. " e > "/dev/stderr";
  }
}
END {
  print n " sample scores, " (bad + 0) " mismatches";
  exit (bad || (!n));
}
