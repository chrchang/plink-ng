# Checks a plink2 --pgen-diff 'dosage[=tol]' report (default columns: ID, IID,
# DS1, DS2) against differences computed from DS-force VCFs of the two
# filesets (exported with CHROM relabeled to an autosome, so DS is on the
# 0..2 scale).
# Usage: awk -v tol=<tolerance> -v include_missing=<0|1> [-v x_males=<file>]
#   -f pdiff_compare.awk <.vcf 1> <.vcf 2> <.pdiff>
# Alleles are matched by string; the second fileset's alleles must be a
# subset of the first's, so the report's DS columns follow the first
# fileset's ALT order.  A sample differs iff some allele's dosage (REF
# included) differs by more than the tolerance (floor(floor(tol * 32768) / 2)
# in 16384ths of a copy; floor(tol * 32768) for the samples listed in
# x_males, whose DS values are reported on the 0..1 scale).  With
# include-missing, a sample missing in exactly one fileset differs too.
function absf(x) { return (x < 0)? -x : x }
BEGIN {
  FS = "\t";
  hap_tol = int(tol * 32768 * (1 + 1e-6));
  dip_tol = int(hap_tol / 2);
  if (x_males != "") {
    while ((getline line < x_males) > 0) is_male[line] = 1;
  }
}
FNR == 1 {++file_idx}
(file_idx == 1) || (file_idx == 2) {
  if (/^##/) next;
  if (/^#CHROM/) {
    for (i = 10; i <= NF; ++i) sid[file_idx, i] = $i;
    ns[file_idx] = NF;
    next;
  }
  v = $3;
  if (file_idx == 1) { vorder[++variant_ct] = v; alts1[v] = $5 }
  else alts2[v] = $5;
  na = split($5, a, ",");
  for (i = 10; i <= ns[file_idx]; ++i) {
    split($i, f, ":");
    s = sid[file_idx, i];
    if (f[2] == ".") { miss[file_idx, v, s] = 1; continue }
    split(f[2], d, ",");
    tot = 0;
    for (j = 1; j <= na; ++j) {
      u = int(d[j] * 16384 + 0.5);
      dos[file_idx, v, s, a[j]] = u;
      tot += u;
    }
    dos[file_idx, v, s, "\001REF"] = 32768 - tot;
  }
  next;
}
file_idx == 3 {
  if (/^#/) next;
  ++line_ct;
  key = $1 SUBSEP $2;
  if (key in reported) { if (bad++ < 10) print "duplicate report line: " $1 " " $2 > "/dev/stderr" }
  reported[key] = $3 "\t" $4;
  next;
}
END {
  for (k = 1; k <= variant_ct; ++k) {
    v = vorder[k];
    na = split(alts1[v], a, ",");
    a[0] = "\001REF";
    for (i = 10; i <= ns[1]; ++i) {
      s = sid[1, i];
      m1 = ((1, v, s) in miss);
      m2 = ((2, v, s) in miss);
      diff = 0;
      if (m1 || m2) {
        diff = include_missing && (!(m1 && m2));
      } else {
        cur_tol = (s in is_male)? hap_tol : dip_tol;
        for (j = 0; j <= na; ++j) {
          if (absf(dos[1, v, s, a[j]] - dos[2, v, s, a[j]]) > cur_tol) diff = 1;
        }
      }
      key = v SUBSEP s;
      if (!diff) {
        if (key in reported) { if (bad++ < 10) print "unexpected report line: " v " " s > "/dev/stderr" }
        continue;
      }
      ++expected_ct;
      if (!(key in reported)) { if (bad++ < 10) print "missing report line: " v " " s > "/dev/stderr"; continue }
      split(reported[key], ds_cols, "\t");
      for (fi = 1; fi <= 2; ++fi) {
        if (((fi == 1)? m1 : m2)) {
          if (ds_cols[fi] != ".") { if (bad++ < 10) print "DS mismatch: " v " " s " " ds_cols[fi] " vs. ." > "/dev/stderr" }
          continue;
        }
        if (split(ds_cols[fi], g, ",") != na) { if (bad++ < 10) print "DS count mismatch: " v " " s > "/dev/stderr"; continue }
        for (j = 1; j <= na; ++j) {
          e = dos[fi, v, s, a[j]] / ((s in is_male)? 32768 : 16384);
          if (absf(e - g[j]) > 6.2e-5) { if (bad++ < 10) print "DS mismatch: " v " " s " " ds_cols[fi] > "/dev/stderr"; break }
        }
      }
    }
  }
  print (line_ct + 0) " report lines, " (expected_ct + 0) " expected, " (bad + 0) " mismatches";
  exit (bad || (line_ct != expected_ct));
}
