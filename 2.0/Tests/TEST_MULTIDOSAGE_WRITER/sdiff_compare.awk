# Checks plink2 --sample-diff 'dosage[=tol]' output against values computed
# from a DS-force VCF of the same (all-diploid) variants.
# Usage: awk -v tol=<tolerance> -v include_missing=<0|1> [-v haploid=1]
#   -f sdiff_compare.awk <DS-force .vcf> <.sdiff.summary> [<main .sdiff report>]
# With haploid=1 (e.g. pairs of chrX males, with the .vcf exported from an
# autosome-relabeled copy so that DS stays on the 0..2 scale), dosages are
# compared against floor(tol * 32768), IBS isn't counted, and reported DS
# values are on the 0..1 scale.
# Each sample's allele dosages (16384 = 1 copy) come from DS (REF = 2 copies
# minus the ALT dosages); a sample is missing iff DS is '.'.  A pair differs
# at a variant iff some allele's dosage differs by more than the tolerance
# (floor(floor(tol * 32768) / 2) in these units, as plink2 computes it).
# With include-missing, variants where one or both samples are missing count
# toward OBS_CT, and half-missing ones toward DIFF_CT (and are reported).
# Checked: OBS_CT, HALFMISS_CT (if present), DIFF_CT, IBS0_CT and IBS1_CT
# (from GT, over variants where both dosages are present and the hardcalls
# differ), and, if given, the main report's (variant, pair) lines and their
# DS values.
function absf(x) { return (x < 0)? -x : x }
# unphased allele pair, smaller code first ("." if missing)
function parse_gt(gt, out,    parts, n, t) {
  n = split(gt, parts, /[\/|]/);
  out[1] = parts[1];
  out[2] = (n == 2)? parts[2] : parts[1];
  if ((out[1] != ".") && (out[1] + 0 > out[2] + 0)) {
    t = out[1]; out[1] = out[2]; out[2] = t;
  }
}
BEGIN {
  FS = "\t";
  hap_tol = int(tol * 32768 * (1 + 1e-6));
  cur_tol = haploid? hap_tol : int(hap_tol / 2);
}
FNR == 1 {++file_idx}
file_idx == 1 {
  if (/^##/) next;
  if (/^#CHROM/) {
    for (i = 10; i <= NF; ++i) col_of[$i] = i;
    next;
  }
  ++variant_ct;
  vid[variant_ct] = $3;
  for (i = 10; i <= NF; ++i) {
    split($i, f, ":");
    gt[variant_ct, i] = f[1];
    ds[variant_ct, i] = f[2];
  }
  next;
}
file_idx == 2 {
  if (FNR == 1) {
    for (i = 1; i <= NF; ++i) {
      h = $i;
      sub(/^#/, "", h);
      sc[h] = i;
    }
    next;
  }
  ++pair_ct;
  id1 = $sc["IID1"];
  id2 = $sc["IID2"];
  c1 = col_of[id1];
  c2 = col_of[id2];
  pair_key[id1, id2] = 1;
  obs = 0; halfmiss = 0; diff = 0; ibs0 = 0; ibs1 = 0;
  for (v = 1; v <= variant_ct; ++v) {
    m1 = (ds[v, c1] == ".");
    m2 = (ds[v, c2] == ".");
    if (m1 && m2) {
      obs += include_missing;
      continue;
    }
    if (m1 || m2) {
      ++halfmiss;
      if (include_missing) {
        ++obs;
        # (half-missing variants count as differences here)
        ++diff;
        expected_line[vid[v], id1, id2] = 1;
      }
      continue;
    }
    ++obs;
    na = split(ds[v, c1], d1, ",");
    split(ds[v, c2], d2, ",");
    ref1 = 32768; ref2 = 32768; max_diff = 0;
    for (j = 1; j <= na; ++j) {
      u1 = int(d1[j] * 16384 + 0.5);
      u2 = int(d2[j] * 16384 + 0.5);
      ref1 -= u1;
      ref2 -= u2;
      if (absf(u1 - u2) > max_diff) max_diff = absf(u1 - u2);
    }
    if (absf(ref1 - ref2) > max_diff) max_diff = absf(ref1 - ref2);
    if (max_diff > cur_tol) {
      ++diff;
      expected_line[vid[v], id1, id2] = 1;
    }
    if (haploid) continue;
    parse_gt(gt[v, c1], g1);
    parse_gt(gt[v, c2], g2);
    if ((g1[1] == g2[1]) && (g1[2] == g2[2])) continue;
    if ((g1[1] == ".") || (g2[1] == ".")) continue;
    shared = (g1[1] == g2[1]) || (g1[1] == g2[2]) || (g1[2] == g2[1]) || (g1[2] == g2[2]);
    if (shared) ++ibs1; else ++ibs0;
  }
  # (HALFMISS_CT is only reported with include-missing)
  if (($sc["OBS_CT"] != obs) || (("HALFMISS_CT" in sc) && ($sc["HALFMISS_CT"] != halfmiss)) || ($sc["DIFF_CT"] != diff) || ($sc["IBS0_CT"] != ibs0) || ($sc["IBS1_CT"] != ibs1)) {
    print "summary mismatch for " id1 " " id2 ": " $sc["OBS_CT"] " " (("HALFMISS_CT" in sc)? $sc["HALFMISS_CT"] : "-") " " $sc["DIFF_CT"] " " $sc["IBS0_CT"] " " $sc["IBS1_CT"] " vs. " obs " " halfmiss " " diff " " ibs0 " " ibs1 > "/dev/stderr";
    ++bad;
  }
  next;
}
file_idx == 3 {
  if (/^#/) {
    for (i = 1; i <= NF; ++i) {
      h = $i;
      sub(/^#/, "", h);
      rc[h] = i;
    }
    next;
  }
  ++line_ct;
  key = $rc["ID"] SUBSEP $rc["IID1"] SUBSEP $rc["IID2"];
  if (!(key in expected_line)) {
    if (bad++ < 10) print "unexpected report line: " $rc["ID"] " " $rc["IID1"] " " $rc["IID2"] > "/dev/stderr";
    next;
  }
  seen[key] = 1;
  # DS values must match the VCF's (to printing precision)
  for (v = 1; v <= variant_ct; ++v) if (vid[v] == $rc["ID"]) break;
  ng = split($(NF - 1) "\t" $NF, gcols, "\t");
  for (k = 1; k <= 2; ++k) {
    vds = ds[v, col_of[$rc[(k == 1)? "IID1" : "IID2"]]];
    if (vds == "." || gcols[k] == ".") {
      if (vds != gcols[k]) { if (bad++ < 10) print "DS mismatch: " $rc["ID"] ": " gcols[k] " vs. " vds > "/dev/stderr" }
      continue;
    }
    na = split(vds, e, ",");
    if (split(gcols[k], g, ",") != na) { if (bad++ < 10) print "DS count mismatch: " $rc["ID"] > "/dev/stderr"; continue }
    for (j = 1; j <= na; ++j) if (absf((haploid? (e[j] / 2) : e[j]) - g[j]) > 6.2e-5) { if (bad++ < 10) print "DS mismatch: " $rc["ID"] ": " gcols[k] " vs. " vds > "/dev/stderr"; break }
  }
}
END {
  if (file_idx >= 3) {
    for (key in expected_line) {
      if (!(key in seen)) {
        split(key, kk, SUBSEP);
        if (bad++ < 10) print "missing report line: " kk[1] " " kk[2] " " kk[3] > "/dev/stderr";
      }
    }
  }
  print pair_ct " pairs, " (line_ct + 0) " report lines, " (bad + 0) " mismatches";
  exit (bad || (!pair_ct));
}
