# Writes tmp_df.ped/.map, tmp_df_within.txt and tmp_df_remove.txt.
#   60 nuclear families: both parents genotyped, one affected child, and
#     (when run with -v extra=1) one more, unaffected, child;
#   60 sibships without parents, 2-4 children with mixed status;
#   150 unrelated samples in 5 --within clusters.
# 40 variants with no missing calls; variant 0 raises risk.  The random
# number stream doesn't depend on 'extra', so both versions share all other
# genotypes.  The --remove list has every third unrelated sample and the
# first child of some sibships.
function rnd() { s = (s * 69069 + 1) % 4294967296; return s / 4294967296 }
function hap(k,   v) { for (v = 0; v < nv; ++v) h[k, v] = (rnd() < p[v]) }
function kid(dst, pa, ma,   v) {
  for (v = 0; v < nv; ++v) { h[dst, v] = h[pa + int(rnd() * 2), v]; h[dst + 1, v] = h[ma + int(rnd() * 2), v] }
}
function geno(k,   v, out) {
  out = ""
  for (v = 0; v < nv; ++v) out = out " " (h[k, v]? "A" : "C") " " (h[k + 1, v]? "A" : "C")
  return out
}
function aff(k) { return (rnd() < 0.25 + 0.25 * (h[k, 0] + h[k + 1, 0]))? 2 : 1 }
BEGIN {
  s = 31; nv = 40
  for (v = 0; v < nv; ++v) p[v] = 0.15 + 0.35 * rnd()
  for (v = 0; v < nv; ++v) print "1 snp" v " 0 " (1000 * (v + 1)) > "tmp_df.map"
  for (f = 0; f < 60; ++f) {
    hap(0); hap(1); hap(2); hap(3)
    print "T" f " P 0 0 1 1" geno(0) > "tmp_df.ped"
    print "T" f " M 0 0 2 1" geno(2) > "tmp_df.ped"
    kid(4, 0, 2); print "T" f " C0 P M 1 2" geno(4) > "tmp_df.ped"
    kid(6, 0, 2); if (extra) print "T" f " U1 P M 2 1" geno(6) > "tmp_df.ped"
  }
  for (f = 0; f < 60; ++f) {
    hap(0); hap(1); hap(2); hap(3)
    for (c = 0; c < 2 + f % 3; ++c) {
      kid(4, 0, 2)
      print "S" f " C" c " 0 0 " (1 + c % 2) " " aff(4) geno(4) > "tmp_df.ped"
    }
  }
  for (i = 0; i < 150; ++i) {
    hap(0); hap(1)
    print "U" i " U 0 0 " (1 + i % 2) " " aff(0) geno(0) > "tmp_df.ped"
    print "U" i, "U", "K" (i % 5) > "tmp_df_within.txt"
    if (i % 3 == 0) print "U" i, "U" > "tmp_df_remove.txt"
  }
  for (f = 0; f < 60; f += 11) print "S" f, "C0" > "tmp_df_remove.txt"
}
