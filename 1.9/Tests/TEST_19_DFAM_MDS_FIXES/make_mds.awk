# Writes tmp_mds.ped/.map: 12 subpopulations (3 regions x 4) of 6-18
# samples each, well enough separated that --cluster --K 12 recovers them,
# 400 variants.
function rnd() { s = (s * 69069 + 1) % 4294967296; return s / 4294967296 }
function clip(x) { return (x < 0.03)? 0.03 : ((x > 0.97)? 0.97 : x) }
# approximately normal (sum of 4 uniforms), sd ~ sd
function nrm(sd) { return (rnd() + rnd() + rnd() + rnd() - 2) * sd * 1.732 }
BEGIN {
  s = 11; nv = 400
  split("6 8 10 12 14 16 18 6 9 11 13 17", sizes, " ")
  for (v = 0; v < nv; ++v) {
    anc = 0.2 + 0.6 * rnd()
    for (r = 0; r < 3; ++r) reg[r, v] = clip(anc + nrm(0.12))
    print "1 m" v " 0 " (1000 * (v + 1)) > "tmp_mds.map"
  }
  for (p = 1; p <= 12; ++p) {
    for (v = 0; v < nv; ++v) pf[v] = clip(reg[int((p - 1) / 4), v] + nrm(0.2))
    for (i = 0; i < sizes[p]; ++i) {
      line = "q" p "_" i " q" p "_" i " 0 0 1 -9"
      for (v = 0; v < nv; ++v) {
        g = (rnd() < pf[v]) + (rnd() < pf[v])
        line = line ((g == 0)? " C C" : ((g == 1)? " A C" : " A A"))
      }
      print line > "tmp_mds.ped"
    }
  }
}
