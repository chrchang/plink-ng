# Writes tmp_ld.ped and tmp_ld.map: 300 samples (the first 150 male) and
# variants on chr1, chrX and chr2.  Haplotypes are copied from 6 founders, with
# a 4% switch rate per variant and 1% mutations, so nearby variants are in LD.
# Males are hemizygous on chrX, and 2% of calls are missing.
function rand01() { seed = (seed * 69069 + 1) % 4294967296; return seed / 4294967296 }
function randint(n) { return int(rand01() * n) }
BEGIN {
  seed = 2026
  n = 300; nmale = 150; nfounder = 6
  split("1 X 2", chroms, " ")
  split("40 40 20", nvars, " ")
  vtot = 0
  for (c = 1; c <= 3; ++c) {
    for (j = 0; j < nvars[c]; ++j) {
      p = 0.1 + 0.5 * rand01()
      for (f = 0; f < nfounder; ++f) founder[f, j] = (rand01() < p)
    }
    for (h = 0; h < 2 * n; ++h) {
      f = randint(nfounder)
      for (j = 0; j < nvars[c]; ++j) {
        if (rand01() < 0.04) f = randint(nfounder)
        hap[h, vtot + j] = founder[f, j] != (rand01() < 0.01)
      }
    }
    for (j = 0; j < nvars[c]; ++j) {
      print chroms[c] "\t" chroms[c] "_" j "\t0\t" (10000 + 1000 * j) > "tmp_ld.map"
      chrom_of[vtot + j] = chroms[c]
    }
    vtot += nvars[c]
  }
  for (s = 0; s < n; ++s) {
    line = "f" s " i" s " 0 0 " ((s < nmale)? 1 : 2) " -9"
    for (v = 0; v < vtot; ++v) {
      a1 = hap[2 * s, v]
      a2 = (chrom_of[v] == "X" && s < nmale)? a1 : hap[2 * s + 1, v]
      if (rand01() < 0.02) {
        line = line " 0 0"
      } else {
        line = line " " (a1? "G" : "A") " " (a2? "G" : "A")
      }
    }
    print line > "tmp_ld.ped"
  }
}
