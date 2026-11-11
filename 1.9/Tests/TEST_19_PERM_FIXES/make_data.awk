# Writes tmp_<out>.ped and tmp_<out>.map.
#   n       sample count; sexes alternate (male first), and so do case/control
#           pairs
#   chrs    comma-separated chromosome codes, one variant per entry
#   seed    random seed
#   qt      if 1, the phenotype is quantitative
#   sparse  index of a variant with only one nonmissing call (0 for none)
#   mdiff   if 1, every third variant has a higher missing rate in cases
# Allele frequencies cycle through 0.01, 0.3 and 0.1, and every fifth variant
# (starting with the third) is monomorphic.  Female chrY calls are missing,
# and male chrY calls are homozygous.  Random numbers come from a 32-bit LCG,
# which awk's doubles evaluate exactly.
function rnd() {
  s = (s * 69069 + 1) % 4294967296
  return s / 4294967296
}
BEGIN {
  s = seed + 0
  nv = split(chrs, chr, ",")
  for (j = 1; j <= nv; ++j) {
    f[j] = (j % 5 == 3)? 0 : ((j % 3 == 1)? 0.01 : ((j % 3 == 2)? 0.3 : 0.1))
    print chr[j] " v" j " 0 " (1000 * j) > ("tmp_" out ".map")
  }
  for (i = 0; i < n; ++i) {
    sex = 1 + (i % 2)
    cc = 1 + (int(i / 2) % 2)
    pheno = qt? sprintf("%.3f", rnd() * 4) : cc
    line = "f" i " i" i " 0 0 " sex " " pheno
    for (j = 1; j <= nv; ++j) {
      g1 = (rnd() < f[j])? "A" : "C"
      g2 = (chr[j] == "Y")? g1 : ((rnd() < f[j])? "A" : "C")
      m = 0.03
      if (mdiff && (j % 3 == 0)) {
        m += 0.1 * (cc - 1)
      }
      if ((rnd() < m) || (j == sparse && i > 0) || (chr[j] == "Y" && sex == 2)) {
        g1 = "0"
        g2 = "0"
      }
      line = line " " g1 " " g2
    }
    print line > ("tmp_" out ".ped")
  }
}
