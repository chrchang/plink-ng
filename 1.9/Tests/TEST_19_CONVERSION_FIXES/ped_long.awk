# Converts a .map (first file) + .ped (second file) to one
# "FID IID variant position genotype" line per call, with each genotype's
# alleles in sorted order.  After sorting, filesets can be compared
# regardless of variant order and allele orientation.
NR == FNR { id[FNR] = $2; pos[FNR] = $4; nv = FNR; next }
{
  for (v = 1; v <= nv; ++v) {
    a = $(5 + 2 * v); b = $(6 + 2 * v)
    if (a > b) { t = a; a = b; b = t }
    print $1, $2, id[v], pos[v], a "/" b
  }
}
