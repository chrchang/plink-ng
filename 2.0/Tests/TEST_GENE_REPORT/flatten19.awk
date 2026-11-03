# Reduces plink 1.9's block-structured .range.report to the
# (gene, chromosome, gene start, gene end, variant ID) tuples that plink 2.0's
# --gene-report emits, in the same order.
/ -- /{
  hdr = $0
  gene = $1
  sub(/^.* -- /, "", hdr)
  sub(/ \(.*$/, "", hdr)
  split(hdr, a, ":")
  chrom = a[1]
  sub(/^chr/, "", chrom)
  ranges = a[2]
  gsub(/,/, " ", ranges)
  n = split(ranges, r, " ")
  split(r[1], s1, "\\.\\.")
  gene_start = s1[1]
  split(r[n], s2, "\\.\\.")
  gene_end = s2[2]
  next
}
/^ *DIST /{next}
NF == 0{next}
{ print gene "\t" chrom "\t" gene_start "\t" gene_end "\t" $4 }
