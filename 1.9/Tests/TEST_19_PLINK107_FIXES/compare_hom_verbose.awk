# Compares PLINK 1.9 --homozyg group-verbose reports with PLINK 1.07
# --homozyg-group --homozyg-verbose reports.  The two versions form pools
# differently (1.07 reports every overlapping subset), so files are not
# compared whole.  Instead:
#   * every genotype cell (sample, variant) printed by both must agree;
#   * the trailing consensus allele list must agree for each variant on which
#     all 1.07 files agree.
# Samples are keyed by FID, which must be unique: 1.07's "Group" member lists
# print the FID in place of the IID.
# Usage: awk -v n107=<number of 1.07 files> -v reversed=<file of variant IDs
#        whose .bim A1 is the major allele> -f compare_hom_verbose.awk
#        <1.07 files...> <1.9 files...>
function strip(g) { gsub(/[][]/, "", g); return g }
FNR == 1 { ++fileno; nmem = 0; is107 = (fileno <= n107) }
$1 == "Group" { nmem = 0; next }
# member lists: "  3)   f2   i2   1"
$1 ~ /^[0-9]+\)$/ { mem[++nmem] = $2; next }
NF == 0 || $1 == "SNP" || $1 == "FID" { next }
{
  if (NF == 2 && $2 !~ /\//) {
    # trailing consensus allele list
    if (is107) {
      if (!($1 in cons107)) cons107[$1] = $2
      else if (cons107[$1] != $2) conflict[$1] = 1
    } else {
      cons19[$1] = $2
    }
    next
  }
  # grid line, possibly with a consensus-allele column after the ID
  first = ($2 !~ /\//)? 3 : 2
  for (k = first; k <= NF; ++k) {
    key = mem[k - first + 1] " " $1
    g = strip($k)
    if (is107) {
      cell107[key] = g
    } else if (key in cell107) {
      ++ncell
      if (cell107[key] != g) { if (++nbad <= 5) print "cell differs: " key " " cell107[key] " vs " g; bad = 1 }
    }
  }
}
END {
  while ((getline v < reversed) > 0) rev[v] = 1
  for (v in cons19) {
    if (!(v in cons107) || (v in conflict)) continue
    ++ncons
    if (v in rev) ++nrev
    if (cons107[v] != cons19[v]) { if (++nbad <= 10) print "consensus differs: " v " " cons107[v] " vs " cons19[v]; bad = 1 }
  }
  if (ncell < 100 || ncons < 50 || nrev < 10) {
    print "too little overlap: " ncell " cells, " ncons " consensus alleles (" nrev " reversed)"
    bad = 1
  }
  if (!bad) print ncell " cells and " ncons " consensus alleles (" nrev " reversed) matched"
  exit bad
}
