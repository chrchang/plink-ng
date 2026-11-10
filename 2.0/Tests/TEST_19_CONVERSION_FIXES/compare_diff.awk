# Normalizes a PLINK 1.9 --merge-mode 6/7 .diff report and a plink2
# --pgen-diff report to "ID FID IID NEW OLD" lines, with each genotype's
# alleles in sorted order.  Sort the output before comparing.
# Usage:
#   awk -v strip=<ID prefix to remove> -f pdiff.awk <1.9 .diff>
#   awk -f pdiff.awk <--pgen-diff fileset's .bim> <plink2 .pdiff>
# (GT1, from the main fileset, is OLD and GT2 is NEW; --pgen-diff codes both
# against the --pgen-diff fileset's REF/ALT, i.e. its .bim A2/A1.)
function norm(g,   p) { split(g, p, "/"); return (p[1] <= p[2])? p[1] "/" p[2] : p[2] "/" p[1] }
function letters(v, g,   p) { split(g, p, "/"); return norm(al[v, p[1]] "/" al[v, p[2]]) }
FILENAME ~ /\.bim$/ { al[$2, 0] = $6; al[$2, 1] = $5; next }
FILENAME ~ /\.pdiff$/ { if (/^#/) next; print $1, $2, $3, letters($1, $5), letters($1, $4); next }
FNR > 1 { id = $1; if (strip != "" && index(id, strip) == 1) id = substr(id, length(strip) + 1); print id, $2, $3, norm($4), norm($5) }
