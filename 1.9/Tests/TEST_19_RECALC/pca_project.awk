# Recomputes --pca --pca-cluster-names/--pca-clusters eigenvectors for every
# sample by projection.  The PCA is run on the cluster samples only, and each
# sample's PC k score is sum_j z_j w_jk over variants j, where w is the subset
# PCA's variant weights (--pca var-wts) and z_j = (A1 dosage - 2p) /
# sqrt(2p(1-p)) with p the A1 frequency from a --freqx report.  Scores are
# then scaled per PC to the subset PCA's eigenvectors by least squares.
# Prints "FID IID e1 ... eK" for each sample in the .raw file, in order.
#
# Usage: awk -v k=<PC count> -f pca_project.awk x.frqx sub.eigenvec.var
#            sub.eigenvec x.raw

FILENAME == ARGV[1] {
    if (FNR > 1) p[$2] = (2 * $5 + $6) / (2 * ($5 + $6 + $7))
    next
}
FILENAME == ARGV[2] {
    for (c = 1; c <= k; ++c) w[$2, c] = $(4 + c)
    next
}
FILENAME == ARGV[3] {
    for (c = 1; c <= k; ++c) sub_ev[$1 " " $2, c] = $(2 + c)
    next
}
FNR == 1 {
    for (i = 7; i <= NF; ++i) {
        name = $i
        sub(/_[^_]*$/, "", name)
        col[i] = name
    }
    next
}
{
    id = $1 " " $2
    order[++ns] = id
    for (c = 1; c <= k; ++c) sc[id, c] = 0
    for (i = 7; i <= NF; ++i) {
        q = p[col[i]]
        z = ($i - 2 * q) / sqrt(2 * q * (1 - q))
        for (c = 1; c <= k; ++c) sc[id, c] += z * w[col[i], c]
    }
}
END {
    for (c = 1; c <= k; ++c) {
        num = den = 0
        for (s = 1; s <= ns; ++s) {
            id = order[s]
            if ((id, c) in sub_ev) {
                num += sub_ev[id, c] * sc[id, c]
                den += sc[id, c] * sc[id, c]
            }
        }
        f[c] = num / den
    }
    for (s = 1; s <= ns; ++s) {
        id = order[s]
        line = id
        for (c = 1; c <= k; ++c) line = line " " sprintf("%.6g", f[c] * sc[id, c])
        print line
    }
}
