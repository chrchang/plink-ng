# Reimplements PLINK 1.9's greedy --rel-cutoff pruning on a lower-triangular
# relationship matrix (one value per line, row by row: (2,1), (3,1), (3,2),
# ...), and prints the 0-based indices of the samples it keeps.
#   While some sample has exactly one remaining too-close relation, take the
#   first such sample and prune its partner.  Otherwise prune the first sample
#   with the most remaining too-close relations.  Stop when none remain.
#
# Usage: awk -v n=<sample count> -v cutoff=<x> -f rel_cutoff.awk values.txt

{
    v[NR - 1] = $1 + 0
}

END {
    k = 0
    for (i = 1; i < n; ++i) {
        for (j = 0; j < i; ++j) {
            if (v[k] > cutoff) {
                close_[i, j] = close_[j, i] = 1
                ++ct[i]
                ++ct[j]
            }
            ++k
        }
    }
    while (1) {
        pick = -1
        for (i = 0; i < n; ++i) {
            if (!(i in gone) && ct[i] == 1) {
                for (j = 0; j < n; ++j) {
                    if (j != i && ((i, j) in close_)) {
                        pick = j
                        break
                    }
                }
                break
            }
        }
        if (pick == -1) {
            best = 0
            for (i = 0; i < n; ++i) {
                if (!(i in gone) && ct[i] > best) {
                    best = ct[i]
                    pick = i
                }
            }
            if (pick == -1) break
        }
        gone[pick] = 1
        for (j = 0; j < n; ++j) {
            if ((pick, j) in close_) {
                delete close_[pick, j]
                delete close_[j, pick]
                --ct[j]
            }
        }
        ct[pick] = 0
    }
    for (i = 0; i < n; ++i) if (!(i in gone)) print i
}
