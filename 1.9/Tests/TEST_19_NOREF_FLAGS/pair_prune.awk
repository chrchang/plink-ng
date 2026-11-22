# Recomputes PLINK 1.9 --indep-pairwise/--indep-pairphase pruning from a
# table of pairwise r^2 values, printing the IDs of the variants kept (the
# .prune.in contents) in .bim order.
#
# Inputs, in order:
#   1. --freq ... counts output (.frq.counts), for minor allele frequencies;
#   2. --r2 table output covering every same-chromosome pair (columns CHR_A
#      BP_A SNP_A CHR_B BP_B SNP_B R2 ...);
#   3. the .bim file.
# Variables: win (window size, in variants, or in kb with kb=1), step (step
# size, in variants) and thr (r^2 threshold).
#
# Within each window, PLINK repeatedly sweeps the unpruned variants in order:
# for each variant i, the first later unpruned variant j with r^2(i, j) > thr
# causes the lower-MAF one of the two to be pruned (j on ties), and the sweep
# moves on to the next i.  Sweeps repeat until one prunes nothing.  Windows
# start every step variants, and contain win variants (or with kb=1, the
# variants within win kb of the first).  Monomorphic variants are pruned
# immediately.
#
# Also sets the exit status to 2 if some r^2 value is within 2e-6 of thr,
# since rounding in the r^2 table could then change the result.

FILENAME == ARGV[1] {
    if (FNR > 1) {
        c1 = $5
        c2 = $6
        maf[$2] = (c1 + c2)? ((c1 < c2)? c1 : c2) / (c1 + c2) : 0
    }
    next
}

FILENAME == ARGV[2] {
    if ((FNR > 1) && ($7 + 0 == $7)) {
        r2[$3, $6] = $7
        d = $7 - thr
        if ((d < 2e-6) && (d > -2e-6)) close_call = 1
    }
    next
}

{
    if ($1 != cur_chr) {
        if (cur_chr != "") prune_chrom()
        cur_chr = $1
        n = 0
    }
    id[n] = $2
    pos[n] = $4
    ++n
}

function window_end(start,    e) {
    if (!kb) {
        e = start + win
        return (e > n)? n : e
    }
    for (e = start; (e < n) && (pos[e] <= pos[start] + 1000 * win); ++e) {}
    return e
}

function prune_chrom(    i, j, start, end, any, mi, mj) {
    for (i = 0; i < n; ++i) pruned[i] = (maf[id[i]] == 0)
    for (start = 0; start < n; start += step) {
        end = window_end(start)
        do {
            any = 0
            for (i = start; i < end; ++i) {
                if (pruned[i]) continue
                for (j = i + 1; j < end; ++j) {
                    if (pruned[j]) continue
                    if (((id[i], id[j]) in r2) && (r2[id[i], id[j]] > thr)) {
                        any = 1
                        mi = maf[id[i]]
                        mj = maf[id[j]]
                        # (PLINK's tolerance is negligible next to the
                        # difference between two distinct allele counts.)
                        if (mi < (1 - 1e-9) * mj) {
                            pruned[i] = 1
                        } else {
                            pruned[j] = 1
                        }
                        break
                    }
                }
            }
        } while (any)
    }
    for (i = 0; i < n; ++i) {
        if (!pruned[i]) print id[i]
    }
}

END {
    if (cur_chr != "") prune_chrom()
    if (close_call) exit 2
}
