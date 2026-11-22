# Writes a .ped and .map of n unrelated samples.  Chromosome 1 holds n_block
# LD blocks of block_len variants, 1 kb apart: each block has 4 ancestral
# haplotypes, and each sample haplotype copies one of them with 4% of its
# alleles flipped.  Chromosome 2 holds n_indep independent variants.  Calls are
# missing at rate miss.  The phenotype depends on causal variants 3 and
# 2 * block_len + 5 (chromosome 1) and the first chromosome 2 variant: with
# qt=1 it is their weighted allele count plus noise; otherwise those samples
# whose score exceeds a cutoff are cases, with 10% of statuses then flipped.
#
# Usage: awk -v n=<n> -v n_block=<b> -v block_len=<l> -v n_indep=<m>
#            -v seed=<s> -v miss=<rate> [-v qt=1] -v map=<.map file>
#            -f make_blocks.awk > x.ped

function rnd() {
    s = (s * 69069 + 1) % 4294967296
    return s / 4294967296
}

BEGIN {
    s = seed
    m1 = n_block * block_len
    m = m1 + n_indep
    for (b = 0; b < n_block; ++b) {
        for (h = 0; h < 4; ++h) {
            for (k = 1; k <= block_len; ++k) {
                anc[b, h, k] = (rnd() < 0.5)? "C" : "A"
            }
        }
    }
    for (j = m1 + 1; j <= m; ++j) freq[j] = 0.1 + 0.4 * rnd()
    c1 = 3
    c2 = 2 * block_len + 5
    c3 = m1 + 1
    for (j = 1; j <= m; ++j) {
        if (j <= m1) {
            printf "1 v%d 0 %d\n", j, 1000 * j > map
        } else {
            printf "2 v%d 0 %d\n", j, 1000 * j > map
        }
    }
    for (i = 1; i <= n; ++i) {
        for (a = 1; a <= 2; ++a) {
            for (b = 0; b < n_block; ++b) {
                h = int(rnd() * 4)
                for (k = 1; k <= block_len; ++k) {
                    al = anc[b, h, k]
                    if (rnd() < 0.04) al = (al == "C")? "A" : "C"
                    hap[a, b * block_len + k] = al
                }
            }
            for (j = m1 + 1; j <= m; ++j) {
                hap[a, j] = (rnd() < freq[j])? "C" : "A"
            }
        }
        score = (hap[1, c1] == "C") + (hap[2, c1] == "C") + 0.7 * ((hap[1, c2] == "C") + (hap[2, c2] == "C")) + 0.5 * ((hap[1, c3] == "C") + (hap[2, c3] == "C"))
        if (qt) {
            pheno = sprintf("%.4f", score + 2 * (rnd() + rnd() + rnd() - 1.5))
        } else {
            pheno = 1 + (score + rnd() * 2 > 2.6)
            if (rnd() < 0.1) pheno = 3 - pheno
        }
        line = "f" i " i" i " 0 0 " (1 + (i % 2)) " " pheno
        for (j = 1; j <= m; ++j) {
            if (rnd() < miss) {
                line = line " 0 0"
            } else {
                line = line " " hap[1, j] " " hap[2, j]
            }
        }
        print line
    }
}
