# Fixture for --flip-scan dprime.  make_vcf.awk builds every block from two
# haplotypes, so nearly every pair inside a block has |D'| = 1 and D' is not
# really being estimated.  Here each block draws from a pool of four to seven
# haplotypes (each one a few mutations away from an earlier one), and each
# allele is occasionally redrawn at random, so the within-block pairs span the
# whole range of D' and the double heterozygotes actually need phasing.
BEGIN {
    OFS = "\t";
    sample_ct = 400;
    block_ct = 16;
    per_block = 12;
    # Every 9th variant has its alleles swapped in the cases (the second half
    # of the samples).
    flip_every = 9;
    seed = 20260922;

    print "##fileformat=VCFv4.2";
    print "##contig=<ID=1>";
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">";
    line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (s = 0; s < sample_ct; ++s) {
        line = line "\tS" s;
    }
    print line;

    for (b = 0; b < block_ct; ++b) {
        pool_ct = 4 + int(rnd() * 4);
        for (j = 0; j < per_block; ++j) {
            hap[0, j] = 0;
        }
        for (h = 1; h < pool_ct; ++h) {
            parent = int(rnd() * h);
            for (j = 0; j < per_block; ++j) {
                hap[h, j] = hap[parent, j];
                if (rnd() < 0.35) {
                    hap[h, j] = 1 - hap[h, j];
                }
            }
        }
        # Unequal pool frequencies: haplotype h is drawn with weight h + 1.
        wtot = pool_ct * (pool_ct + 1) / 2;
        for (s = 0; s < sample_ct; ++s) {
            draw1[s] = pick(pool_ct, wtot);
            draw2[s] = pick(pool_ct, wtot);
        }
        for (j = 0; j < per_block; ++j) {
            vidx = b * per_block + j;
            line = "1" OFS (b * 100000 + j * 2000 + 1) OFS "b" b "v" j OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS "GT";
            for (s = 0; s < sample_ct; ++s) {
                if (rnd() < 0.02) {
                    gt = "./.";
                } else {
                    a1 = hap[draw1[s], j];
                    a2 = hap[draw2[s], j];
                    if (rnd() < 0.04) {
                        a1 = (rnd() < 0.5)? 0 : 1;
                    }
                    if (rnd() < 0.04) {
                        a2 = (rnd() < 0.5)? 0 : 1;
                    }
                    if (((vidx % flip_every) == 0) && (s >= sample_ct / 2)) {
                        a1 = 1 - a1;
                        a2 = 1 - a2;
                    }
                    gt = a1 "/" a2;
                }
                line = line OFS gt;
            }
            print line;
        }
    }
}

function pick(ct, wtot,    r, h, acc) {
    r = rnd() * wtot;
    acc = 0;
    for (h = 0; h < ct; ++h) {
        acc += h + 1;
        if (r < acc) {
            return h;
        }
    }
    return ct - 1;
}

function rnd() {
    seed = (seed * 16807) % 2147483647;
    return seed / 2147483647;
}
