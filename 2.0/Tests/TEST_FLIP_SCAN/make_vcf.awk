# Deterministic fileset with planted haplotype blocks: each block of 12
# variants is built from two haplotypes, so the pairs inside it are in strong
# LD and the pairs across blocks are not.
BEGIN {
    OFS = "\t";
    sample_ct = 400;
    block_ct = 20;
    per_block = 12;
    # Every 7th variant has its alleles swapped in the second half of the
    # samples, which are the cases: a strand inconsistency for --flip-scan to
    # find.
    flip_every = 7;
    seed = 20250905;

    print "##fileformat=VCFv4.2";
    print "##contig=<ID=1>";
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">";
    line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (s = 0; s < sample_ct; ++s) {
        line = line "\tS" s;
    }
    print line;

    for (b = 0; b < block_ct; ++b) {
        for (h = 0; h < 2; ++h) {
            for (j = 0; j < per_block; ++j) {
                hap[h, j] = (rnd() < 0.5)? 0 : 1;
            }
        }
        for (s = 0; s < sample_ct; ++s) {
            draw1[s] = (rnd() < 0.5)? 0 : 1;
            draw2[s] = (rnd() < 0.5)? 0 : 1;
        }
        for (j = 0; j < per_block; ++j) {
            line = "1" OFS (b * 100000 + j * 2000 + 1) OFS "b" b "v" j OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS "GT";
            for (s = 0; s < sample_ct; ++s) {
                if (rnd() < 0.02) {
                    gt = "./.";
                } else {
                    a1 = hap[draw1[s], j];
                    a2 = hap[draw2[s], j];
                    if ((flip_every > 0) && (((b * per_block + j) % flip_every) == 0) && (s >= sample_ct / 2)) {
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

function rnd() {
    seed = (seed * 16807) % 2147483647;
    return seed / 2147483647;
}
