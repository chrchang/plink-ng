# Phased fileset for the known-phase checks.  Each block of 12 variants is
# built from four haplotypes, so the double-het phase EM has real ambiguity to
# resolve and known phase changes the answer.
#
# mode=phased:  one phased sample per individual
# mode=dup:     every phased sample written twice
# mode=haphom:  every haplotype written as a homozygous sample
# The last two have the same haplotype counts, and no unknown-phase
# double-hets, so --blocks must report the same blocks on both.
BEGIN {
    OFS = "\t";
    sample_ct = 150;
    block_ct = 12;
    per_block = 12;
    seed = 20261014;
    if (mode == "") {
        mode = "phased";
    }

    print "##fileformat=VCFv4.2";
    print "##contig=<ID=1>";
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">";
    line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (s = 0; s < sample_ct; ++s) {
        if (mode == "phased") {
            line = line "\tS" s;
        } else {
            line = line "\tS" s "a\tS" s "b";
        }
    }
    print line;

    for (b = 0; b < block_ct; ++b) {
        for (h = 0; h < 4; ++h) {
            for (j = 0; j < per_block; ++j) {
                hap[h, j] = (rnd() < 0.5)? 0 : 1;
            }
        }
        for (s = 0; s < sample_ct; ++s) {
            draw1[s] = int(rnd() * 4);
            draw2[s] = int(rnd() * 4);
            # occasional recombinant haplotype, so that not every pair is in
            # perfect LD
            brk1[s] = (rnd() < 0.15)? int(rnd() * per_block) : per_block;
            alt1[s] = int(rnd() * 4);
        }
        for (j = 0; j < per_block; ++j) {
            line = "1" OFS (b * 100000 + j * 2000 + 1) OFS "b" b "v" j OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS "GT";
            for (s = 0; s < sample_ct; ++s) {
                a1 = (j < brk1[s])? hap[draw1[s], j] : hap[alt1[s], j];
                a2 = hap[draw2[s], j];
                missing = (rnd() < 0.02);
                if (mode == "phased") {
                    gt = missing? "./." : (a1 "|" a2);
                } else if (mode == "dup") {
                    gt = missing? "./.\t./." : (a1 "|" a2 "\t" a1 "|" a2);
                } else {
                    gt = missing? "./.\t./." : (a1 "/" a1 "\t" a2 "/" a2);
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
