# Deterministic single-chromosome VCF with planted runs of homozygosity.
# Every third sample carries one long homozygous stretch; the rest are
# heterozygous at the usual rate, so the scan has something to find and
# something to reject.
BEGIN {
    OFS = "\t";
    sample_ct = 30;
    variant_ct = 6000;
    seed = 20250905;

    print "##fileformat=VCFv4.2";
    print "##contig=<ID=1>";
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">";
    line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (s = 0; s < sample_ct; ++s) {
        line = line "\tS" s;
        if (s % 3 == 0) {
            roh_start[s] = (s * 137) % 3000;
            roh_end[s] = roh_start[s] + 1200;
        } else {
            roh_start[s] = -1;
            roh_end[s] = -1;
        }
    }
    print line;

    for (v = 0; v < variant_ct; ++v) {
        line = "1" OFS ((v + 1) * 1000) OFS "v" v OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS "GT";
        for (s = 0; s < sample_ct; ++s) {
            r = rnd();
            if (r < 0.005) {
                gt = "./.";
            } else if (v >= roh_start[s] && v < roh_end[s]) {
                gt = (rnd() < 0.5)? "0/0" : "1/1";
            } else {
                r = rnd();
                gt = (r < 0.25)? "0/0" : ((r < 0.75)? "0/1" : "1/1");
            }
            line = line OFS gt;
        }
        print line;
    }
}

function rnd() {
    seed = (seed * 16807) % 2147483647;
    return seed / 2147483647;
}
