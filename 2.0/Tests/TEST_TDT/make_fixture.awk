# 150 trios (father, mother, affected child) with Mendelian transmission, plus
# missing calls and a few deliberate Mendel errors, written as a VCF plus the
# .fam columns the pedigree needs.
BEGIN {
    OFS = "\t";
    trio_ct = 150;
    variant_ct = 400;
    seed = 20250905;

    print "##fileformat=VCFv4.2" > "tmp_trios.vcf";
    print "##contig=<ID=1>" > "tmp_trios.vcf";
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">" > "tmp_trios.vcf";
    line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (t = 0; t < trio_ct; ++t) {
        line = line "\tF" t "_pa\tF" t "_ma\tF" t "_kid";
    }
    print line > "tmp_trios.vcf";

    for (v = 0; v < variant_ct; ++v) {
        maf = 0.1 + 0.35 * rnd();
        line = "1" OFS (v * 5000 + 1) OFS "v" v OFS "A" OFS "G" OFS "." OFS "." OFS "." OFS "GT";
        for (t = 0; t < trio_ct; ++t) {
            pa1 = (rnd() < maf)? 1 : 0; pa2 = (rnd() < maf)? 1 : 0;
            ma1 = (rnd() < maf)? 1 : 0; ma2 = (rnd() < maf)? 1 : 0;
            kid1 = (rnd() < 0.5)? pa1 : pa2;
            kid2 = (rnd() < 0.5)? ma1 : ma2;
            r = rnd();
            if (r < 0.01) {
                # deliberate Mendel error: flip the child's paternal allele
                kid1 = 1 - kid1;
            }
            pa = pa1 "/" pa2; ma = ma1 "/" ma2; kid = kid1 "/" kid2;
            if (rnd() < 0.01) { pa = "./." }
            if (rnd() < 0.01) { ma = "./." }
            if (rnd() < 0.01) { kid = "./." }
            line = line OFS pa OFS ma OFS kid;
        }
        print line > "tmp_trios.vcf";
    }

    for (t = 0; t < trio_ct; ++t) {
        print "F" t, "pa", 0, 0, 1, 1 > "tmp_ped.txt";
        print "F" t, "ma", 0, 0, 2, 1 > "tmp_ped.txt";
        print "F" t, "kid", "pa", "ma", (t % 2) + 1, 2 > "tmp_ped.txt";
    }
}

function rnd() {
    seed = (seed * 16807) % 2147483647;
    return seed / 2147483647;
}
