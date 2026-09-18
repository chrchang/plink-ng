# Writes tmp_multi.vcf, where every other variant is triallelic, and
# tmp_pooled.vcf, the same calls with each variant's nonmajor alleles merged
# into one.  --test-mishap handles multiallelic variants as major vs. rest, so
# the two reports must carry the same numbers.
BEGIN {
    OFS = "\t";
    sample_ct = 400;
    variant_ct = 30;
    seed = 20260916;

    header = "##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (s = 0; s < sample_ct; ++s) {
        header = header "\tS" s;
    }
    print header > "tmp_multi.vcf";
    print header > "tmp_pooled.vcf";
    split("A C G", allele_names, " ");

    for (v = 0; v < variant_ct; ++v) {
        allele_ct = (v % 2)? 2 : 3;
        weight_sum = 0;
        for (a = 0; a < allele_ct; ++a) {
            weight[a] = 0.1 + rnd();
            weight_sum += weight[a];
            allele_cts[a] = 0;
        }
        for (s = 0; s < sample_ct; ++s) {
            if (rnd() < 0.08) {
                missing[s] = 1;
                continue;
            }
            missing[s] = 0;
            call1[s] = draw();
            call2[s] = draw();
            ++allele_cts[call1[s]];
            ++allele_cts[call2[s]];
        }
        maj = 0;
        for (a = 1; a < allele_ct; ++a) {
            if (allele_cts[a] > allele_cts[maj]) {
                maj = a;
            }
        }
        alts = allele_names[2];
        if (allele_ct == 3) {
            alts = alts "," allele_names[3];
        }
        # The pooled file puts the major allele first and the pooled rest second.
        multi_line = "1" OFS ((v + 1) * 1000) OFS "v" v OFS allele_names[1] OFS alts OFS "." OFS "." OFS "." OFS "GT";
        pooled_line = "1" OFS ((v + 1) * 1000) OFS "v" v OFS allele_names[maj + 1] OFS "T" OFS "." OFS "." OFS "." OFS "GT";
        for (s = 0; s < sample_ct; ++s) {
            if (missing[s]) {
                multi_line = multi_line OFS "./.";
                pooled_line = pooled_line OFS "./.";
            } else {
                multi_line = multi_line OFS call1[s] "/" call2[s];
                pooled_line = pooled_line OFS (call1[s] != maj) "/" (call2[s] != maj);
            }
        }
        print multi_line > "tmp_multi.vcf";
        print pooled_line > "tmp_pooled.vcf";
    }
}

function draw(    r, a) {
    r = rnd() * weight_sum;
    for (a = 0; a < allele_ct - 1; ++a) {
        if (r < weight[a]) {
            return a;
        }
        r -= weight[a];
    }
    return allele_ct - 1;
}

function rnd() {
    seed = (seed * 16807) % 2147483647;
    return seed / 2147483647;
}
