# Splits a VCF from make_vcf_dprime.awk for the --flip-scan-ref-pfile check.
# The first half of the samples go to tmp_ref_main.vcf, the second half to
# tmp_ref_panel.vcf, and tmp_ref_merged.vcf gets what merging the two should
# produce.  The reference copy is made awkward on purpose:
#   * every 5th variant has REF and ALT swapped there;
#   * every 7th variant gains a third allele there (some ALT calls become T),
#     so it is biallelic in one fileset and triallelic in the other;
#   * every 11th variant is missing from it entirely.
# None of this may change the answer relative to the merged fileset.
BEGIN {
    OFS = "\t";
}
/^##/ {
    print > "tmp_ref_main.vcf";
    print > "tmp_ref_panel.vcf";
    print > "tmp_ref_merged.vcf";
    next;
}
/^#CHROM/ {
    half = (NF - 9) / 2;
    main_line = $1;
    panel_line = $1;
    for (i = 2; i <= 9; ++i) {
        main_line = main_line OFS $i;
        panel_line = panel_line OFS $i;
    }
    for (i = 10; i < 10 + half; ++i) {
        main_line = main_line OFS $i;
    }
    for (i = 10 + half; i <= NF; ++i) {
        panel_line = panel_line OFS $i;
    }
    print main_line > "tmp_ref_main.vcf";
    print panel_line > "tmp_ref_panel.vcf";
    print > "tmp_ref_merged.vcf";
    next;
}
{
    ++vidx;
    swap = (vidx % 5 == 0);
    tri = (!swap) && (vidx % 7 == 0);
    absent = (vidx % 11 == 0);
    main_line = $1;
    for (i = 2; i <= 9; ++i) {
        main_line = main_line OFS $i;
    }
    for (i = 10; i < 10 + half; ++i) {
        main_line = main_line OFS $i;
    }
    print main_line > "tmp_ref_main.vcf";

    merged_line = $1 OFS $2 OFS $3 OFS $4 OFS (tri? "G,T" : $5);
    for (i = 6; i <= 9; ++i) {
        merged_line = merged_line OFS $i;
    }
    for (i = 10; i < 10 + half; ++i) {
        merged_line = merged_line OFS $i;
    }
    panel_line = $1 OFS $2 OFS $3 OFS (swap? $5 : $4) OFS (swap? $4 : (tri? "G,T" : $5));
    for (i = 6; i <= 9; ++i) {
        panel_line = panel_line OFS $i;
    }
    for (i = 10 + half; i <= NF; ++i) {
        gt = $i;
        if (tri && ((i % 3) == 0)) {
            gsub(/1/, "2", gt);
        }
        merged_line = merged_line OFS (absent? "./." : gt);
        if (swap) {
            gsub(/0/, "x", gt);
            gsub(/1/, "0", gt);
            gsub(/x/, "1", gt);
        }
        panel_line = panel_line OFS gt;
    }
    print merged_line > "tmp_ref_merged.vcf";
    if (!absent) {
        print panel_line > "tmp_ref_panel.vcf";
    }
}
