# Compares PLINK 1.9's fixed-width --twolocus report against plink2's table.
#
# 1.9 writes, per group, a count matrix and then a frequency matrix, both with
# '*/*' marginal rows and columns; plink2 writes one row per (group, genotype,
# genotype) cell and leaves the marginals to the reader.  Genotype labels are
# also not consistently ordered between the two loci in 1.9's output, so both
# sides are normalized by sorting the two alleles.
function norm(gt,   a) {
    if (gt == "0/0" || gt == ".") { return "." }
    split(gt, a, "/");
    return (a[1] <= a[2])? (a[1] "/" a[2]) : (a[2] "/" a[1]);
}
FNR == NR {
    if ($0 ~ /^All individuals$/) { group = "ALL"; matrix_idx = 0; in_matrix = 0; next }
    if ($0 ~ /^Cases$/) { group = "CASE"; matrix_idx = 0; in_matrix = 0; next }
    if ($0 ~ /^Controls$/) { group = "CTRL"; matrix_idx = 0; in_matrix = 0; next }
    if (group == "") { next }
    if (NF == 1 && $1 == id2) { ++matrix_idx; expect_header = 1; in_matrix = 0; next }
    if (expect_header) {
        for (i = 1; i <= NF; ++i) { col[i] = norm($i) }
        col_ct = NF;
        expect_header = 0;
        in_matrix = 1;
        next
    }
    if (!in_matrix) { next }
    if (NF == 0) { in_matrix = 0; next }
    # Only the first matrix in each group holds counts.
    if (matrix_idx != 1) { next }
    first = 1;
    if ($1 == id1) { first = 2 }
    row_gt = norm($first);
    if (row_gt == "*/*") { next }
    for (i = 1; i <= col_ct; ++i) {
        if (col[i] == "*/*") { continue }
        want[group "|" row_gt "|" col[i]] = $(first + i);
        ++n1;
    }
    next
}
/^#/ { next }
{
    ++n2;
    key = $1 "|" norm($3) "|" norm($5);
    if (!(key in want)) { print "cell missing from the PLINK 1.9 report: " key; failed = 1; exit 1 }
    if (want[key] != $6) { print "count differs on " key ": " want[key] " vs " $6; failed = 1; exit 1 }
    delete want[key];
}
END {
    if (failed) { exit 1 }
    if (n1 != n2) { print "cell count mismatch: " n1 " vs " n2; exit 1 }
    if (n1 == 0) { print "no cells compared"; exit 1 }
    print n1 " cells matched"
}
