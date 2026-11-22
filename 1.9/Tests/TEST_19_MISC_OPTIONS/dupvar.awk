# Writes the --list-duplicate-vars report expected for a sorted .bim file:
# groups of 2+ variants with the same chromosome, position and alleles, in
# file order, with IDs in file order.
#
# Usage: awk [-v same_ref=1] [-v ids_only=1] [-v suppress_first=1]
#            -f dupvar.awk x.bim
#   same_ref: alleles must match as (A2, A1) pairs instead of as sets; the
#             report then has REF (A2) and ALT (A1) columns.
#   ids_only: print just the IDs.
#   suppress_first: leave out each group's first ID.

{
    if (same_ref) {
        key = $1 " " $4 " " $6 " " $5
        alleles = $6 "\t" $5
    } else if ($5 < $6) {
        key = $1 " " $4 " " $5 " " $6
        alleles = $5 "," $6
    } else {
        key = $1 " " $4 " " $6 " " $5
        alleles = $6 "," $5
    }
    if (!(key in ct)) {
        order[++key_ct] = key
        head[key] = $1 "\t" $4 "\t" alleles
    }
    ids[key, ++ct[key]] = $2
}

END {
    if (!ids_only) print same_ref? "CHR\tPOS\tREF\tALT\tIDS" : "CHR\tPOS\tALLELES\tIDS"
    for (k = 1; k <= key_ct; ++k) {
        key = order[k]
        if (ct[key] < 2) continue
        line = ""
        for (i = suppress_first? 2 : 1; i <= ct[key]; ++i) line = line ((line == "")? "" : " ") ids[key, i]
        print ids_only? line : head[key] "\t" line
    }
}
