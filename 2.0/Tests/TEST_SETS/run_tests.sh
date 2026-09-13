#!/bin/bash

# --set/--make-set and the rest of PLINK 1.x's set-manipulation flags.
# PLINK 1.9's --write-set writes the same format --set reads, so most of these
# checks are direct diffs against 1.9's output; the exceptions are noted where
# they appear.

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

plink --simulate simulate.txt --simulate-ncases 20 --simulate-ncontrols 20 --out tmp_data > /dev/null

# sv_0..sv_29, at positions 1..30 of chromosome 1.

# geneA is three scattered variants, geneB two adjacent ones plus an ID that is
# not in the fileset, geneC is empty, and geneD spans a whole run, so both
# setdef forms get exercised.
{
    printf 'geneA\nsv_0\nsv_1\nsv_20\nEND\n\n'
    printf 'geneB\nsv_2\nsv_3\nnosuchvariant\nEND\n\n'
    printf 'geneC\nEND\n\n'
    printf 'geneD\n'
    for i in $(seq 5 19); do printf 'sv_%u\n' $i; done
    printf 'END\n'
} > tmp_sets.txt

# ranges for --make-set: two groups, and a one-variant set.
{
    printf '1 2 4 geneA grp1\n'
    printf '1 10 12 geneB grp1\n'
    printf '1 20 20 geneC grp2\n'
} > tmp_ranges.txt

# 1. Round trip: the report is the same sets, in file order, with the unknown
#    ID dropped and the empty set kept.  1.9 writes the same thing.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --write-set --out plink2_rt
{
    printf 'geneA\nsv_0\nsv_1\nsv_20\nEND\n\n'
    printf 'geneB\nsv_2\nsv_3\nEND\n\n'
    printf 'geneC\nEND\n\n'
    printf 'geneD\n'
    for i in $(seq 5 19); do printf 'sv_%u\n' $i; done
    printf 'END\n\n'
} > tmp_expected.set
diff -q tmp_expected.set plink2_rt.set
plink --bfile tmp_data --set tmp_sets.txt --write-set --out plink19_rt > /dev/null
diff -q plink19_rt.set plink2_rt.set

# 2. --set-names keeps the named sets, and only those, in file order.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names geneC,geneA --write-set --out plink2_names
printf 'geneA\nsv_0\nsv_1\nsv_20\nEND\n\ngeneC\nEND\n\n' > tmp_names.set
diff -q tmp_names.set plink2_names.set

# 3. --subset takes the same list from a file, and the two combine.
printf 'geneC\ngeneA\n' > tmp_subset.txt
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --subset tmp_subset.txt --write-set --out plink2_subset
diff -q tmp_names.set plink2_subset.set
plink --bfile tmp_data --set tmp_sets.txt --subset tmp_subset.txt --write-set --out plink19_subset > /dev/null
diff -q plink19_subset.set plink2_subset.set
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --subset tmp_subset.txt --set-names geneB --write-set --out plink2_subset2
plink --bfile tmp_data --set tmp_sets.txt --subset tmp_subset.txt --set-names geneB --write-set --out plink19_subset2 > /dev/null
diff -q plink19_subset2.set plink2_subset2.set

# 4. --set-collapse-all replaces every set with their union under one name, so
#    the members come out sorted and deduplicated.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-collapse-all EVERYTHING --write-set --out plink2_all
head -n 1 plink2_all.set | grep -qx 'EVERYTHING'
awk 'NF && $1 != "EVERYTHING" && $1 != "END"' plink2_all.set > tmp_all_members.txt
awk 'NF && $1 != "END" && $1 != "nosuchvariant" && $1 !~ /^gene/' tmp_sets.txt | sort -u > tmp_union.txt
diff <(sort -u tmp_all_members.txt) tmp_union.txt
test "$(wc -l < tmp_all_members.txt)" -eq "$(wc -l < tmp_union.txt)"

# 5. A set is defined over the variants that survived filtering, so a filter
#    takes its variants out of the sets too.
printf 'sv_0\nsv_2\n' > tmp_exclude.txt
$1/plink2 $2 $3 --bfile tmp_data --exclude tmp_exclude.txt --set tmp_sets.txt --write-set --out plink2_filt
fails grep -qx 'sv_0' plink2_filt.set
fails grep -qx 'sv_2' plink2_filt.set
grep -qx 'sv_1' plink2_filt.set
grep -qx 'sv_3' plink2_filt.set

# 6. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --write-set zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.set.zst > plink2_zs.set
diff -q plink2_rt.set plink2_zs.set

# 7. --make-set defines the same kind of sets from bp ranges, with the set
#    names sorted and deduplicated.  --make-set-border stretches each range,
#    and --make-set-collapse-group groups by the fifth column instead.
for extra_flags in "" "--make-set-border 0.002" "--make-set-collapse-group" "--make-set-collapse-group --make-set-border 0.002"; do
    $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt $extra_flags --write-set --out plink2_ms
    plink --bfile tmp_data --make-set tmp_ranges.txt $extra_flags --write-set --out plink19_ms > /dev/null
    diff -q plink19_ms.set plink2_ms.set
done
printf 'geneA\nsv_1\nsv_2\nsv_3\nEND\n\ngeneB\nsv_9\nsv_10\nsv_11\nEND\n\ngeneC\nsv_19\nEND\n\n' > tmp_ms_expected.set
$1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt --write-set --out plink2_ms
diff -q tmp_ms_expected.set plink2_ms.set

# 8. --subset and --set-names also select --make-set sets, by set ID even when
#    the sets themselves are groups.
printf 'geneA\n' > tmp_subset_ms.txt
for extra_flags in "--subset tmp_subset_ms.txt" "--set-names geneB" "--make-set-collapse-group --subset tmp_subset_ms.txt"; do
    $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt $extra_flags --write-set --out plink2_ms_sub
    plink --bfile tmp_data --make-set tmp_ranges.txt $extra_flags --write-set --out plink19_ms_sub > /dev/null
    diff -q plink19_ms_sub.set plink2_ms_sub.set
done

# 9. --complement-sets inverts each set and renames it, and
#    --make-set-complement-all/-group are the collapsing variants of that.
for extra_flags in "--complement-sets" "--make-set-complement-all NOTME" "--make-set-complement-group" "--set-collapse-all ALL"; do
    $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt $extra_flags --write-set --out plink2_comp
    plink --bfile tmp_data --make-set tmp_ranges.txt $extra_flags --write-set --out plink19_comp > /dev/null
    diff -q plink19_comp.set plink2_comp.set
done

# 10. --complement-sets with --set is checked against the complement computed
#     here instead of against 1.9, for two reasons: 1.9 does not apply the
#     'C_' prefix on this path (it only does so for --make-set, though its help
#     text promises the prefix outright), and it drops the first variant from
#     each complement whose set does not contain it.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --complement-sets --write-set --out plink2_setcomp
grep -qx 'C_geneA' plink2_setcomp.set
grep -qx 'C_geneD' plink2_setcomp.set
awk '/^C_geneB$/{flag=1; next} /^END$/{flag=0} flag' plink2_setcomp.set > tmp_c_geneb.txt
{ printf 'sv_0\nsv_1\n'; for i in $(seq 4 29); do printf 'sv_%u\n' $i; done; } > tmp_c_geneb_expected.txt
diff -q tmp_c_geneb_expected.txt tmp_c_geneb.txt
# C_geneC is the complement of an empty set, so it is every variant.
test "$(awk '/^C_geneC$/{flag=1; next} /^END$/{flag=0} flag' plink2_setcomp.set | wc -l)" -eq 30

# 11. --set-table writes the variant-by-set membership table 1.9 writes.
for set_flags in "--set tmp_sets.txt" "--make-set tmp_ranges.txt" "--make-set tmp_ranges.txt --complement-sets"; do
    $1/plink2 $2 $3 --bfile tmp_data $set_flags --set-table --out plink2_tbl
    plink --bfile tmp_data $set_flags --set-table --out plink19_tbl > /dev/null
    diff -q plink19_tbl.set.table plink2_tbl.set.table
done
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-table zs --out plink2_tbl_zs
$1/plink2 $2 $3 --zst-decompress plink2_tbl_zs.set.table.zst > plink2_tbl_zs.set.table
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-table --out plink2_tbl
diff -q plink2_tbl.set.table plink2_tbl_zs.set.table

# 12. Malformed files are rejected: an unmatched END, and a set with no END.
printf 'geneA\nsv_0\nEND\nEND\n' > tmp_extra_end.txt
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_extra_end.txt --write-set --out plink2_bad
printf 'geneA\nsv_0\n' > tmp_no_end.txt
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_no_end.txt --write-set --out plink2_bad
# a --make-set line without a set ID, and one naming an unknown chromosome
printf '1 2 4\n' > tmp_short_range.txt
fails $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_short_range.txt --write-set --out plink2_bad
printf '1 2 4 geneA\nnosuchchr 2 4 geneB\n' > tmp_bad_chr.txt
fails $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_bad_chr.txt --write-set --out plink2_bad
# --make-set-collapse-group needs the fifth column
printf '1 2 4 geneA\n' > tmp_nogroup.txt
fails $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_nogroup.txt --make-set-collapse-group --write-set --out plink2_bad

# 13. Rejected: selections that match nothing, and the flag dependencies.
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names nosuchset --write-set --out plink2_bad
printf 'nosuchset\n' > tmp_nosuch.txt
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --subset tmp_nosuch.txt --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt --subset tmp_nosuch.txt --write-set --out plink2_bad
: > tmp_empty.txt
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --subset tmp_empty.txt --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set-table --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set-names geneA --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set-collapse-all ALL --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --subset tmp_subset.txt --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --complement-sets --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --make-set-complement-all ALL --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --make-set-border 5 --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --make-set-collapse-group --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --make-set tmp_ranges.txt --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt --complement-sets --make-set-complement-all ALL --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt --make-set-collapse-group --make-set-complement-group --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --make-set tmp_ranges.txt --make-set-collapse-group --set-collapse-all ALL --write-set --out plink2_bad
