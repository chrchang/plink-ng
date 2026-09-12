#!/bin/bash

# --set/--set-names/--set-collapse-all/--write-set, checked by round-tripping
# through --write-set, since PLINK 1.9's --write-set output is a different
# format (it writes a variant-per-line table under --write-set-r2 style
# options, and its plain --write-set has no plink2 counterpart yet).

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

plink --simulate simulate.txt --simulate-ncases 20 --simulate-ncontrols 20 --out tmp_data > /dev/null

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

# 1. Round trip: the report is the same sets, in file order, with the unknown
#    ID dropped and the empty set kept.
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

# 2. --set-names keeps the named sets, and only those, in file order.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names geneC,geneA --write-set --out plink2_names
printf 'geneA\nsv_0\nsv_1\nsv_20\nEND\n\ngeneC\nEND\n\n' > tmp_names.set
diff -q tmp_names.set plink2_names.set

# 3. --set-collapse-all replaces every set with their union under one name, so
#    the members come out sorted and deduplicated.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-collapse-all EVERYTHING --write-set --out plink2_all
head -n 1 plink2_all.set | grep -qx 'EVERYTHING'
awk 'NF && $1 != "EVERYTHING" && $1 != "END"' plink2_all.set > tmp_all_members.txt
awk 'NF && $1 != "END" && $1 != "nosuchvariant" && $1 !~ /^gene/' tmp_sets.txt | sort -u > tmp_union.txt
diff <(sort -u tmp_all_members.txt) tmp_union.txt
test "$(wc -l < tmp_all_members.txt)" -eq "$(wc -l < tmp_union.txt)"

# 4. A set is defined over the variants that survived filtering, so a filter
#    takes its variants out of the sets too.
printf 'sv_0\nsv_2\n' > tmp_exclude.txt
$1/plink2 $2 $3 --bfile tmp_data --exclude tmp_exclude.txt --set tmp_sets.txt --write-set --out plink2_filt
fails grep -qx 'sv_0' plink2_filt.set
fails grep -qx 'sv_2' plink2_filt.set
grep -qx 'sv_1' plink2_filt.set
grep -qx 'sv_3' plink2_filt.set

# 5. 'zs' is the same report, compressed.
$1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --write-set zs --out plink2_zs
$1/plink2 $2 $3 --zst-decompress plink2_zs.set.zst > plink2_zs.set
diff -q plink2_rt.set plink2_zs.set

# 6. Malformed files are rejected: an unmatched END, and a set with no END.
printf 'geneA\nsv_0\nEND\nEND\n' > tmp_extra_end.txt
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_extra_end.txt --write-set --out plink2_bad
printf 'geneA\nsv_0\n' > tmp_no_end.txt
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_no_end.txt --write-set --out plink2_bad

# 7. Rejected: a --set-names that matches nothing, and the flag dependencies.
fails $1/plink2 $2 $3 --bfile tmp_data --set tmp_sets.txt --set-names nosuchset --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --write-set --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set-names geneA --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --set-collapse-all ALL --out plink2_bad
