#!/bin/bash

# --neighbour: the PC-space outlier statistic from section 3.4 of Prive et al.
# (2020).

set -exo pipefail

fails() {
    "$@" && exit 1 || true
}

python3 simulate.py
plink --file tmp_nb --make-bed --out tmp_nb > /dev/null

# 1. The report agrees with a direct recomputation from the .eigenvec it was
#    derived from, at the default K.
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour --out plink2_nb
test -s plink2_nb.nearest
python3 oracle.py plink2_nb.eigenvec plink2_nb.nearest 5

# 2. K is adjustable, and 5 is the default.
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour 5 --out plink2_k5
diff -q plink2_nb.nearest plink2_k5.nearest
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour 12 --out plink2_k12
python3 oracle.py plink2_k12.eigenvec plink2_k12.nearest 12
fails diff -q plink2_nb.nearest plink2_k12.nearest

# 3. The PC count carries through: a different --pca gives a different report,
#    and the oracle still matches it.
$1/plink2 $2 $3 --bfile tmp_nb --pca 2 --neighbour --out plink2_pc2
python3 oracle.py plink2_pc2.eigenvec plink2_pc2.nearest 5
fails diff -q plink2_nb.nearest plink2_pc2.nearest

# 4. The planted outliers come out on top.  Each is homozygous across its own
#    block of variants, so each sits alone in PC space; a tight group of
#    identical samples would read as a small cluster instead.
grep -v '^#' plink2_nb.nearest | sort -g -r -k5,5 | head -4 | cut -f2 | sort > tmp_top4.txt
sort tmp_nb_outliers.txt > tmp_want.txt
diff -q tmp_top4.txt tmp_want.txt

# 5. cols= drops columns, and the remaining values are unchanged.
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour cols=stat --out plink2_cols
head -1 plink2_cols.nearest | grep -qx '#IID	STAT'
python3 oracle.py plink2_nb.eigenvec plink2_cols.nearest 5
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour cols=fid,sid,distself --out plink2_cols2
head -1 plink2_cols2.nearest | grep -qx '#FID	IID	SID	DIST_SELF'

# 6. 'zs' compresses the report, and the contents survive the round trip.
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour zs --out plink2_zs
test -s plink2_zs.nearest.zst
$1/plink2 $2 $3 --zst-decompress plink2_zs.nearest.zst > tmp_unzs.nearest
diff -q tmp_unzs.nearest plink2_nb.nearest

# 7. --neighbor is the same flag.
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbor --out plink2_alias
diff -q plink2_alias.nearest plink2_nb.nearest

# 8. It needs PC coordinates, and says so.  (--neighbour is parsed before --pca
#    and --read-eigvec, so this check cannot live in the parser.)
fails $1/plink2 $2 $3 --bfile tmp_nb --neighbour --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_nb --neighbour --out plink2_bad > tmp_err.txt 2>&1 || true
grep -q 'must be used with --pca or --read-eigvec' tmp_err.txt

# 8b. --read-eigvec replays a previous --pca, and lands on the same report.
#     The coordinates are re-read from text, so they arrive rounded; compare
#     against the same recomputation the oracle does rather than byte for byte.
$1/plink2 $2 $3 --bfile tmp_nb --read-eigvec plink2_nb.eigenvec --neighbour --out plink2_re
python3 oracle.py plink2_nb.eigenvec plink2_re.nearest 5
$1/plink2 $2 $3 --bfile tmp_nb --read-eigvec plink2_nb.eigenvec --neighbour 12 --out plink2_re12
python3 oracle.py plink2_nb.eigenvec plink2_re12.nearest 12

# 8c. PLINK 1.9's .eigenvec is headerless and space-delimited, and works too.
grep -v '^#' plink2_nb.eigenvec | tr '\t' ' ' > tmp_19.eigenvec
$1/plink2 $2 $3 --bfile tmp_nb --read-eigvec tmp_19.eigenvec --neighbour --out plink2_re19
diff -q plink2_re19.nearest plink2_re.nearest

# 8d. Samples not in the current fileset are ignored, so one .eigenvec can be
#     replayed against a subset.
head -30 tmp_nb.fam | awk '{print $1 "\t" $2}' > tmp_keep.txt
$1/plink2 $2 $3 --bfile tmp_nb --keep tmp_keep.txt --read-eigvec plink2_nb.eigenvec --neighbour --out plink2_re_keep
test $(grep -vc '^#' plink2_re_keep.nearest) -eq 30

# 8e. A sample with no entry is an error, and so is a duplicate entry.
head -1 plink2_nb.eigenvec > tmp_short.eigenvec
grep -v '^#' plink2_nb.eigenvec | head -30 >> tmp_short.eigenvec
fails $1/plink2 $2 $3 --bfile tmp_nb --read-eigvec tmp_short.eigenvec --neighbour --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_nb --read-eigvec tmp_short.eigenvec --neighbour --out plink2_bad > tmp_err3.txt 2>&1 || true
grep -q 'lack a --read-eigvec entry' tmp_err3.txt
cat plink2_nb.eigenvec > tmp_dup.eigenvec
grep -v '^#' plink2_nb.eigenvec | head -1 >> tmp_dup.eigenvec
fails $1/plink2 $2 $3 --bfile tmp_nb --read-eigvec tmp_dup.eigenvec --neighbour --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_nb --read-eigvec tmp_dup.eigenvec --neighbour --out plink2_bad > tmp_err4.txt 2>&1 || true
grep -q 'Duplicate sample ID' tmp_err4.txt

# 8f. A ragged or non-numeric file is rejected.
sed '3s/\t[^\t]*$//' plink2_nb.eigenvec > tmp_ragged.eigenvec
fails $1/plink2 $2 $3 --bfile tmp_nb --read-eigvec tmp_ragged.eigenvec --neighbour --out plink2_bad
sed '3s/\t[^\t]*$/\tabc/' plink2_nb.eigenvec > tmp_nonnum.eigenvec
fails $1/plink2 $2 $3 --bfile tmp_nb --read-eigvec tmp_nonnum.eigenvec --neighbour --out plink2_bad

# 8g. The two sources are mutually exclusive, and --read-eigvec needs
#     --neighbour.
fails $1/plink2 $2 $3 --bfile tmp_nb --pca 4 --read-eigvec plink2_nb.eigenvec --neighbour --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_nb --read-eigvec plink2_nb.eigenvec --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_nb --read-eigvec plink2_nb.eigenvec --out plink2_bad > tmp_err5.txt 2>&1 || true
grep -q 'must be used with --neighbour' tmp_err5.txt

# 9. K has to leave room for K neighbours.
fails $1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour 60 --out plink2_bad
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour 60 --out plink2_bad > tmp_err2.txt 2>&1 || true
grep -q 'must be smaller than the number of' tmp_err2.txt
$1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour 59 --out plink2_k59
test -s plink2_k59.nearest

# 10. Bad arguments are rejected.
fails $1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour 0 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour abc --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour cols=nonsense --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_nb --pca 4 --neighbour cols=stat cols=distnn --out plink2_bad
