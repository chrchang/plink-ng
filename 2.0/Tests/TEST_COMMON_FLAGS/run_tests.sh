#!/bin/bash

set -exo pipefail

# Covers common flags which no other test mentions.  Most comparisons are with
# plink 1.9, or with values known from how make_data.py built the data.
#
# The end of the file covers some less common flags: --vscore-col-nums,
# --loop-cats, --error-on-freq-calc, --warning-errcode, and
# --allow-misleading-out-arg.

fails() {
    "$@" && exit 1 || true
}

P2="$1/plink2 $2 $3"

python3 make_data.py
plink --file tmp_data --make-bed --out tmp_data
plink --file tmp_sex --make-bed --out tmp_sex

# --hwe: the 8 all-heterozygous variants, v7, v47, ..., are gone.
plink --bfile tmp_data --hwe 0.01 --write-snplist --out plink_hwe
$P2 --bfile tmp_data --hwe 0.01 0 --write-snplist --out plink2_hwe
diff plink_hwe.snplist plink2_hwe.snplist
test $(wc -l < plink2_hwe.snplist) -gt 200
test $(grep -cE '^v(7|47|87|127|167|207|247|287)$' plink2_hwe.snplist) -eq 0

# --indep-pairwise: each near-copy pair loses one member.
plink --bfile tmp_data --indep-pairwise 50 5 0.5 --out plink_prune
$P2 --bfile tmp_data --indep-pairwise 50 5 0.5 --out plink2_prune
diff plink_prune.prune.in plink2_prune.prune.in
diff plink_prune.prune.out plink2_prune.prune.out
test $(cat plink2_prune.prune.in plink2_prune.prune.out | sort -u | wc -l) -eq 300
test $(awk 'NR == FNR { kept[$1] = 1; next } ($1 in kept) && ($2 in kept)' plink2_prune.prune.in tmp_ld_pairs.txt | wc -l) -eq 0

# --mendel: exactly the injected errors, all with the code for 0/0 x 0/0 -> 1/1.
plink --bfile tmp_data --mendel --out plink_mendel
$P2 --bfile tmp_data --mendel --out plink2_mendel
diff <(sort tmp_mendel_truth.txt) <(tail -n +2 plink2_mendel.mendel | cut -f 1,4 | tr '\t' ' ' | sort)
test $(tail -n +2 plink2_mendel.mendel | cut -f 5 | sort -u) = 5
diff <(tail -n +2 plink_mendel.lmendel | awk '{print $2, $3}') <(tail -n +2 plink2_mendel.lmendel | cut -f 2,3 | tr '\t' ' ')
diff <(tail -n +2 plink_mendel.imendel | awk '{print $1, $2, $3}') <(tail -n +2 plink2_mendel.imendel | tr '\t' ' ')

# --keep-founders
$P2 --bfile tmp_data --keep-founders --make-just-psam --out plink2_founders
diff <(sort tmp_founders.txt) <(tail -n +2 plink2_founders.psam | cut -f 1,2 | tr '\t' ' ' | sort)

# --impute-sex: fixes the two swapped samples and fills in the missing one.
plink --bfile tmp_sex --impute-sex 0.2 0.8 --make-bed --out plink_sex
$P2 --bfile tmp_sex --impute-sex max-female-xf=0.2 min-male-xf=0.8 --make-just-psam --out plink2_sex
diff tmp_sex_truth.txt <(tail -n +2 plink2_sex.psam | tr '\t' ' ')
diff tmp_sex_truth.txt <(cut -d ' ' -f 1,2,5 plink_sex.fam)
test $(grep -c PROBLEM plink2_sex.sexcheck) -eq 3

# --make-grm-bin
plink --bfile tmp_data --make-grm-bin --out plink_grm
$P2 --bfile tmp_data --make-grm-bin --out plink2_grm
diff plink_grm.grm.id plink2_grm.grm.id
python3 check.py grm

# --update-ids
printf 't0 pa T0 PA\nt0 ma T0 MA\nu0 s U0 S\nnope x y z\n' > tmp_ids.txt
plink --bfile tmp_data --update-ids tmp_ids.txt --make-just-fam --out plink_ids
$P2 --bfile tmp_data --update-ids tmp_ids.txt --make-bed --out plink2_ids
diff <(tr ' ' '\t' < plink_ids.fam) plink2_ids.fam
test "$(head -n 1 plink2_ids.fam | cut -f 1,2)" = "T0	PA"

# --update-alleles: v2's old alleles do not match, so it is reported instead.
printf 'v0 A C G T\nv1 C A T G\nv2 A G T T\n' > tmp_alleles.txt
plink --bfile tmp_data --update-alleles tmp_alleles.txt --make-just-bim --out plink_alleles
$P2 --bfile tmp_data --update-alleles tmp_alleles.txt --make-bed --out plink2_alleles
diff plink_alleles.bim plink2_alleles.bim
diff plink_alleles.allele.no.snp plink2_alleles.allele.no.snp
test "$(head -n 3 plink2_alleles.bim | cut -f 5,6 | tr '\t\n' '  ')" = "T G T G C A "
test "$(cut -f 1 plink2_alleles.allele.no.snp)" = v2

# --fst: Hudson is recomputed by hand, and WC matches plink 1.9's --fst.
awk '{print $1, $2, (NR % 2)? "P1" : "P2"}' tmp_data.fam > tmp_within.txt
{ echo '#FID IID POP'; cat tmp_within.txt; } > tmp_pop.txt
$P2 --bfile tmp_data --pheno tmp_pop.txt --fst POP report-variants --out plink2_fst
python3 check.py fst
plink --bfile tmp_data --within tmp_within.txt --fst --out plink_fst
$P2 --bfile tmp_data --pheno tmp_pop.txt --fst POP method=wc report-variants --out plink2_fstwc
paste <(tail -n +2 plink_fst.fst | cut -f 2,5) <(tail -n +2 plink2_fstwc.P1.P2.fst.var | cut -f 3,5) \
  | awk '$1 != $3 { exit 1 } ($2 - $4 > 1e-4) || ($4 - $2 > 1e-4) { exit 1 }'

# --score-list
printf 'v0 C 1\nv2 C 2\n' > tmp_score1.txt
printf 'v2 C 3\nv3 A 1\n' > tmp_score2.txt
printf 'tmp_score1.txt\ntmp_score2.txt\n' > tmp_scores.txt
$P2 --bfile tmp_data --score-list tmp_scores.txt --out plink2_score
python3 check.py score

# Four samples, small enough to work out by hand.  IDs a and c appear twice: a
# with identical genotypes, c with different ones.
{
  printf '##fileformat=VCFv4.3\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n'
  printf '1\t10\ta\tA\tC\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/0\n'
  printf '1\t20\tb\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\t0/0\t./.\n'
  printf '1\t30\ta\tA\tC\t.\t.\t.\tGT\t0/0\t0/1\t1/1\t0/0\n'
  printf '1\t40\tc\tC\tG\t.\t.\t.\tGT\t1/1\t1/1\t0/1\t0/0\n'
  printf '1\t50\tc\tC\tG\t.\t.\t.\tGT\t0/0\t1/1\t0/1\t0/0\n'
} > tmp_small.vcf

pvar_pos() {
    grep -v '^#' "$1" | cut -f 2 | tr '\n' ' '
}

# --rm-dup: the default errors out.  'a' counts as a mismatch too, because its
# two records have different positions.
fails $P2 --vcf tmp_small.vcf --rm-dup --make-just-pvar --out plink2_dup_err
test "$(tr '\n' ' ' < plink2_dup_err.rmdup.mismatch)" = "a c "
$P2 --vcf tmp_small.vcf --rm-dup force-first list --make-just-pvar --out plink2_dup_first
test "$(pvar_pos plink2_dup_first.pvar)" = "10 20 40 "
test "$(tr '\n' ' ' < plink2_dup_first.rmdup.list)" = "a c "
$P2 --vcf tmp_small.vcf --rm-dup exclude-mismatch --make-just-pvar --out plink2_dup_ex
test "$(pvar_pos plink2_dup_ex.pvar)" = "20 "
$P2 --vcf tmp_small.vcf --rm-dup retain-mismatch --make-just-pvar --out plink2_dup_keep
test "$(pvar_pos plink2_dup_keep.pvar)" = "10 20 30 40 50 "

# --sample-diff
$P2 --vcf tmp_small.vcf --sample-diff base=s1 s2 s3 --out plink2_diff
test "$(tail -n +2 plink2_diff.sdiff.summary | tr '\t\n' '  ')" = "s1 s2 5 3 s1 s3 5 5 "
test "$(awk '$7 == "s2" {print $2}' plink2_diff.s1.sdiff | tr '\n' ' ')" = "10 30 50 "

# --variant-score: b's missing s4 is imputed as 2 * 2/6.
printf '#IID\tw\ns1\t1\ns2\t2\ns3\t0.5\ns4\t-1\n' > tmp_weights.txt
$P2 --vcf tmp_small.vcf --variant-score tmp_weights.txt --out plink2_vscore
test "$(tail -n +2 plink2_vscore.vscore | cut -f 6 | tr '\n' ' ')" = "3 2.33333 3 6.5 4.5 "

# --ref-from-fa: pos 10 has C as the reference, so REF and ALT swap; pos 20 stays.
python3 -c "
seq = ['A'] * 100
seq[9] = 'C'
seq[19] = 'G'
print('>1')
print(''.join(seq))" > tmp_ref.fa
head -n 6 tmp_small.vcf > tmp_fa.vcf
fails $P2 --vcf tmp_fa.vcf --fa tmp_ref.fa --ref-from-fa --make-pgen --out plink2_fa_err
$P2 --vcf tmp_fa.vcf --fa tmp_ref.fa --ref-from-fa force --make-pgen --out plink2_fa
test "$(grep -v '^#' plink2_fa.pvar | cut -f 2,4,5 | tr '\t\n' '  ')" = "10 C A 20 G T "
# Genotypes must flip too: C was 3/8 of the alleles at pos 10.
$P2 --pfile plink2_fa --freq --out plink2_fa_freq
test "$(tail -n +2 plink2_fa_freq.afreq | cut -f 5 | tr '\n' ' ')" = "0.625 0.333333 "

# --vscore-col-nums picks weight columns by position (the IID column is 1).
# u's values are 0/1, so its scores are the ALT dosages of s2 and s4 (with b's
# missing s4 imputed as 2/3).
printf '#IID\tw\tu\tz\ns1\t1\t0\t2\ns2\t2\t1\t0\ns3\t0.5\t0\t1\ns4\t-1\t1\t0\n' > tmp_weights3.txt
$P2 --vcf tmp_small.vcf --variant-score tmp_weights3.txt --vscore-col-nums 2,4 --out plink2_vscore24
test "$(head -n 1 plink2_vscore24.vscore | cut -f 6-)" = "$(printf 'w\tz')"
test "$(tail -n +2 plink2_vscore24.vscore | cut -f 6 | tr '\n' ' ')" = "3 2.33333 3 6.5 4.5 "
test "$(tail -n +2 plink2_vscore24.vscore | cut -f 7 | tr '\n' ' ')" = "2 2 2 5 1 "
$P2 --vcf tmp_small.vcf --variant-score tmp_weights3.txt --vscore-col-nums 3-4 --out plink2_vscore34
test "$(head -n 1 plink2_vscore34.vscore | cut -f 6-)" = "$(printf 'u\tz')"
test "$(tail -n +2 plink2_vscore34.vscore | cut -f 6 | tr '\n' ' ')" = "1 1.66667 1 2 2 "
fails $P2 --vcf tmp_small.vcf --vscore-col-nums 2 --freq --out plink2_vscore_bad

# --loop-cats: one run per category of GRP, skipping s4's NONE.  The --maf
# filter is recomputed for each category: x (s1, s3) keeps only the two 'a'
# variants, y (s2) keeps a, b, a.  The x category's .afreq matches --keep.
printf '#IID\tGRP\tQ\ns1\tx\t1\ns2\ty\t2\ns3\tx\t3\ns4\tNONE\t4\n' > tmp_grp.txt
rm -f plink2_loop.*
$P2 --vcf tmp_small.vcf --pheno tmp_grp.txt --loop-cats GRP --freq --maf 0.3 --write-snplist allow-dups --out plink2_loop
test "$(ls plink2_loop.*.afreq | tr '\n' ' ')" = "plink2_loop.x.afreq plink2_loop.y.afreq "
test "$(tr '\n' ' ' < plink2_loop.x.snplist)" = "a a "
test "$(tr '\n' ' ' < plink2_loop.y.snplist)" = "a b a "
printf 's1\ns3\n' > tmp_keep_x.txt
$P2 --vcf tmp_small.vcf --keep tmp_keep_x.txt --freq --out plink2_keep_x
diff plink2_keep_x.afreq plink2_loop.x.afreq
grep -q "Processing category 'y' (1 sample)" plink2_loop.log
fails $P2 --vcf tmp_small.vcf --pheno tmp_grp.txt --loop-cats Q --freq --out plink2_loop_bad
fails $P2 --vcf tmp_small.vcf --pheno tmp_grp.txt --loop-cats NOPE --freq --out plink2_loop_bad

# Run-control flags, on a copy of tmp_small with unique variant IDs.
sed '/^1	30	/d; /^1	50	/d' tmp_small.vcf > tmp_uniq.vcf
$P2 --vcf tmp_uniq.vcf --make-pgen --out tmp_uniq

# --error-on-freq-calc: --maf needs frequencies, unless --read-freq supplies
# them.  b (ALT freq 1/3) is filtered out; a and c (3/8 and 5/8) stay.
fails $P2 --pfile tmp_uniq --maf 0.35 --error-on-freq-calc --make-just-pvar --out plink2_efc_bad
$P2 --pfile tmp_uniq --freq --out plink2_efc_freq
$P2 --pfile tmp_uniq --read-freq plink2_efc_freq.afreq --maf 0.35 --error-on-freq-calc --make-just-pvar --out plink2_efc
test "$(grep -v '^#' plink2_efc.pvar | cut -f 3 | tr '\n' ' ')" = "a c "
$P2 --pfile tmp_uniq --error-on-freq-calc --make-just-pvar --out plink2_efc_none

# --warning-errcode: exit code 61 iff a warning was printed.
set +e
$P2 --pfile tmp_uniq --mendel --warning-errcode --out plink2_wec
wec_rc=$?
set -e
test $wec_rc -eq 61
grep -q 'Skipping --me/--mendel since there are no trios' plink2_wec.log
$P2 --pfile tmp_uniq --mendel --out plink2_wec_off
$P2 --pfile tmp_uniq --freq --warning-errcode --out plink2_wec_none

# --allow-misleading-out-arg: an --out value ending in a fileset extension
# is rejected, unless the flag says it's really meant as a prefix.
rm -f plink2_out.vcf.*
fails $P2 --pfile tmp_uniq --make-pgen --out plink2_out.vcf
test ! -e plink2_out.vcf.pgen
$P2 --pfile tmp_uniq --make-pgen --allow-misleading-out-arg --out plink2_out.vcf
test -e plink2_out.vcf.pgen
$P2 --pfile tmp_uniq --make-pgen --out plink2_out.vcf2
