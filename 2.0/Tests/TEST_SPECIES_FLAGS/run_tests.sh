#!/bin/bash

# Chromosome-set flags with no other test coverage: --cow, --dog, --horse,
# --mouse, --rice, --sheep, --human, --autosome-num, and --chr-override.
#
# Also a regression test for ##chrSet header parsing: 'PAR2' used to be
# stored in the PAR1 slot, so after reading "...,PAR1,PAR2>" (which plink2
# writes itself, including into the temporary .pvar of a VCF import), PAR2
# was undefined and PAR1 pointed at the PAR2 code.  PAR2 variants were then
# rejected, PAR1 variants crashed plink2, and PAR2 vanished from rewritten
# ##chrSet lines.

set -exo pipefail

plink2="$1/plink2 $2 $3"

# $1: output VCF; remaining args: chromosome codes, one het variant each.
make_vcf() {
    local out=$1
    shift
    {
        printf '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM\tF\n'
        for c in "$@"; do
            printf "$c\t100\tc$c\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n"
        done
    } > $out
}
printf '#IID\tSEX\nM\t1\nF\t2\n' > tmp.psam

# $1: VCF; remaining args: plink2 flags.  The import must fail.
check_import_fails() {
    local vcf=$1
    shift
    if $plink2 --vcf $vcf --psam tmp.psam "$@" --make-pgen --out tmp_fail > /dev/null 2>&1; then
        exit 1
    fi
}

# $1: species flag; $2: expected ##chrSet line; $3: expected "CHROM:
# HOM_REF HET HOM_ALT MISSING" rows, '|'-joined; remaining args: chromosome
# codes.  Both samples are het everywhere, so male X/Y calls become missing,
# and females are left out of chrY counts.
check_species() {
    local flag=$1
    local chrset=$2
    local counts=$3
    shift 3
    make_vcf tmp_$flag.vcf "$@"
    $plink2 --vcf tmp_$flag.vcf --psam tmp.psam --$flag --make-pgen --out tmp_$flag > /dev/null
    test "$(grep '^##chrSet' tmp_$flag.pvar)" = "$chrset"
    $plink2 --pfile tmp_$flag --geno-counts --out tmp_$flag > /dev/null
    test "$(tail -n +2 tmp_$flag.gcount | cut -f 1,5,6,7,10 | tr '\n' '|')" = "$(printf "$counts")|"
}

auto='0\t2\t0\t0'
xrow='0\t1\t0\t1'
yrow='0\t0\t0\t1'
# Codes: 1, last autosome, then X, Y, (XY,) MT as numbers.
check_species cow '##chrSet=<ID=1,autosomePairCt=29,X,Y,M>' "1\t$auto|29\t$auto|X\t$xrow|Y\t$yrow|MT\t$auto" 1 29 30 31 33
check_species dog '##chrSet=<ID=1,autosomePairCt=38,X,Y,XY,M,PAR1,PAR2>' "1\t$auto|38\t$auto|X\t$xrow|Y\t$yrow|XY\t$auto|MT\t$auto" 1 38 39 40 41 42
check_species horse '##chrSet=<ID=1,autosomePairCt=31,X,Y>' "1\t$auto|31\t$auto|X\t$xrow|Y\t$yrow" 1 31 32 33
check_species mouse '##chrSet=<ID=1,autosomePairCt=19,X,Y>' "1\t$auto|19\t$auto|X\t$xrow|Y\t$yrow" 1 19 20 21
check_species sheep '##chrSet=<ID=1,autosomePairCt=26,X,Y>' "1\t$auto|26\t$auto|X\t$xrow|Y\t$yrow" 1 26 27 28
# Rice chromosomes are all haploid, so every het call is missing.
check_species rice '##chrSet=<ID=1,haploidAutosomeCt=12>' '1\t0\t0\t0\t2|12\t0\t0\t0\t2' 1 12
# Named codes work too.
make_vcf tmp_cow_named.vcf 1 X Y MT
$plink2 --vcf tmp_cow_named.vcf --psam tmp.psam --cow --make-pgen --out tmp_cow_named > /dev/null
test "$(grep -v '^#' tmp_cow_named.pvar | cut -f 1 | tr '\n' ' ')" = "1 X Y MT "

# Codes outside each set are rejected: cow has no XY (32), nothing past MT,
# and the human default stops at 26.
make_vcf tmp_bad32.vcf 1 32
make_vcf tmp_bad34.vcf 1 34
make_vcf tmp_bad13.vcf 1 13
check_import_fails tmp_bad32.vcf --cow
check_import_fails tmp_bad34.vcf --cow
check_import_fails tmp_cow.vcf
check_import_fails tmp_bad13.vcf --rice
check_import_fails tmp_cow.vcf --cow --dog
check_import_fails tmp_cow.vcf --cow --chr-set 29

# --output-chr renders the species' codes.
$plink2 --pfile tmp_cow --output-chr 26 --make-just-pvar --out tmp_cow26 > /dev/null
test "$(grep -v '^#' tmp_cow26.pvar | cut -f 1 | tr '\n' ' ')" = "1 29 30 31 33 "
$plink2 --pfile tmp_cow --output-chr chrM --make-just-pvar --out tmp_cowchr > /dev/null
test "$(grep -v '^#' tmp_cowchr.pvar | cut -f 1 | tr '\n' ' ')" = "chr1 chr29 chrX chrY chrM "

# --autosome-num n: n autosomes plus X; no Y/XY/MT.
make_vcf tmp_an.vcf 1 5 6
$plink2 --vcf tmp_an.vcf --psam tmp.psam --autosome-num 5 --make-pgen --out tmp_an > /dev/null
test "$(grep '^##chrSet' tmp_an.pvar)" = '##chrSet=<ID=1,autosomePairCt=5,X>'
test "$(grep -v '^#' tmp_an.pvar | cut -f 1 | tr '\n' ' ')" = "1 5 X "
$plink2 --pfile tmp_an --geno-counts --out tmp_an > /dev/null
test "$(awk -F '\t' '$1 == "X"' tmp_an.gcount | cut -f 5,6,7,10)" = "$(printf "$xrow")"
make_vcf tmp_an7.vcf 1 7
check_import_fails tmp_an7.vcf --autosome-num 5

# --human: only it adds ##chrSet for the human set.
make_vcf tmp_h.vcf 1 22
$plink2 --vcf tmp_h.vcf --human --make-pgen --out tmp_h > /dev/null
test "$(grep '^##chrSet' tmp_h.pvar)" = '##chrSet=<ID=1,autosomePairCt=22,X,Y,XY,M,PAR1,PAR2>'
$plink2 --vcf tmp_h.vcf --make-pgen --out tmp_h0 > /dev/null
test "$(grep -c '^##chrSet' tmp_h0.pvar || true)" -eq 0

# A file ##chrSet line conflicting with the command line errors out, unless
# --chr-override says which one wins.
if $plink2 --pfile tmp_cow --dog --make-just-pvar --out tmp_conflict > /dev/null 2>&1; then
    exit 1
fi
$plink2 --pfile tmp_cow --dog --chr-override --make-just-pvar --out tmp_co > /dev/null
test "$(grep '^##chrSet' tmp_co.pvar)" = '##chrSet=<ID=1,autosomePairCt=38,X,Y,XY,M,PAR1,PAR2>'
# Under the dog set, the cow's 29 is an autosome and X/Y/MT keep their names.
test "$(grep -v '^#' tmp_co.pvar | cut -f 1 | tr '\n' ' ')" = "1 29 X Y MT "
$plink2 --pfile tmp_cow --dog --chr-override file --make-just-pvar --out tmp_cof > /dev/null
test "$(grep '^##chrSet' tmp_cof.pvar)" = '##chrSet=<ID=1,autosomePairCt=29,X,Y,M>'
# Same set on both sides: no conflict, no override needed.
$plink2 --pfile tmp_cow --cow --make-just-pvar --out tmp_cowcow > /dev/null

# --- ##chrSet PAR1/PAR2 parsing regression ---
chrset_human='##chrSet=<ID=1,autosomePairCt=22,X,Y,XY,M,PAR1,PAR2>'
printf "$chrset_human\n#CHROM\tPOS\tID\tREF\tALT\n1\t100\ta\tA\tG\nPAR1\t200\tp1\tA\tG\nX\t300000000\tx\tA\tG\nPAR2\t400\tp2\tA\tG\n" > tmp_par.pvar
$plink2 --pvar tmp_par.pvar --make-just-pvar --out tmp_par_out > /dev/null
test "$(grep -v '^#' tmp_par_out.pvar | cut -f 1,3 | tr '\n' ' ')" = "$(printf '1\ta PAR1\tp1 X\tx PAR2\tp2 ')"
test "$(grep '^##chrSet' tmp_par_out.pvar)" = "$chrset_human"
# PAR1 only (this used to crash).
grep -v '^PAR2' tmp_par.pvar > tmp_par1.pvar
$plink2 --pvar tmp_par1.pvar --make-just-pvar --out tmp_par1_out > /dev/null
test "$(grep -v '^#' tmp_par1_out.pvar | cut -f 1 | tr '\n' ' ')" = "1 PAR1 X "
# VCF import goes through a temporary .pvar with a ##chrSet line.
for flag in human dog; do
    make_vcf tmp_par_$flag.vcf 1 PAR2
    $plink2 --vcf tmp_par_$flag.vcf --$flag --make-pgen --out tmp_par_$flag > /dev/null
    test "$(grep -v '^#' tmp_par_$flag.pvar | cut -f 1 | tr '\n' ' ')" = "1 PAR2 "
    grep -q '^##chrSet=<.*,PAR1,PAR2>$' tmp_par_$flag.pvar
done
