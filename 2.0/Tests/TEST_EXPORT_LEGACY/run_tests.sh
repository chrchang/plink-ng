#!/bin/bash

# The PLINK 1.x --recode text formats that "--export" now covers, checked
# against PLINK 1.9 output.

set -exo pipefail

BUILD=$1
EXTRA1=$2
EXTRA2=$3

plink --simulate simulate.txt --simulate-missing 0.03 --simulate-ncases 60 --simulate-ncontrols 60 --out tmp_data > /dev/null

# Whitespace is the only formatting difference, so both sides are squeezed
# before comparing.
same() {
    diff -q <(tr -s ' \t' ' ' < "$1") <(tr -s ' \t' ' ' < "$2")
}

# 1. The formats PLINK 1.9 also writes.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export lgen --out plink2_lgen
plink --bfile tmp_data --recode lgen --out plink19_lgen
same plink19_lgen.lgen plink2_lgen.lgen
same plink19_lgen.map plink2_lgen.map

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export list --out plink2_list
plink --bfile tmp_data --recode list --out plink19_list
same plink19_list.list plink2_list.list

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export rlist --out plink2_rlist
plink --bfile tmp_data --recode rlist --out plink19_rlist
same plink19_rlist.rlist plink2_rlist.rlist

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export structure --out plink2_structure
plink --bfile tmp_data --recode structure --out plink19_structure
same plink19_structure.recode.strct_in plink2_structure.strct_in

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export fastphase --out plink2_fastphase
plink --bfile tmp_data --recode fastphase --out plink19_fastphase
same plink19_fastphase.chr-1.recode.phase.inp plink2_fastphase.chr-1.phase.inp

$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export HV --out plink2_hv
plink --bfile tmp_data --recode HV --out plink19_hv
same plink19_hv.chr-1.info plink2_hv.chr-1.info
same plink19_hv.chr-1.ped plink2_hv.chr-1.ped

# 2. The '-1chr' variants drop the per-chromosome suffix.  The fixture is
#    single-chromosome, so the contents have to be the same.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export fastphase-1chr --out plink2_fp1
same plink2_fastphase.chr-1.phase.inp plink2_fp1.phase.inp
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export HV-1chr --out plink2_hv1
same plink2_hv.chr-1.ped plink2_hv1.ped
same plink2_hv.chr-1.info plink2_hv1.info

# 3. 'lgen-ref' omits the calls that match the reference allele, so it is
#    smaller, and it writes the .ref file that makes it readable back.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export lgen-ref --out plink2_lgenref
test -s plink2_lgenref.ref
test "$(wc -c < plink2_lgenref.lgen)" -lt "$(wc -c < plink2_lgen.lgen)"
plink --bfile tmp_data --recode lgen-ref --out plink19_lgenref
same plink19_lgenref.lgen plink2_lgenref.lgen
same plink19_lgenref.ref plink2_lgenref.ref

# 4. The BEAGLE and mean-genotype formats are plink2's own redesign, so they
#    are checked for shape rather than against 1.9: one line per variant, and
#    the sample count implied by the header.
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export beagle-unphased --out plink2_beagle
test "$(grep -c '^M' plink2_beagle.dat)" -eq "$(wc -l < tmp_data.bim)"
grep -q '^P FID' plink2_beagle.dat
grep -q '^I IID' plink2_beagle.dat
grep -q '^A ' plink2_beagle.dat
$BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export mgf --out plink2_mgf
test "$(wc -l < plink2_mgf.mgf)" -eq "$(wc -l < tmp_data.bim)"
test "$(wc -l < plink2_mgf.pos.txt)" -eq "$(wc -l < tmp_data.bim)"

# 5. The retired chromosome-split BEAGLE form has to say so rather than write
#    something unexpected.
if $BUILD/plink2 $EXTRA1 $EXTRA2 --bfile tmp_data --export beagle --out plink2_beagle_old 2> tmp_beagle_err.txt; then
    echo "\"--export beagle\" unexpectedly succeeded"
    exit 1
fi
grep -q "retired" tmp_beagle_err.txt
