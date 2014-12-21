#! /usr/bin/env python

import subprocess
import sys

def main():
    bfile_names = ['dummy_cc1', 'dummy_cc2', 'dummy1', 'dummy2', 'trio']
    bfile_names_cc = ['dummy_cc1', 'dummy_cc2'] # there must be at least 2
    bfile_names_qt = ['dummy1', 'dummy2']
    bfile_names_fam = ['trio']

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --recode --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --file/--maf/--nonfounders/--recode/--write-snplist test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --file/--maf/--nonfounders/--recode/--write-snplist test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.ped test2.ped', shell=True)
        if not retval == 0:
            print '--file/--maf/--nonfounders/--recode/--write-snplist test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.map test2.map', shell=True)
        if not retval == 0:
            print '--file/--maf/--nonfounders/--recode/--write-snplist test failed.'
            sys.exit(1)

        retval = subprocess.call('plink1 --file test1 --silent --nonfounders --maf 0.04999 --write-snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --file/--maf/--nonfounders/--recode/--write-snplist test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --file test2 --silent --nonfounders --maf 0.04999 --write-snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --file/--maf/--nonfounders/--recode/--write-snplist test.'
            sys.exit(1)
        subprocess.call('diff -q test1.snplist test2.snplist', shell=True)
        if not retval == 0:
            print '--file/--maf/--nonfounders/--recode/--write-snplist test failed.'
            sys.exit(1)
    print '--file/--maf/--nonfounders/--recode/--write-snplist test passed.'

    # --hwe takes case/control status into account, and VCF doesn't store that,
    # so this test is quantitative-trait only
    for bfn in bfile_names_qt:
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode vcf --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --hwe/--recode vcf/--vcf test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --hwe 0.009999 --write-snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --hwe/--recode vcf/--vcf test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --vcf test2.vcf --silent --hwe 0.009999 --write-snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --hwe/--recode vcf/--vcf test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.snplist test2.snplist', shell=True)
        if not retval == 0:
            print '--hwe/--recode vcf/--vcf test failed.'
            sys.exit(1)
    print '--hwe/--recode vcf/--vcf test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --recode --transpose --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --missing/--recode transpose/--tfile test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode transpose --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --missing/--recode transpose/--tfile test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.tped test2.tped', shell=True)
        if not retval == 0:
            print '--missing/--recode transpose/--tfile test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.tfam test2.tfam', shell=True)
        if not retval == 0:
            print '--missing/--recode transpose/--tfile test failed.'
            sys.exit(1)

        retval = subprocess.call('plink1 --tfile test1 --silent --missing --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --missing/--recode transpose/--tfile test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --tfile test2 --silent --missing --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --missing/--recode transpose/--tfile test.'
            sys.exit(1)
        subprocess.call('diff -q test1.imiss test2.imiss', shell=True)
        if not retval == 0:
            print '--missing/--recode transpose/--tfile test failed.'
            sys.exit(1)
        subprocess.call('diff -q test1.lmiss test2.lmiss', shell=True)
        if not retval == 0:
            print '--missing/--recode transpose/--tfile test failed.'
            sys.exit(1)
    print '--missing/--recode transpose/--tfile test passed.'

    retval = subprocess.call('plink1 --bfile ' + bfile_names_cc[0] + ' --silent --bmerge ' + bfile_names_cc[1] + '.bed ' + bfile_names_cc[1] + '.bim ' + bfile_names_cc[1] + '.fam --max-maf 0.4999 --make-bed --out test1', shell=True)
    if not retval == 0:
        print 'Unexpected error in --bmerge/--make-bed/--max-maf test.'
        sys.exit(1)
    retval = subprocess.call('plink2 --bfile ' + bfile_names_cc[0] + ' --silent --bmerge ' + bfile_names_cc[1] + ' --max-maf 0.4999 --make-bed --out test2', shell=True)
    if not retval == 0:
        print 'Unexpected error in --bmerge/--make-bed/--max-maf test.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bed test2.bed', shell=True)
    if not retval == 0:
        print '--bmerge/--make-bed/--max-maf test failed.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
    if not retval == 0:
        print '--bmerge/--make-bed/--max-maf test failed.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
    if not retval == 0:
        print '--bmerge/--make-bed/--max-maf test failed.'
        sys.exit(1)
    print '--bmerge/--make-bed/--max-maf test passed.'

    subprocess.call('rm test2.bim', shell=True)
    retval = subprocess.call('plink2 --bfile ' + bfile_names_cc[0] + ' --silent --bmerge ' + bfile_names_cc[1] + ' --max-maf 0.4999 --make-just-bim --out test2', shell=True)
    if not retval == 0:
        print 'Unexpected error in --make-just-bim test.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
    if not retval == 0:
        print '--make-just-bim test failed.'
        sys.exit(1)
    print '--make-just-bim test passed.'

    subprocess.call('rm test2.fam', shell=True)
    retval = subprocess.call('plink2 --bfile ' + bfile_names_cc[0] + ' --silent --bmerge ' + bfile_names_cc[1] + ' --make-just-fam --out test2', shell=True)
    if not retval == 0:
        print 'Unexpected error in --make-just-bim test.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.fam test2.fam', shell=True)
    if not retval == 0:
        print '--make-just-fam test failed.'
        sys.exit(1)
    print '--make-just-fam test passed.'

    for bfn in bfile_names_cc:
        # force at least 14 to be missing out of 513...
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --geno 0.06049 --max-maf 0.4999 --write-snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --extract/--freq/--geno test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --freq --extract test2.snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --extract/--freq/--geno test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --freq --extract test2.snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --extract/--freq/--geno test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.frq test2.frq', shell=True)
        if not retval == 0:
            print '--extract/--freq/--geno test failed.'
            sys.exit(1)
    print '--extract/--freq/--geno test passed.'

    for bfn in bfile_names_fam:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --covar ' + bfn + '.fam --covar-number 3 --filter-founders --write-covar --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --covar/--covar-number/--filter-founders/--write-covar test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --covar ' + bfn + '.fam --covar-number 3 --filter-founders --write-covar --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --covar/--covar-number/--filter-founders/--write-covar test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.cov test2.cov', shell=True)
        if not retval == 0:
            print '--covar/--covar-number/--filter-founders/--write-covar test failed.'
            sys.exit(1)
    print '--covar/--covar-number/--filter-founders/--write-covar test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --cluster --K 2 --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --cluster/--filter-females/--within/--write-cluster test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --cluster only2 old-tiebreaks --K 2 --out ' + bfn, shell=True)
        if not retval == 0:
            print 'Unexpected error in --cluster/--filter-females/--within/--write-cluster test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.cluster2 ' + bfn + '.cluster2', shell=True)
        if not retval == 0:
            print '--cluster/--filter-females/--within/--write-cluster test failed.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --within ' + bfn + '.cluster2 --filter-females --write-cluster --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --cluster/--filter-females/--within/--write-cluster test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --within ' + bfn + '.cluster2 --filter-females --write-cluster --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --cluster/--filter-females/--within/--write-cluster test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.clst test2.clst', shell=True)
        if not retval == 0:
            print '--cluster/--filter-females/--within/--write-cluster test failed.'
            sys.exit(1)
    print '--cluster/--filter-females/--within/--write-cluster test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --make-set set.txt --write-set --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --make-set/--set/--set-table/--write-set test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --make-set set.txt --write-set --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --make-set/--set/--set-table/--write-set test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.set test2.set', shell=True)
        if not retval == 0:
            print '--make-set/--set/--set-table/--write-set test failed.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --make-set set.txt --set-table --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --make-set/--set/--set-table/--write-set test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --set test2.set --set-table --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --make-set/--set/--set-table/--write-set test.'
            sys.exit(1)
        # PLINK 1.07 generates a .set.table file with a double-tab after the
        # third column, which we deliberately don't replicate.  Remove it and
        # the header row (which doesn't have the double-tab) before diffing.
        retval = subprocess.call('cat test1.set.table | tail -n +2 | cut -f 1-3,5- > test1.set.table2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --make-set/--set/--set-table/--write-set test.'
            sys.exit(1)
        retval = subprocess.call('cat test2.set.table | tail -n +2 > test2.set.table2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --make-set/--set/--set-table/--write-set test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.set.table2 test2.set.table2', shell=True)
        if not retval == 0:
            print '--make-set/--set/--set-table/--write-set test failed.'
            sys.exit(1)
    print '--make-set/--set/--set-table/--write-set test passed.'
 
    retval = subprocess.call('plink2 --bfile ' + bfile_names_cc[1] + ' --silent --recode --out test2', shell=True)
    if not retval == 0:
        print 'Unexpected error in --merge test.'
        sys.exit(1)
    retval = subprocess.call('plink2 --bfile ' + bfile_names_cc[0] + ' --silent --merge test2 --max-maf 0.4999 --make-bed --out test2', shell=True)
    if not retval == 0:
        print 'Unexpected error in --merge test.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bed test2.bed', shell=True)
    if not retval == 0:
        print '--merge test failed.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
    if not retval == 0:
        print '--merge test failed.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
    if not retval == 0:
        print '--merge test failed.'
        sys.exit(1)
    print '--merge test passed.'

    subprocess.call('rm merge_list.txt', shell=True)
    retval = subprocess.call('echo ' + bfile_names_cc[0] + '.bed ' + bfile_names_cc[0] + '.bim ' + bfile_names_cc[0] + '.fam > merge_list.txt', shell=True)
    if not retval == 0:
        print 'Unexpected error in --merge-list test.'
        sys.exit(1)
    retval = subprocess.call('echo ' + bfile_names_cc[1] + '.bed ' + bfile_names_cc[1] + '.bim ' + bfile_names_cc[1] + '.fam >> merge_list.txt', shell=True)
    if not retval == 0:
        print 'Unexpected error in --merge-list test.'
        sys.exit(1)
    retval = subprocess.call('plink2 --merge-list merge_list.txt --silent --max-maf 0.4999 --make-bed --out test2', shell=True)
    if not retval == 0:
        print 'Unexpected error in --merge-list test.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bed test2.bed', shell=True)
    if not retval == 0:
        print '--merge-list test failed.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
    if not retval == 0:
        print '--merge-list test failed.'
        sys.exit(1)
    retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
    if not retval == 0:
        print '--merge-list test failed.'
        sys.exit(1)
    subprocess.call('rm merge_list.txt', shell=True)
    print '--merge-list test passed.'

    for bfn in bfile_names_cc:
        # diff error likely at MAF = 0.4875 due to rounding
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --flip-scan --mind 0.05399 --maf 0.4876 --max-maf 0.4999 --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --flip-scan/--mind test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --flip-scan --mind 0.05399 --maf 0.4876 --max-maf 0.4999 --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --flip-scan/--mind test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.flipscan test2.flipscan', shell=True)
        if not retval == 0:
            print '--flip-scan/--mind test failed.'
            sys.exit(1)
    print '--flip-scan/--mind test passed.'


    # don't bother with --test-mishap for now due to EM phasing algorithm
    # change

    for bfn in bfile_names_cc:
        # force at least 14 to be missing out of 513...
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --geno 0.0254 --write-snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --exclude/--hardy test.'
            sys.exit(1)
        # ...but no more than 31 out of 512.

        # --hardy + --geno order of operations has changed, for good reason.
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --exclude test1.snplist --geno 0.06049 --write-snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --exclude/--hardy test.'
            sys.exit(1)

        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --extract test2.snplist --hardy --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --exclude/--hardy test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --extract test2.snplist --hardy --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --exclude/--hardy test.'
            sys.exit(1)
        # skip column 8 due to likelihood of floating point error
        # note that leading space(s) cause it to be column 9 after tr
        retval = subprocess.call("cat test1.hwe | tr -s ' ' '\t' | cut -f 1-8,10 > test1.hwe2", shell=True)
        retval = subprocess.call("cat test2.hwe | tr -s ' ' '\t' | cut -f 1-8,10 > test2.hwe2", shell=True)
        retval = subprocess.call('diff -q test1.hwe2 test2.hwe2', shell=True)
        if not retval == 0:
            print '--exclude/--hardy test failed.'
            sys.exit(1)
    print '--exclude/--hardy test passed.'

    for bfn in bfile_names_fam:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --mendel --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --mendel test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --mendel --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --mendel test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.mendel test2.mendel', shell=True)
        if not retval == 0:
            print '--mendel test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.imendel test2.imendel', shell=True)
        if not retval == 0:
            print '--mendel test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.fmendel test2.fmendel', shell=True)
        if not retval == 0:
            print '--mendel test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.lmendel test2.lmendel', shell=True)
        if not retval == 0:
            print '--mendel test failed.'
            sys.exit(1)
    print '--mendel test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode oxford --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --data/--het/--recode oxford test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --het --nonfounders --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --data/--het/--recode oxford test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --data test2 --silent --het --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --data/--het/--recode oxford test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.het test2.het', shell=True)
        if not retval == 0:
            print '--data/--het/--recode oxford test failed.'
            sys.exit(1)
    print '--data/--het/--recode oxford test passed.'

    for bfn in bfile_names_cc:
        retval = subprocess.call('rm test1.bim.tmp', shell=True)
        retval = subprocess.call('cat ' + bfn + ".bim | sed 's/^1/23/' > test1.bim.tmp", shell=True)
        if not retval == 0:
            print 'Unexpected error in --check-sex/--impute-sex/--read-freq test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bed ' + bfn + '.bed --bim test1.bim.tmp --fam ' + bfn + '.fam --silent --freq --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --check-sex/--impute-sex/--read-freq test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bed ' + bfn + '.bed --bim test1.bim.tmp --fam ' + bfn + '.fam --silent --read-freq test1.frq --check-sex --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --check-sex/--impute-sex/--read-freq test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bed ' + bfn + '.bed --bim test1.bim.tmp --fam ' + bfn + '.fam --silent --read-freq test1.frq --check-sex --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --check-sex/--impute-sex/--read-freq test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.sexcheck test2.sexcheck', shell=True)
        if not retval == 0:
            print '--check-sex/--impute-sex/--read-freq test failed.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bed ' + bfn + '.bed --bim test1.bim.tmp --fam ' + bfn + '.fam --silent --read-freq test1.frq --maf 0.4876 --max-maf 0.4999 --impute-sex --make-bed --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --check-sex/--impute-sex/--read-freq test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bed ' + bfn + '.bed --bim test1.bim.tmp --fam ' + bfn + '.fam --silent --read-freq test1.frq --maf 0.4876 --max-maf 0.4999 --impute-sex --make-bed --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --check-sex/--impute-sex/--read-freq test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.sexcheck test2.sexcheck', shell=True)
        if not retval == 0:
            print '--check-sex/--impute-sex/--read-freq test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.bed test2.bed', shell=True)
        if not retval == 0:
            print '--check-sex/--impute-sex/--read-freq test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.bim test2.bim', shell=True)
        if not retval == 0:
            print '--check-sex/--impute-sex/--read-freq test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.fam test2.fam', shell=True)
        if not retval == 0:
            print '--check-sex/--impute-sex/--read-freq test failed.'
            sys.exit(1)
    print '--check-sex/--impute-sex/--read-freq test passed.'

    # Skip --indep-pairwise for now since harmless minor differences are
    # expected.

    # floating-point error is more likely to affect --r, so we just test --r2
    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --r2 --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --r2 test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --r2 --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --r2 test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.ld test2.ld', shell=True)
        if not retval == 0:
            print '--r2 test failed.'
            sys.exit(1)
    print '--r2 test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --make-set set.txt --write-set --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --gene/--show-tags test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --set test1.set --gene set2 --show-tags all --tag-r2 0.031 --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --gene/--show-tags test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --make-set set.txt --gene set2 --show-tags all --tag-r2 0.031 --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --gene/--show-tags test.'
            sys.exit(1)
        # ignore column 7 since we changed the interval length definition
        retval = subprocess.call("cat test1.tags.list | sed 's/^[[:space:]]*//g' | tr -s ' ' '\t' | cut -f 1-6,8 > test1.tags.list2", shell=True)
        if not retval == 0:
            print 'Unexpected error in --gene/--show-tags test.'
            sys.exit(1)
        retval = subprocess.call("cat test2.tags.list | sed 's/^[[:space:]]*//g' | tr -s ' ' '\t' | cut -f 1-6,8 > test2.tags.list2", shell=True)
        if not retval == 0:
            print 'Unexpected error in --gene/--show-tags test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.tags.list2 test2.tags.list2', shell=True)
        if not retval == 0:
            print '--gene/--show-tags test failed.'
            sys.exit(1)
    print '--gene/--show-tags test passed.'

    # skip --blocks for now since test files lack the necessary LD

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --cluster --distance-matrix --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --distance-matrix test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --distance-matrix --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --distance-matrix test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.mdist test2.mdist', shell=True)
        if not retval == 0:
            print '--distance-matrix test failed.'
            sys.exit(1)
    print '--distance-matrix test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --genome --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --genome test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --genome --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --genome test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.genome test2.genome', shell=True)
        if not retval == 0:
            print '--genome test failed.'
            sys.exit(1)
    print '--genome test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --homozyg --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --homozyg test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --homozyg subtract-1-from-lengths --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --homozyg test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.hom test2.hom', shell=True)
        if not retval == 0:
            print '--homozyg test failed.'
            sys.exit(1)
    print '--homozyg test passed.'

    # skip --neighbour test for now due to likelihood of ties.  (Caught a
    # command-line parsing bug when trying to write that test, though.)

    print 'All tests passed.'
    subprocess.call('rm test1.*', shell=True)
    subprocess.call('rm test2.*', shell=True)


if __name__ == '__main__':
    main()
