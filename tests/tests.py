#! /usr/bin/env python

import subprocess
import sys

def main():
    bfile_names = ['dummy_cc1', 'dummy_cc2', 'dummy1', 'dummy2']
    bfile_names_cc = ['dummy_cc1', 'dummy_cc2'] # there must be at least 2
    bfile_names_qt = ['dummy1', 'dummy2']

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --recode --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --recode test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --recode test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.ped test2.ped', shell=True)
        if not retval == 0:
            print '--recode test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.map test2.map', shell=True)
        if not retval == 0:
            print '--recode test failed.'
            sys.exit(1)

        retval = subprocess.call('plink1 --file test1 --silent --maf 0.04999 --write-snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --file/--maf/--write-snplist test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --file test2 --silent --maf 0.04999 --write-snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --file/--maf/--write-snplist test.'
            sys.exit(1)
        subprocess.call('diff -q test1.snplist test2.snplist', shell=True)
        if not retval == 0:
            print '--file/--maf/--write-snplist test failed.'
            sys.exit(1)
    print '--file/--maf/--recode/--write-snplist test passed.'

    # --hwe takes case/control status into account, and VCF doesn't store that,
    # so this test is quantitative-trait only
    for bfn in bfile_names_qt:
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode vcf --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --recode vcf test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --hwe 0.009999 --write-snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --hwe/--vcf test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --vcf test2.vcf --silent --hwe 0.009999 --write-snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --hwe/--vcf test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.snplist test2.snplist', shell=True)
        if not retval == 0:
            print '--hwe/--recode vcf/--vcf test failed.'
            sys.exit(1)
    print '--hwe/--recode vcf/--vcf test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --recode --transpose --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --recode transpose test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode transpose --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --recode transpose test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.tped test2.tped', shell=True)
        if not retval == 0:
            print '--recode transpose test failed.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.tfam test2.tfam', shell=True)
        if not retval == 0:
            print '--recode transpose test failed.'
            sys.exit(1)

        retval = subprocess.call('plink1 --tfile test1 --silent --missing --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --missing/--tfile test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --tfile test2 --silent --missing --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --missing/--tfile test.'
            sys.exit(1)
        subprocess.call('diff -q test1.imiss test2.imiss', shell=True)
        if not retval == 0:
            print '--missing/--tfile test failed.'
            sys.exit(1)
        subprocess.call('diff -q test1.lmiss test2.lmiss', shell=True)
        if not retval == 0:
            print '--missing/--tfile test failed.'
            sys.exit(1)
    print '--missing/--recode transpose/--tfile test passed.'

    for bfn in bfile_names:
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --recode oxford --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --recode oxford test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --het --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --data/--het test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --data test2 --silent --het small-sample --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --data/--het test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.het test2.het', shell=True)
        if not retval == 0:
            print '--data/--het/--recode oxford test failed.'
            sys.exit(1)
    print '--data/--het/--recode oxford test passed.'

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

    for bfn in bfile_names_cc:
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --freq --geno 0.5399 --max-maf 0.4999 --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --freq/--geno test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --freq --geno 0.5399  --max-maf 0.4999 --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --freq/--geno test.'
            sys.exit(1)
        retval = subprocess.call('diff -q test1.frq test2.frq', shell=True)
        if not retval == 0:
            print '--freq/--geno test failed.'
            sys.exit(1)
    print '--freq/--geno test passed.'

    # don't bother with --test-mishap for now due to EM phasing algorithm
    # change

    for bfn in bfile_names_cc:
        # force at least 14 to be missing out of 513
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --geno 0.0254 --write-snplist --out test2', shell=True)
        if not retval == 0:
            print 'Unexpected error in --exclude/--hardy test.'
            sys.exit(1)
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --exclude test2.snplist --hardy --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --exclude/--hardy test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --bfile ' + bfn + ' --silent --exclude test2.snplist --hardy --out test2', shell=True)
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

    print 'All tests passed.'
    subprocess.call('rm test1.*', shell=True)
    subprocess.call('rm test2.*', shell=True)

if __name__ == '__main__':
    main()
