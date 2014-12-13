#! /usr/bin/env python

import subprocess
import sys

def main():
    bfile_names = ['dummy_cc1', 'dummy_cc2', 'dummy1', 'dummy2']
    bfile_names_cc = ['dummy_cc1', 'dummy_cc2']
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

        retval = subprocess.call('plink1 --file test1 --silent --maf 0.05 --write-snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --file/--maf/--write-snplist test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --file test2 --silent --maf 0.05 --write-snplist --out test2', shell=True)
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
        retval = subprocess.call('plink1 --bfile ' + bfn + ' --silent --hwe 0.01 --write-snplist --out test1', shell=True)
        if not retval == 0:
            print 'Unexpected error in --hwe/--vcf test.'
            sys.exit(1)
        retval = subprocess.call('plink2 --vcf test2.vcf --silent --hwe 0.01 --write-snplist --out test2', shell=True)
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

    print 'All tests passed.'
    subprocess.call('rm test1.*', shell=True)
    subprocess.call('rm test2.*', shell=True)

if __name__ == '__main__':
    main()
