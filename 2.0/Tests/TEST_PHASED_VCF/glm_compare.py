#!/usr/bin/env python3
"""
This compares a plink 1.x .assoc.linear or .assoc.logistic file with a plink
2.0 .glm.linear or .glm.logistic file, verifying that CHR/CHROM, BP/POS,
SNP/ID, A1, TEST, NMISS/OBS_CT, BETA or OR, STAT, and P columns match (or have
min(symmetric absolute percentage error, absolute error) within a given
tolerance) when both results aren't NA.

todo:
verify plink2 NA frequency isn't too high
"""

import argparse
import csv
import math
import sys

def parse_commandline_args():
    """
    Standard command-line parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    requiredarg = parser.add_argument_group('Required Arguments')
    requiredarg.add_argument('-1', '--plink1', type=str, required=True,
                             help="Full plink1 .assoc.linear or .assoc.logistic filename.")
    requiredarg.add_argument('-2', '--plink2', type=str, required=True,
                             help="Full plink2 .glm.linear or .glm.logistic filename.")
    requiredarg.add_argument('-t', '--tolerance', type=float, required=True,
                             help="Maximum allowed min(SMAE, absolute error).")
    cmd_args = parser.parse_args()
    return cmd_args


# See https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python .
def eprint(*args):
    print(*args, file=sys.stderr)


def float_compare_ok(val1, val2, tol):
    delta = abs(val1 - val2)
    if delta < tol:
        return True
    # SMAE = (2 * delta) / (abs(val1) + abs(val2))
    return (2 * delta) < (abs(val1) + abs(val2)) * tol


def main():
    cmd_args = parse_commandline_args()
    pc_ct = 0
    tol = cmd_args.tolerance
    with open(cmd_args.plink1, 'r') as plink1_file, open(cmd_args.plink2, 'r') as plink2_file:
        reader1 = csv.reader(plink1_file, delimiter=' ', skipinitialspace=True)
        reader2 = csv.reader(plink2_file, delimiter='\t')
        first_line1 = next(reader1)
        first_line2 = next(reader2)
        if first_line2.index('#CHROM') != 0:
            eprint('Unsupported plink2 association file format (this script requires #CHROM in front, even though plink2 can omit it).')
            sys.exit(1)
        pos_col2 = first_line2.index('POS')
        id_col2 = first_line2.index('ID')
        a1_col2 = first_line2.index('A1')
        test_col2 = first_line2.index('TEST')
        obsct_col2 = first_line2.index('OBS_CT')
        betaor_col2 = -1
        stat_col2 = -1
        is_odds_ratio = False
        try:
            betaor_col2 = first_line2.index('BETA')
        except ValueError:
            betaor_col2 = first_line2.index('OR')
            is_odds_ratio = True
        if not is_odds_ratio:
            try:
                stat_col2 = first_line2.index('T_STAT')
            except ValueError:
                stat_col2 = first_line2.index('T_OR_F_STAT')
        else:
            try:
                stat_col2 = first_line2.index('Z_STAT')
            except ValueError:
                stat_col2 = first_line2.index('Z_OR_F_STAT')
        p_col2 = first_line2.index('P')

        if (first_line1.index('CHR') != 0) or \
           (first_line1.index('BP') != 2) or \
           (first_line1.index('SNP') != 1) or \
           (first_line1.index('A1') != 3) or \
           (first_line1.index('TEST') != 4) or \
           (first_line1.index('NMISS') != 5) or \
           (first_line1.index('STAT') != 7) or \
           (first_line1.index('P') != 8):
            eprint('Unexpected plink1 association file format.')
            sys.exit(1)

        for row1 in reader1:
            row2 = next(reader2)
            if row1[0] != row2[0] or \
               row1[2] != row2[pos_col2] or \
               row1[1] != row2[id_col2] or \
               row1[3] != row2[a1_col2] or \
               row1[4] != row2[test_col2] or \
               row1[5] != row2[obsct_col2]:
                eprint('Header column mismatch between association files.')
                sys.exit(1)
            if row1[6] != 'NA' and row2[betaor_col2] != 'NA':
                val1 = float(row1[6])
                val2 = float(row2[betaor_col2])
                if is_odds_ratio:
                    # this is a more appropriate scale
                    val1 = math.log(val1)
                    val2 = math.log(val2)
                if not float_compare_ok(val1, val2, tol):
                    eprint('BETA/OR mismatch.')
                    sys.exit(1)
            if row1[7] != 'NA' and row1[7] != 'inf' and row2[stat_col2] != 'NA':
                val1 = float(row1[7])
                val2 = float(row2[stat_col2])
                if row1[4][-2:] != 'DF':
                    # skip GENO_2DF, USER_xDF due to transition to F-statistics
                    if not float_compare_ok(val1, val2, tol):
                        eprint('STAT mismatch.')
                        sys.exit(1)
            if row1[8] != 'NA' and row2[p_col2] != 'NA':
                val1 = math.log(float(row1[8]))
                val2 = math.log(float(row2[p_col2]))
                if not float_compare_ok(val1, val2, tol):
                    eprint('P-value mismatch.')
                    sys.exit(1)


if __name__ == '__main__':
    main()
