#!/usr/bin/python
"""
This compares two sets of {.eigenval, .eigenvec, .eigenvec.var} files, and
verifies that symmetric absolute percentage error between the two sets is less
than the given tolerance.  (To avoid problems spurious test failures due to
isolated near-zero values, we only compute a single error statistic on sums of
absolute values, instead of taking means across all samples/PCs.)
"""

from __future__ import print_function
import argparse
import sys

def parse_commandline_args():
    """
    Standard command-line parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    requiredarg = parser.add_argument_group('Required Arguments')
    requiredarg.add_argument('-1', '--plink1', type=str, required=True,
                             help="Filename prefix for plink1 PCA values.")
    requiredarg.add_argument('-2', '--plink2', type=str, required=True,
                             help="Filename prefix for plink2 PCA values to validate.")
    requiredarg.add_argument('-t', '--tolerance', type=float, required=True,
                             help="Maximum allowed SMAE.")
    cmd_args = parser.parse_args()
    return cmd_args


# See https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python .
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def main():
    cmd_args = parse_commandline_args()
    pc_ct = 0
    with open(cmd_args.plink1 + '.eigenval', 'r') as eigval1_file, open(cmd_args.plink2 + '.eigenval', 'r') as eigval2_file:
        sum_absdiff = 0.0
        sum_avgabs_x2 = 0.0
        for line1 in eigval1_file:
            line1.rstrip('\n')
            eigval1 = float(line1)
            line2 = eigval2_file.readline()
            line2.rstrip('\n')
            eigval2 = float(line2)
            sum_absdiff += abs(eigval2 - eigval1)
            sum_avgabs_x2 += abs(eigval2) + abs(eigval1)
            pc_ct += 1
        if sum_absdiff > cmd_args.tolerance * sum_avgabs_x2 * 0.5:
            eprint('Eigenvalue mismatch.')
            sys.exit(1)
    # todo: fill in other two comparison functions
    with open(cmd_args.plink1 + '.eigenvec', 'r') as eigvec1_file, open(cmd_args.plink2 + '.eigenvec', 'r') as eigvec2_file:
        pass
    with open(cmd_args.plink1 + '.eigenvec.var', 'r') as varwt1_file, open(cmd_args.plink2 + '.eigenvec.var', 'r') as varwt2_file:
        pass


if __name__ == '__main__':
    main()
