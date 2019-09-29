#!/usr/bin/env python3
"""
This compares two sets of {.eigenval, .eigenvec, .eigenvec.var} files, and
verifies that symmetric absolute percentage error between the two sets is less
than the given tolerance.  (To avoid spurious test failures due to isolated
near-zero values, we only compute a single error statistic on sums of absolute
values, instead of taking means across all samples/PCs.)
"""

import argparse
import csv
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
def eprint(*args):
    print(*args, file=sys.stderr)


def main():
    cmd_args = parse_commandline_args()
    pc_ct = 0
    tol = cmd_args.tolerance
    with open(cmd_args.plink1 + '.eigenval', 'r') as eigval1_file, open(cmd_args.plink2 + '.eigenval', 'r') as eigval2_file:
        absdiff_sum = 0.0
        avgabs_x2_sum = 0.0
        for line1 in eigval1_file:
            line1.rstrip('\n')
            eigval1 = float(line1)
            line2 = eigval2_file.readline()
            line2.rstrip('\n')
            eigval2 = float(line2)
            absdiff_sum += abs(eigval2 - eigval1)
            avgabs_x2_sum += abs(eigval2) + abs(eigval1)
            pc_ct += 1
        if absdiff_sum > tol * avgabs_x2_sum * 0.5:
            eprint('Eigenvalue mismatch.')
            sys.exit(1)
    with open(cmd_args.plink1 + '.eigenvec', 'r') as eigvec1_file, open(cmd_args.plink2 + '.eigenvec', 'r') as eigvec2_file:
        reader1 = csv.reader(eigvec1_file, delimiter='\t')
        reader2 = csv.reader(eigvec2_file, delimiter='\t')
        next(reader1)
        next(reader2)
        absdiff1_sums = [0.0] * pc_ct
        absdiff2_sums = [0.0] * pc_ct
        avgabs_x2_sums = [0.0] * pc_ct
        for row1 in reader1:
            row2 = next(reader2)
            if row1[0] != row2[0] or row1[1] != row2[1]:
                eprint('Sample ID mismatch between .eigenvec files.')
                sys.exit(1)
            row1pc = row1[2:]
            row2pc = row2[2:]
            for pc_idx in range(pc_ct):
                val1 = float(row1pc[pc_idx])
                val2 = float(row2pc[pc_idx])
                absdiff1_sums[pc_idx] += abs(val2 - val1)
                # all signs may be flipped
                absdiff2_sums[pc_idx] += abs(val2 + val1)
                avgabs_x2_sums[pc_idx] += abs(val2) + abs(val1)
        for pc_idx in range(pc_ct):
            if min(absdiff1_sums[pc_idx], absdiff2_sums[pc_idx]) > tol * avgabs_x2_sums[pc_idx] * 0.5:
                eprint('Eigenvector mismatch.')
                sys.exit(1)
    with open(cmd_args.plink1 + '.eigenvec.var', 'r') as varwt1_file, open(cmd_args.plink2 + '.eigenvec.var', 'r') as varwt2_file:
        reader1 = csv.reader(varwt1_file, delimiter='\t')
        reader2 = csv.reader(varwt2_file, delimiter='\t')
        next(reader1)
        next(reader2)
        absdiff1_sums = [0.0] * pc_ct
        absdiff2_sums = [0.0] * pc_ct
        avgabs_x2_sums = [0.0] * pc_ct
        for row1 in reader1:
            row2 = next(reader2)
            if row1[0] != row2[0]:
                eprint('Chromosome mismatch between .eigenvec.var files.')
                sys.exit(1)
            if row1[1] != row2[1]:
                eprint('Variant ID mismatch between .eigenvec.var files.')
                sys.exit(1)
            if not (((row1[2] == row2[2]) and (row1[3] == row2[3])) or
                    ((row1[2] == row2[3]) and (row1[3] == row2[2]))):
                eprint('Allele mismatch between .eigenvec.var files.')
                sys.exit(1)
            row1pc = row1[4:]
            row2pc = row2[4:]
            for pc_idx in range(pc_ct):
                val1 = float(row1pc[pc_idx])
                val2 = float(row2pc[pc_idx])
                absdiff1_sums[pc_idx] += abs(val2 - val1)
                # all signs may be flipped
                absdiff2_sums[pc_idx] += abs(val2 + val1)
                avgabs_x2_sums[pc_idx] += abs(val2) + abs(val1)
        for pc_idx in range(pc_ct):
            if min(absdiff1_sums[pc_idx], absdiff2_sums[pc_idx]) > tol * avgabs_x2_sums[pc_idx] * 0.5:
                eprint('Variant weight mismatch.')
                sys.exit(1)


if __name__ == '__main__':
    main()
