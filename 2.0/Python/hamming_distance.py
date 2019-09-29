#!/usr/bin/env python3
import pgenlib
import numpy as np
import sys

def main():
    arg_ct = len(sys.argv)
    if arg_ct < 4:
        print("Usage: python hamming_distance.py <.bed/.pgen> <sample index 1> <sample index 2> [raw_sample_ct=<val>]")
        print("* raw_sample_ct is required for .bed files.")
        return
    specified_sample_ct = None
    sample_ct = 2
    sample_subset = np.empty(sample_ct, np.uint32)
    sample_subset[0] = int(sys.argv[2])
    sample_subset[1] = int(sys.argv[3])
    if arg_ct > 4:
        if not sys.argv[4].startswith("raw_sample_ct="):
            print("Error: Invalid raw_sample_ct parameter.")
            return
        specified_sample_ct = int(sys.argv[4][14:])
    with pgenlib.PgenReader(bytes(sys.argv[1], 'utf8'), raw_sample_ct = specified_sample_ct, sample_subset = sample_subset) as pf:
        variant_ct = pf.get_variant_ct()
        hamming_distance = 0
        buf = np.empty(sample_ct, np.int8)
        for vidx in range(variant_ct):
            pf.read(vidx, buf)
            geno0 = buf[0]
            geno1 = buf[1]
            if geno0 != -9 and geno1 != -9:
                hamming_distance += abs(geno0 - geno1)
        sys.stdout.write("Hamming distance: " + str(hamming_distance))
        sys.stdout.write('\n')

if __name__ == "__main__":
    main()
