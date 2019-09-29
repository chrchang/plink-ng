#!/usr/bin/env python3
import pgenlib
import numpy as np
import sys

def main():
    arg_ct = len(sys.argv)
    if arg_ct < 3:
        print("Usage: python3 extract_haps.py <.bed/.pgen> <output filename> [raw_sample_ct=<val>] [sample idx(s)...]")
        print("* raw_sample_ct is required for .bed files.")
        return
    sample_subset = None
    specified_sample_ct = None
    sample_ct = None
    if arg_ct > 3:
        offset = 3
        if sys.argv[3].startswith("raw_sample_ct="):
            specified_sample_ct = int(sys.argv[3][14:])
        if arg_ct > offset:
            sample_ct = arg_ct - offset
            sample_subset = np.empty(sample_ct, np.uint32)
            for idx in range(sample_ct):
                sample_subset[idx] = int(sys.argv[offset + idx])
    with pgenlib.PgenReader(bytes(sys.argv[1], 'utf8'), raw_sample_ct = specified_sample_ct, sample_subset = sample_subset) as infile:
        raw_sample_ct = infile.get_raw_sample_ct()
        if sample_ct is None:
            sample_ct = raw_sample_ct
        variant_ct = infile.get_variant_ct()
        allele_code_buf = np.empty([sample_ct * 2, variant_ct], dtype=np.int32)
        infile.read_alleles_range(0, variant_ct, allele_code_buf, 1)
        with open(sys.argv[2], 'w') as outfile:
            for vidx in range(variant_ct):
                for sidx in range(sample_ct):
                    if sidx != 0:
                        outfile.write(' ')
                    outfile.write(str(allele_code_buf[2 * sidx, vidx]) + "|" + str(allele_code_buf[2 * sidx + 1, vidx]))
                outfile.write('\n')

if __name__ == "__main__":
    main()
