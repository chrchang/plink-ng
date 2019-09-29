#!/usr/bin/env python3
import pgenlib
import numpy as np
import sys

def main():
    arg_ct = len(sys.argv)
    if arg_ct < 2:
        print("Usage: python3 het_count.py <.bed/.pgen> [raw_sample_ct=<val>] [sample idx(s)...]")
        print("* raw_sample_ct is required for .bed files.")
        print("* sample indexes must be in increasing order.")
        print("* The Python side of this is really slow if you don't explicitly specify sample indexes.")
        return
    sample_subset = None
    specified_sample_ct = None
    sample_ct = None
    if arg_ct > 2:
        offset = 2
        if sys.argv[2].startswith("raw_sample_ct="):
            specified_sample_ct = int(sys.argv[2][14:])
            offset = 3
        if arg_ct > offset:
            sample_ct = arg_ct - offset
            sample_subset = np.empty(sample_ct, np.uint32)
            for idx in range(sample_ct):
                sample_subset[idx] = int(sys.argv[offset + idx])
    with pgenlib.PgenReader(bytes(sys.argv[1], 'utf8'), raw_sample_ct = specified_sample_ct, sample_subset = sample_subset) as pf:
        raw_sample_ct = pf.get_raw_sample_ct()
        if sample_ct is None:
            sample_ct = raw_sample_ct
        variant_ct = pf.get_variant_ct()
        tot_hets = np.zeros((sample_ct,), dtype=np.uint32)
        buf = np.empty(sample_ct, np.int8)
        for vidx in range(variant_ct):
            pf.read(vidx, buf)
            # this is horribly inefficient...
            for sample_idx in range(sample_ct):
                if buf[sample_idx] == 1:
                    tot_hets[sample_idx] += 1
        for sample_idx in range(sample_ct):
            sys.stdout.write(str(tot_hets[sample_idx]))
            sys.stdout.write(' ')
        sys.stdout.write('\n')

if __name__ == "__main__":
    main()
