#!/usr/bin/env python3
import pgenlib
import numpy as np
import sys

def main():
    arg_ct = len(sys.argv)
    if arg_ct < 3:
        print("Usage: python single_variant_test.py <.bed/.pgen> <variant idx> [raw_sample_ct=<val>] [sample idx(s)...]")
        print("* raw_sample_ct is required for .bed files.")
        print("* sample indexes must be in increasing order.")
        return
    sample_subset = None
    specified_sample_ct = None
    sample_ct = None
    if arg_ct > 3:
        offset = 3
        if sys.argv[3].startswith("raw_sample_ct="):
            specified_sample_ct = int(sys.argv[3][14:])
            offset = 4
        if arg_ct > offset:
            sample_ct = arg_ct - offset
            sample_subset = np.empty(sample_ct, np.uint32)
            for idx in range(sample_ct):
                sample_subset[idx] = int(sys.argv[offset + idx])
    vidx = int(sys.argv[2])
    with pgenlib.PgenReader(bytes(sys.argv[1], 'utf8'), raw_sample_ct = specified_sample_ct, sample_subset = sample_subset) as pf:
        raw_sample_ct = pf.get_raw_sample_ct()
        if sample_ct is None:
            sample_ct = raw_sample_ct
        # variant_ct = pf.get_variant_ct()
        buf = np.empty(raw_sample_ct * 2, np.int32)
        buf2 = np.empty(raw_sample_ct, np.bool_)
        pf.read_alleles_and_phasepresent(vidx, buf, buf2)
        for sample_idx in range(sample_ct):
            sys.stdout.write(str(buf[2 * sample_idx]))
            if buf2[sample_idx]:
                sys.stdout.write('|')
            else:
                sys.stdout.write('/')
            sys.stdout.write(str(buf[2 * sample_idx + 1]))
            sys.stdout.write(' ')
        sys.stdout.write('\n')
        pf.change_sample_subset()
        pf.read_alleles_and_phasepresent(vidx, buf, buf2)
        for sample_idx in range(raw_sample_ct):
            sys.stdout.write(str(buf[2 * sample_idx]))
            if buf2[sample_idx]:
                sys.stdout.write('|')
            else:
                sys.stdout.write('/')
            sys.stdout.write(str(buf[2 * sample_idx + 1]))
            sys.stdout.write(' ')
        sys.stdout.write('\n')

if __name__ == "__main__":
    main()
