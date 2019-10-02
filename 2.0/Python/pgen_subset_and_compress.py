#!/usr/bin/env python3
import pgenlib
import numpy as np
import sys

def main():
    arg_ct = len(sys.argv)
    if arg_ct < 3:
        print("Usage: python3 pgen_subset_and_compress.py <input .bed/.pgen> <output .pgen name> [raw_sample_ct=<val>] [sample idx(s)...]")
        print("* raw_sample_ct is required for .bed files.")
        print("* sample indexes must be in increasing order.")
        print("* This currently assumes that A2 alleles are always reference.")
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
    with pgenlib.PgenReader(bytes(sys.argv[1], 'utf8'), raw_sample_ct = specified_sample_ct, sample_subset = sample_subset) as infile:
        raw_sample_ct = infile.get_raw_sample_ct()
        if sample_ct is None:
            sample_ct = raw_sample_ct
        variant_ct = infile.get_variant_ct()
        hardcall_phase_present = infile.hardcall_phase_present()
        geno_buf = None
        allele_code_buf = None
        phasepresent_buf = None
        if not hardcall_phase_present:
            geno_buf = np.empty(sample_ct, np.int8)
        else:
            allele_code_buf = np.empty(sample_ct * 2, np.int32)
            phasepresent_buf = np.empty(sample_ct, np.bool_)
        with pgenlib.PgenWriter(bytes(sys.argv[2], 'utf8'), sample_ct, variant_ct, False, hardcall_phase_present=hardcall_phase_present) as outfile:
            if not hardcall_phase_present:
                for vidx in range(variant_ct):
                    infile.read(vidx, geno_buf)
                    outfile.append_biallelic(geno_buf)
            else:
                for vidx in range(variant_ct):
                    infile.read_alleles_and_phasepresent(vidx, allele_code_buf, phasepresent_buf)
                    outfile.append_partially_phased(allele_code_buf, phasepresent_buf)
    return


if __name__ == "__main__":
    main()
