Pgenlib Python API specification, v0.7
("import numpy as np" assumed)

class PgenReader:
* PgenReader(filename, raw_sample_ct = None, variant_ct = None,
             sample_subset = None)
  Constructor, opens the .pgen or .bed file.  Requires a filename (exception is
  thrown if file doesn't exist or it's invalid).
  - raw_sample_ct is required for a .bed file (otherwise an exception will be
    thrown), and optional for a .pgen file.  If it's provided for a .pgen file,
    an exception is thrown if sample_ct does not match the value in the .pgen.
  - variant_ct is always optional.  An exception is thrown if it does not match
    the value implied by the .bed file's size, or the explicitly stored value
    in the .pgen.
  - sample_subset is an optional numpy uint32 array of (0-based) indexes in
    increasing order, telling the reader to only consider the specified samples
    when loading or counting genotypes.  For example, if your .fam file looks
    something like
      40184_187545456 40184_187545456 0       0       1       -9
      40195_187545457 40195_187545457 0       0       1       -9
      40206_187545458 40206_187545458 0       0       1       -9
      40217_187545459 40217_187545459 0       0       2       -9
      40228_187545460 40228_187545460 0       0       1       -9
      40239_187545461 40239_187545461 0       0       1       -9
      40250_187520807 40250_187520807 0       0       1       -9
      40261_187520806 40261_187520806 0       0       1       -9
      40272_187520805 40272_187520805 0       0       2       -9
      40283_187520804 40283_187520804 0       0       1       -9
    then a [1, 4, 5, 6] array specifies samples 40195_187545457,
    40228_187545460, 40239_187545461, and 40250_187520807.
    None indicates that all samples should be included.  Otherwise, an
    exception is thrown if any value is not an integer in 0..(raw_sample_ct -
    1), or the values are not strictly increasing, or the array is empty.
    (Empty array might be permitted in the future.)

* get_raw_sample_ct()
  Returns the number of samples in the .pgen file.

* get_variant_ct()
  Returns the number of variants in the .pgen file.

* hardcall_phase_present()
  Returns True iff phased hardcalls may be present.

* read(uint32_t variant_idx, np.ndarray[np.int{8,32,64}_t] geno_int_out,
       uint32_t allele_idx = 1)
  Takes a (0-based) variant index; fills geno_int_out, which must be a numpy
  int8, int32, or int64 array with at least sample_ct entries, where sample_ct
  (as distinguished from raw_sample_ct) is the size of the current sample
  subset.  These values are in {0, 1, 2, -9}; by default, 0/1/2 represents the
  alternate allele count.  Setting allele_idx = 0 causes reference allele
  counts to be reported instead.  (If we ever use multiallelic variant records,
  setting allele_idx = 2 causes alternate allele 2 counts, etc.)
  - I tried the interface without the required outbuffer, and it could take
    twice as long due to all the extra allocation/deallocation.  The buffer can
    be allocated with
      buf = np.empty(sample_ct, np.int{8,32,64})
  - An exception is thrown if variant_idx or allele_idx is invalid (or there's
    e.g. a file I/O error, or the file isn't open)
  - Note that pgenlib is unaware of chromosome boundaries.  The caller is
    responsible for chrX/chrY/chrM-specific corrections for now.

* read_dosages(uint32_t variant_idx,
               np.ndarray[np.float{32,64}_t] floatarr_out,
	       uint32_t allele_idx = 1)
  Takes a (0-based) variant index; fills floatarr_out, which must be a numpy
  float32 or float64 array with at least sample_ct entries.  Missing entries
  are encoded as -9, everything else is a dosage in [0, 2] (more precisely,
  they'll be a multiple of 2^{-14}).

* read_alleles(uint32_t variant_idx, np.ndarray[np.int32_t] allele_int32_out)
* read_alleles_and_phasepresent(uint32_t variant_idx,
                                np.ndarray[np.int32_t] allele_int32_out,
                                np.ndarray[np.uint8_t,cast=True] phasepresent_out)
  Takes a (0-based) variant index.
  - allele_int32_out must have space for at least (2 * sample_ct) allele
    indexes, where elements 2n and (2n+1) correspond to sample n.  Both indexes
    are -9 if the genotype is missing.  If the genotype is unphased, the lower
    index appears first.
  - For read_alleles_and_phasepresent(), phasepresent_out should be a numpy
    bool_, int8, or uint8 array.  Element n is set to true iff the genotype for
    sample n has known phase.  Most of these values will be true even when the
    raw data is unphased, because homozygous genotypes always have known phase.
    (Missing genotypes are considered to have unknown phase, of course.)
  - If int32_t is problematically large, you probably want to work with a 2-bit
    representation rather than int8_t.  That's outside the scope of this API;
    you may instead want to adapt the appropriate pgenlib.pyx function.

* read_range(uint32_t variant_idx_start, uint32_t variant_idx_end,
             np.ndarray[np.int{8,32,64}_t,mode="c",ndim=2] geno_int_out,
	     uint32_t allele_idx = 1, bint sample_maj = 0)
  read_list(np.ndarray[np.uint32_t] variant_idxs,
            np.ndarray[np.int{8,32,64}_t,mode="c",ndim=2] geno_int_out,
	    uint32_t allele_idx = 1, bint sample_maj = 0)
  read_alleles_range(uint32_t variant_idx_start, uint32_t variant_idx_end,
                     np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out,
		     bint hap_maj = 0)
  read_alleles_list(np.ndarray[np.uint32_t] variant_idxs,
                    np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out,
		    bint hap_maj = 0)
  read_alleles_and_phasepresent_range(uint32_t variant_idx_start,
                                      uint32_t variant_idx_end,
                                      np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out,
			              np.ndarray[np.uint8_t,mode="c",cast=True,ndim=2] phasepresent_out,
			              bint hap_maj = 0)
  read_alleles_and_phasepresent_list(np.ndarray[np.uint32_t] variant_idxs,
                                     np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out,
			             np.ndarray[np.uint8_t,mode="c",cast=True,ndim=2] phasepresent_out,
			             bint hap_maj = 0)
  Read functions which handle multiple variants at once, saving data to 2D
  arrays.  By default, the return arrays are variant-major (# of columns is
  based on sample_ct, # of rows is equal to the block size); set
  sample_maj(/hap_maj) = True to make them sample(/haplotype)-major instead.
  For the _range() functions, the interval is half-open (i.e. variant_idx_start
  is included, variant_idx_end - 1 is included, variant_idx_end is not).  For
  the _list() functions, it's currently okay for the variant indexes to be
  unsorted, or for duplicates to be present, but that may not remain true.
  (read_alleles_and_phasepresent_{range,list} not implemented yet.)
  (Todo: read_dosages_range() and read_dosages_list().)

* count(uint32_t variant_idx, np.ndarray[np.uint32_t] genocount_uint32_out,
        uint32_t allele_idx = 1)
  Takes a (0-based) variant index, fills genocount_uint32_out (which must have
  at least 4 slots).  First element corresponds to the number of genotypes with
  0 alt1 alleles, the second element corresponds to the # with one alt1 allele,
  the third corresponds to homozygous-alt1, and the fourth corresponds to the
  number of missing genotypes.
  - Setting allele_idx = 0 replaces "alt1" with "ref" above.
  - Setting allele_idx = None is functionally equivalent to allele_idx = 1 on
    all-biallelic data, but the length of the array is no longer limited to 4
    for multiallelic variants; instead it has 1 + v(v+1)/2 elements where v is
    the number of alleles, in 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, 0/3, ...,
    (v-1)/(v-1), ./. order.  (This matches VCF FORMAT:GL order.)
  (Todo: check whether multi-variant analogues of read/read_phased/count have
  any value.)

* phased_samples() [not written yet]
  If phase set data is present, or the entire dataset is unphased, returns a
  numpy bool_ array where element n is true iff sample n has at least one phase
  set.  Throws an exception if phased variants are present, but no phase set
  data track is present.  (Note that the phase set data track is not
  implemented yet, but it should be within the next few months.)

* change_sample_subset(sample_subset = None)
  Changes the subset of samples to read.
  - sample_subset format is the same as for the constructor.  (This isn't
    really optimized; it's assumed that you won't be calling this often.  If
    you ever need a higher-performance version of this, let me know.)

* close()
  Closes the file, no further read operations can be performed.
  - Not strictly necessary when e.g. a with-statement is used, but still a good
    idea in longer scripts since otherwise object cleanup/file close may be
    delayed.


class PgenWriter:
* PgenWriter(filename, sample_ct, variant_ct, nonref_flags,
             allele_idx_offsets = None, hardcall_phase_present = False,
	     dosage_present = False, dosage_phase_present = False)
  Constructor, creates a new .pgen file and writes a mostly-empty header (which
  gets filled at the end).
  - sample_ct and variant_ct must be positive (and less than about 2^31).
  - nonref_flags should be True when the data is from an ordinary PLINK 1.x
    .bed file (where the A2 allele is major rather than consistently
    reference), False when the A2 allele *is* consistently reference, or a
    numpy bool_ array of length variant_ct when this is mixed.  (Strictly
    speaking, None is also permitted--that delegates tracking of nonref
    information to the .pvar file--but it's discouraged, since .pgen+.bim+.fam
    is a useful data representation with direct plink2 support.)
  - allele_idx_offsets is an optional numpy intp array of length
    (variant_ct+1), where allele_idx_offsets[0] = 0, and the number of alleles
    for (0-based) variant n is (allele_idx_offsets[n+1]-allele_idx_offsets[n]).
    # of alleles must be at least 2 for each variant.  If allele_idx_offsets
    is None, all variants are assumed to be biallelic.

* append_biallelic(genobytes)
  Takes a numpy int8 array with sample_ct {0, 1, 2, -9} elements, and appends
  it to the .pgen.  Ok to use 3 in place of -9.

* append_alleles(allele_codes, all_phased = False)
  append_partially_phased(allele_codes, np.ndarray[np.uint8_t,cast=True] phasepresent)
  Takes a numpy int32 array with (2 * sample_ct) allele codes (0 = ref,
  1 = alt1, -9 = missing), and appends it to the .pgen.  -9s must occur in
  pairs.
  - With append_alleles(), all genotypes are treated as unphased by default.
    Set all_phased to True to treat all genotypes as phased instead.
  - With append_partially_phased(), phasepresent should be a numpy bool_ (or
    equivalent int8_t/uint8_t; in that case, all values must be 0s and 1s)
    array of length sample_ct.  Entries for non-heterozygous calls are ignored.
  - It's fine for an unphased het call to be stored in 1/0 instead of 0/1
    order.

* append_dosages(floatarr)
  Takes a numpy float32 or float64 array of dosages, and appends it to the
  .pgen.
  - The writer must have been initialized with dosage_present = True.
  - Regular dosages are expected to be in [0, 2] (you'll want to multiply
    haploid dosages by 2).
  - Except for small (2^{-16}) tolerances around 0 and 2, any out-of-range
    value (including -9) is interpreted as missing.

* append_biallelic_batch(np.ndarray[np.int8_t,mode="c",ndim=2] genobytes_batch)
  append_alleles_batch(np.ndarray[np.int32_t,mode="c",ndim=2] allele_codes_batch,
                       all_phased = False)
  append_partially_phased_batch(np.ndarray[np.int32_t,mode="c",ndim=2] allele_codes_batch,
                                np.ndarray[np.{int8,uint8}_t,mode="c",ndim=2] phasepresent_batch)
  append_dosages_batch(np.ndarray[np.float{32,64}_t] floatarr_batch)
  Multi-variant forms of append_biallelic(), append_alleles(),
  append_partially_phased(), and append_dosages().  Input matrices must be
  variant-major.

* close()
  Backfills the header and closes the file; no further append operations can be
  performed.  Throws an exception if the number of written variants is not
  equal to the initially-declared variant_ct.
