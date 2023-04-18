#!/usr/bin/env python3
import pgenlib
import numpy as np
import random
import os

def generate_unphased_biallelic_test_genotypes(nsample, nvariant):
    geno_buf = np.zeros([nvariant, nsample], dtype=np.int8)
    for vidx in range(nvariant):
        approx_freq = 0.0
        if random.randrange(2) == 0:
            approx_freq = random.random() * 0.05
        else:
            approx_freq = random.random()
        # This sampling procedure is expected to be inaccurate for larger
        # values of approx_freq, but that doesn't matter.
        one_attempts = int(2 * approx_freq * nvariant)
        for i in range(one_attempts):
            geno_buf[vidx][random.randrange(nsample)] = 1
        two_attempts = int(approx_freq * approx_freq * nvariant)
        for i in range(two_attempts):
            geno_buf[vidx][random.randrange(nsample)] = 2

        # Throw in some missing genotypes.
        while True:
            if random.randrange(2) == 0:
                break
            geno_buf[vidx][random.randrange(nsample)] = -9

    if nvariant >= 3:
        all_alt_vidx = random.randrange(nvariant)
        geno_buf[all_alt_vidx].fill(2)
        all_missing_vidx = random.randrange(nvariant - 1)
        if all_missing_vidx >= all_alt_vidx:
            all_missing_vidx = all_missing_vidx + 1
        geno_buf[all_missing_vidx].fill(-9)

    return geno_buf


def assert_alleles(want_geno, a0, a1):
    if want_geno == -9:
        assert (a0 == -9) and (a1 == -9)
    elif want_geno == 0:
        assert (a0 == 0) and (a1 == 0)
    elif want_geno == 1:
        assert (a0 == 0) and (a1 == 1)
    else:
        assert (want_geno == 2)
        assert (a0 == 1) and (a1 == 1)


def check_unphased_biallelic_read_concordance(r, raw_nsample, nvariant, test_genotypes, sample_subset):
    want_genotypes = test_genotypes
    cur_nsample = raw_nsample
    if sample_subset is not None:
        cur_nsample = len(sample_subset)
        want_genotypes = np.empty([nvariant, cur_nsample], dtype=np.int8)
        for vidx in range(nvariant):
            for i in range(cur_nsample):
                want_genotypes[vidx][i] = test_genotypes[vidx][sample_subset[i]]
    assert raw_nsample == r.get_raw_sample_ct()
    assert nvariant == r.get_variant_ct()

    geno_int = np.empty([cur_nsample], np.int8)
    vidx = random.randrange(nvariant)
    r.read(vidx, geno_int)
    assert np.array_equal(geno_int, want_genotypes[vidx])
    r.read(vidx, geno_int, allele_idx=0)
    for i in range(cur_nsample):
        assert ((geno_int[i] == -9) and (want_genotypes[vidx][i] == -9)) or (geno_int[i] == 2 - want_genotypes[vidx][i])

    allele_int32 = np.empty([2 * cur_nsample], np.int32)
    vidx = random.randrange(nvariant)
    r.read_alleles(vidx, allele_int32)
    for i in range(cur_nsample):
        assert_alleles(want_genotypes[vidx][i], allele_int32[2*i], allele_int32[2*i+1])

    tmp_nvariant = 0
    if nvariant > 258:
        tmp_nvariant = random.randrange(258, nvariant)
    else:
        tmp_nvariant = random.randrange(1, nvariant + 1)
    vidx_start = 0
    if nvariant > tmp_nvariant:
        vidx_start = random.randrange(nvariant - tmp_nvariant)
    vidx_end = vidx_start + tmp_nvariant
    vidx_list = np.asarray(sorted(random.sample(range(nvariant), k=tmp_nvariant)), np.uint32)
    geno_int_2d = np.empty([tmp_nvariant, cur_nsample], np.int8)
    allele_int32_2d = np.empty([tmp_nvariant, 2 * cur_nsample], np.int32)
    geno_int_t = np.empty([cur_nsample, tmp_nvariant], np.int8)

    r.read_range(vidx_start, vidx_end, geno_int_2d)
    np.array_equal(geno_int_2d, want_genotypes[vidx_start:vidx_end])
    r.read_range(vidx_start, vidx_end, geno_int_t, sample_maj=True)
    np.array_equal(geno_int_2d, np.transpose(geno_int_t))

    r.read_list(vidx_list, geno_int_2d)
    for i in range(tmp_nvariant):
        assert np.array_equal(geno_int_2d[i], want_genotypes[vidx_list[i]])
    r.read_list(vidx_list, geno_int_t, sample_maj=True)
    np.array_equal(geno_int_2d, np.transpose(geno_int_t))

    # fewer variants since we're not testing transpose, and this is relatively
    # slow
    if nvariant > 10:
        tmp_nvariant = random.randrange(5, 10)
    else:
        tmp_nvariant = random.randrange(1, nvariant + 1)
    vidx_start = 0
    if nvariant > tmp_nvariant:
        vidx_start = random.randrange(nvariant - tmp_nvariant)
    vidx_end = vidx_start + tmp_nvariant
    vidx_list = np.asarray(sorted(random.sample(range(nvariant), k=tmp_nvariant)), np.uint32)
    allele_int32_2d = np.empty([tmp_nvariant, 2 * cur_nsample], np.int32)
    r.read_alleles_range(vidx_start, vidx_end, allele_int32_2d)
    for i in range(tmp_nvariant):
        for j in range(cur_nsample):
            assert_alleles(want_genotypes[vidx_start+i][j], allele_int32_2d[i][2*j], allele_int32_2d[i][2*j+1])
    r.read_alleles_list(vidx_list, allele_int32_2d)
    for i in range(tmp_nvariant):
        for j in range(cur_nsample):
            assert_alleles(want_genotypes[vidx_list[i]][j], allele_int32_2d[i][2*j], allele_int32_2d[i][2*j+1])

    genocount_int32 = np.empty([4], np.uint32)
    want_counts = np.empty([4], np.uint32)
    for i in range(tmp_nvariant):
        vidx = vidx_list[i]
        r.count(vidx, genocount_int32)
        want_counts.fill(0)
        for j in range(cur_nsample):
            want_geno = want_genotypes[vidx][j]
            if want_geno == -9:
                want_counts[3] += 1
            else:
                want_counts[want_geno] += 1
        assert np.array_equal(genocount_int32, want_counts)
        r.count(vidx, genocount_int32, allele_idx=0)
        assert genocount_int32[0] == want_counts[2]
        assert genocount_int32[1] == want_counts[1]
        assert genocount_int32[2] == want_counts[0]
        assert genocount_int32[3] == want_counts[3]


def unphased_biallelic_case(tmp_path, case_idx, nsample_min, nsample_limit, nvariant_min, nvariant_limit):
    # 1. Generate a matrix of genotypes, where ~half have alt-freq < 0.05, at
    #    least one has alt-freq = 1, and at least one is all-missing.
    # 2. Use a mix of append_biallelic() and append_biallelic_batch() to write
    #    these genotypes to a .pgen.
    # 3. Verify that get_raw_sample_ct(), get_variant_ct(), read(),
    #    read_alleles(), read_range(), read_list(), read_alleles_range(),
    #    read_alleles_list(), count() return expected results.  Make sure
    #    sample_maj=True cases are covered.
    # 4. Call change_sample_subset(), and repeat step 3.
    random.seed(case_idx)
    nsample = random.randrange(nsample_min, nsample_limit)
    nvariant = random.randrange(nvariant_min, nvariant_limit)
    test_genotypes = generate_unphased_biallelic_test_genotypes(nsample, nvariant)
    test_pgen_path = bytes(tmp_path / ("unphased_biallelic_" + str(case_idx) + ".pgen"))
    with pgenlib.PgenWriter(test_pgen_path, nsample, variant_ct=nvariant, nonref_flags=False) as w:
        for vidx in range(nvariant // 3):
            w.append_biallelic(test_genotypes[vidx])
        w.append_biallelic_batch(test_genotypes[nvariant // 3 : (2 * nvariant) // 3])
        for vidx in range((2 * nvariant) // 3, nvariant):
            w.append_biallelic(test_genotypes[vidx])
    with pgenlib.PgenReader(test_pgen_path) as r:
        check_unphased_biallelic_read_concordance(r, nsample, nvariant, test_genotypes, None)
        sample_subset = sorted(random.sample(range(nsample), k=nsample // 2))
        r.change_sample_subset(np.asarray(sample_subset, np.uint32))
        check_unphased_biallelic_read_concordance(r, nsample, nvariant, test_genotypes, sample_subset)


def test_unphased_biallelic(tmp_path):
    ncase = 5
    # We want some cases to have dimensions > 256, so we can catch
    # transpose bugs.
    for case_idx in range(0, ncase):
        unphased_biallelic_case(tmp_path, case_idx, 260, 340, 260, 340)
    # Also cover smaller possibilities.
    for case_idx in range(ncase, 2*ncase):
        unphased_biallelic_case(tmp_path, case_idx, 1, 100, 260, 340)
    for case_idx in range(2*ncase, 3*ncase):
        unphased_biallelic_case(tmp_path, case_idx, 260, 340, 1, 100)
    for case_idx in range(3*ncase, 4*ncase):
        unphased_biallelic_case(tmp_path, case_idx, 1, 100, 1, 100)


def generate_phased_multiallelic_test_acodes(nsample, nvariant, allele_ct_limit):
    allele_codes = np.zeros([nvariant, 2*nsample], dtype=np.int32)
    phasepresent_bytes = np.zeros([nvariant, nsample], dtype=np.uint8)
    allele_idx_offsets = np.zeros([nvariant+1], dtype=np.uintp)
    phasepresent_bytes[(nvariant // 3):nvariant].fill(1)
    for vidx in range(nvariant):
        approx_freq = 0.0
        if random.randrange(2) == 0:
            approx_freq = random.random() * 0.05
        else:
            approx_freq = random.random()
        attempts = int(2 * approx_freq * nvariant)
        for i in range(attempts):
            allele_codes[vidx][random.randrange(2*nsample)] = random.randrange(1, allele_ct_limit)

        # Throw in some missing genotypes.
        while True:
            if random.randrange(2) == 0:
                break
            sample_idx = random.randrange(nsample)
            allele_codes[vidx][2*sample_idx] = -9
            allele_codes[vidx][2*sample_idx+1] = -9
            # in case it's in the all-phased range
            phasepresent_bytes[vidx][sample_idx] = 0

    if nvariant >= 2:
        all_missing_vidx = random.randrange(nvariant)
        allele_codes[all_missing_vidx].fill(-9)
        phasepresent_bytes[all_missing_vidx].fill(0)

    # Now adjust for phase.
    # First third are all-unphased, second third are all-phased, last third is
    # random.
    for vidx in range(nvariant // 3):
        for i in range(nsample):
            ac0 = allele_codes[vidx][2*i]
            ac1 = allele_codes[vidx][2*i+1]
            if (ac0 == ac1):
                if ac0 != -9:
                    # Even in the all-unphased case,
                    # read_alleles_and_phasepresent() will mark homozygous
                    # genotypes as phased.
                    phasepresent_bytes[vidx][i] = 1
            elif ac0 > ac1:
                allele_codes[vidx][2*i], allele_codes[vidx][2*i+1] = ac1, ac0
    for vidx in range((2 * nvariant) // 3, nvariant):
        attempts = nsample // 2
        for i in range(attempts):
            sample_idx = random.randrange(nsample)
            ac0 = allele_codes[vidx][2*sample_idx]
            ac1 = allele_codes[vidx][2*sample_idx+1]
            if ac0 != ac1:
                phasepresent_bytes[vidx][sample_idx] = 0
                if ac0 > ac1:
                    allele_codes[vidx][2*sample_idx] = ac1
                    allele_codes[vidx][2*sample_idx+1] = ac0
    offset = 0
    for vidx in range(nvariant):
        max_allele_code_p1 = 1 + np.amax(allele_codes[vidx], initial=1)
        if max_allele_code_p1 < allele_ct_limit:
            if random.randrange(10) == 0:
                max_allele_code_p1 += 1
        offset += max_allele_code_p1
        allele_idx_offsets[vidx+1] = offset

    return allele_codes, phasepresent_bytes, allele_idx_offsets


def check_phased_multiallelic_read_concordance(r, raw_nsample, nvariant, test_acodes, test_phasepresent_bytes, sample_subset):
    want_acodes = test_acodes
    want_phasepresent_bytes = test_phasepresent_bytes
    cur_nsample = raw_nsample
    if sample_subset is not None:
        cur_nsample = len(sample_subset)
        want_acodes = np.empty([nvariant, 2*cur_nsample], dtype=np.int32)
        want_phasepresent_bytes = np.empty([nvariant, cur_nsample], dtype=np.uint8)
        for vidx in range(nvariant):
            for i in range(cur_nsample):
                src_i = sample_subset[i]
                want_acodes[vidx][2*i] = test_acodes[vidx][2*src_i]
                want_acodes[vidx][2*i+1] = test_acodes[vidx][2*src_i+1]
                want_phasepresent_bytes[vidx][i] = test_phasepresent_bytes[vidx][src_i]
    assert raw_nsample == r.get_raw_sample_ct()
    assert nvariant == r.get_variant_ct()

    allele_int32 = np.empty([2 * cur_nsample], np.int32)
    phasepresent_bytes = np.empty([cur_nsample], np.uint8)
    vidx = random.randrange(nvariant)
    r.read_alleles_and_phasepresent(vidx, allele_int32, phasepresent_bytes)
    assert np.array_equal(allele_int32, want_acodes[vidx])
    assert np.array_equal(phasepresent_bytes, want_phasepresent_bytes[vidx])


def phased_multiallelic_case(tmp_path, case_idx, nsample_min, nsample_limit, nvariant_min, nvariant_limit, allele_ct_max):
    # 1. Generate allele_codes and phasepresent_bytes, with some all-unphased
    #    and some all-phased variants.
    # 2. Use a mix of append_alleles(), append_partially_phased(),
    #    append_alleles_batch(), and append_partially_phased_batch() to write
    #    these genotypes to a .pgen.
    # 3. Verify that get_raw_sample_ct(), get_variant_ct(),
    #    read_alleles_and_phasepresent(),
    #    read_alleles_and_phasepresent_range(),
    #    read_alleles_and_phasepresent_list(), and count() return expected
    #    results.
    # 4. Call change_sample_subset(), and repeat step 3.
    random.seed(case_idx)
    nsample = random.randrange(nsample_min, nsample_limit)
    nvariant = random.randrange(nvariant_min, nvariant_limit)
    allele_ct_limit = random.randrange(2, allele_ct_max + 1)
    test_acodes, test_phasepresent_bytes, test_allele_idx_offsets = generate_phased_multiallelic_test_acodes(nsample, nvariant, allele_ct_limit)
    test_pgen_path = bytes(tmp_path / ("phased_multiallelic_" + str(case_idx) + ".pgen"))
    with pgenlib.PgenWriter(test_pgen_path, nsample, variant_ct=nvariant, nonref_flags=False, allele_ct_limit=allele_ct_limit, hardcall_phase_present=True) as w:
        for vidx in range(nvariant // 6):
            w.append_alleles(test_acodes[vidx])
        w.append_alleles_batch(test_acodes[nvariant // 6 : nvariant // 3])
        for vidx in range(nvariant // 3, nvariant // 2):
            w.append_alleles(test_acodes[vidx], all_phased=True)
        w.append_alleles_batch(test_acodes[nvariant // 2 : (2 * nvariant) // 3], all_phased=True)
        for vidx in range((2 * nvariant) // 3, (5 * nvariant) // 6):
            w.append_partially_phased(test_acodes[vidx], test_phasepresent_bytes[vidx])
        w.append_partially_phased_batch(test_acodes[(5 * nvariant) // 6 : nvariant], test_phasepresent_bytes[(5 * nvariant) // 6 : nvariant])
    with pgenlib.PgenReader(test_pgen_path, allele_idx_offsets=test_allele_idx_offsets) as r:
        check_phased_multiallelic_read_concordance(r, nsample, nvariant, test_acodes, test_phasepresent_bytes, None)
        sample_subset = sorted(random.sample(range(nsample), k=nsample // 2))
        r.change_sample_subset(np.asarray(sample_subset, np.uint32))
        check_phased_multiallelic_read_concordance(r, nsample, nvariant, test_acodes, test_phasepresent_bytes, sample_subset)


def test_phased_multiallelic(tmp_path):
    ncase = 5
    for case_idx in range(0, ncase):
        phased_multiallelic_case(tmp_path, case_idx, 1, 100, 1, 100, 2)
    for case_idx in range(ncase, 2*ncase):
        phased_multiallelic_case(tmp_path, case_idx, 1, 100, 1, 100, 255)
