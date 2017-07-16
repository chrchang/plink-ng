# from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t, uintptr_t, uint32_t, int32_t, uint16_t, uint8_t, int8_t
from cpython.mem cimport PyMem_Malloc, PyMem_Free
# from cpython.view cimport array as cvarray
import numpy as np
cimport numpy as np
import sys

cdef extern from "../pgenlib_python_support.h":
    # macros aren't namespaced
    uintptr_t DIV_UP(uintptr_t val, uintptr_t divisor)

cdef extern from "../pgenlib_python_support.h" namespace "plink2":
    ctypedef uint32_t boolerr_t
    ctypedef enum pglerr_t:
        kPglRetSuccess
        kPglRetSkipped
        kPglRetNomem
        kPglRetOpenFail
        kPglRetReadFail
        kPglRetWriteFail
        kPglRetMalformedInput
        kPglRetInconsistentInput
        kPglRetInvalidCmdline
        kPglRetHelp
        kPglRetThreadCreateFail
        kPglRetNetworkFail
        kPglRetSampleMajorBed
        kPglRetImproperFunctionCall
        kPglRetNotYetSupported
        kPglRetLongLine
        kPglRetEmptyFile

    boolerr_t cachealigned_malloc(uintptr_t size, void* aligned_pp)
    void aligned_free(void* aligned_ptr)

    void fill_ulong_zero(uintptr_t entry_ct, uintptr_t* ularr)
    # void bitvec_and(const uintptr_t* arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec)
    void fill_interleaved_mask_vec(const uintptr_t* subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec)
    void fill_cumulative_popcounts(const uintptr_t* subset_mask, uint32_t word_ct, uint32_t* cumulative_popcounts)

    ctypedef uintptr_t vul_t
    void transpose_quaterblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, vul_t* vecaligned_buf)
    void transpose_bitblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, vul_t* vecaligned_buf)

    void genovec_invert_unsafe(uint32_t sample_ct, uintptr_t* genovec)
    void biallelic_dosage16_invert(uint32_t dosage_ct, uint16_t* dosage_vals)

    void genoarr_to_bytes_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int8_t* genobytes)
    void genoarr_to_int32s_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* geno_int32)
    void genoarr_to_int64s_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int64_t* geno_int64)
    void genoarr_to_allele_codes(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* allele_codes)
    void genoarr_phased_to_allele_codes(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t sample_ct, uint32_t phasepresent_ct, unsigned char* phasebytes, int32_t* allele_codes)
    void genoarr_phased_to_hap_codes(const uintptr_t* genoarr, const uintptr_t* phaseinfo, uint32_t variant_batch_size, int32_t* hap0_codes_iter, int32_t* hap1_codes_iter)
    void dosage16_to_floats_minus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_vals, uint32_t sample_ct, uint32_t dosage_ct, float* geno_float)
    void dosage16_to_doubles_minus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_vals, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double)
    void bytes_to_bits_unsafe(const uint8_t* boolbytes, uint32_t sample_ct, uintptr_t* bitarr)
    void bytes_to_genoarr_unsafe(const int8_t* genobytes, uint32_t sample_ct, uintptr_t* genoarr)
    void allele_codes_to_genoarr_unsafe(const int32_t* allele_codes, const unsigned char* phasepresent_bytes, uint32_t sample_ct, uintptr_t* genoarr, uintptr_t* phasepresent, uintptr_t* phaseinfo)
    void floats_to_dosage16(const float* floatarr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr)
    void doubles_to_dosage16(const double* doublearr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr)

    cdef enum:
        k1LU
    cdef enum:
        kCacheline
    cdef enum:
        kBitsPerWord
    cdef enum:
        kBitsPerWordD2
    cdef enum:
        kBitsPerVec
    cdef enum:
        kBitsPerCacheline
    cdef enum:
        kQuatersPerVec
    cdef enum:
        kQuatersPerCacheline
    cdef enum:
        kBytesPerVec
    cdef enum:
        kInt32PerVec
    cdef enum:
        kInt32PerCacheline
    cdef enum:
        kWordsPerVec
    cdef enum:
        kPglErrstrBufBlen

    cdef enum:
        kPglQuaterTransposeBatch
    cdef enum:
        kPglQuaterTransposeWords
    cdef enum:
        kPglQuaterTransposeBufbytes
    cdef enum:
        kPglBitTransposeBufbytes

    ctypedef uint32_t pgen_global_flags_t
    cdef enum:
        kfPgenGlobal0
    cdef enum:
        kfPgenGlobalHardcallPhasePresent
    cdef enum:
        kfPgenGlobalDosagePresent

    cdef cppclass pgen_file_info_t:
        uint32_t raw_variant_ct
        uint32_t raw_sample_ct        
        unsigned char* vrtypes
        uint32_t gflags

    ctypedef uint32_t pgen_header_ctrl_t

    void pgfi_preinit(pgen_file_info_t* pgfip)
    
    pglerr_t pgfi_init_phase1(const char* fname, uint32_t raw_variant_ct, uint32_t raw_sample_ct, uint32_t use_mmap, pgen_header_ctrl_t* header_ctrl_ptr, pgen_file_info_t* pgfip, uintptr_t* pgfi_alloc_cacheline_ct_ptr, char* errstr_buf)

    pglerr_t pgfi_init_phase2(pgen_header_ctrl_t header_ctrl, uint32_t allele_cts_already_loaded, uint32_t nonref_flags_already_loaded, uint32_t use_blockload, uint32_t vblock_idx_start, uint32_t vidx_end, uint32_t* max_vrec_width_ptr, pgen_file_info_t* pgfip, unsigned char* pgfi_alloc, uintptr_t* pgr_alloc_cacheline_ct_ptr, char* errstr_buf)
    
    cdef cppclass pgen_reader_t:
        pgen_file_info_t fi
        unsigned char* fread_buf

    void pgr_preinit(pgen_reader_t* pgrp)

    pglerr_t pgr_init(const char* fname, uint32_t max_vrec_width, pgen_file_info_t* pgfip, pgen_reader_t* pgrp, unsigned char* pgr_alloc)

    pglerr_t pgr_read_allele_countvec_subset_unsafe(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, pgen_reader_t* pgrp, uintptr_t* allele_countvec)

    pglerr_t pgr_read_refalt1_genovec_hphase_subset_unsafe(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo, uint32_t* phasepresent_ct_ptr)

    pglerr_t pgr_read_refalt1_genovec_dosage16_subset_unsafe(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uintptr_t* genovec, uintptr_t* dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr, uint32_t* is_explicit_alt1_ptr)
    
    pglerr_t pgr_get_refalt1_genotype_counts(const uintptr_t* sample_include, const uintptr_t* sample_include_interleaved_vec, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, pgen_reader_t* pgrp, uint32_t* genocounts)
    
    boolerr_t pgfi_cleanup(pgen_file_info_t* pgfip)
    boolerr_t pgr_cleanup(pgen_reader_t* pgrp)

    cdef cppclass pgen_writer_common_t:
        uint32_t variant_ct
        uint32_t sample_ct
        uintptr_t* allele_idx_offsets
        uint32_t vidx

    cdef cppclass st_pgen_writer_t:
        pgen_writer_common_t pwc

    pglerr_t spgw_init_phase1(const char* fname, uintptr_t* allele_idx_offsets, uintptr_t* explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, pgen_global_flags_t phase_dosage_gflags, uint32_t nonref_flags_storage, st_pgen_writer_t* spgwp, uintptr_t* alloc_cacheline_ct_ptr, uint32_t* max_vrec_len_ptr)

    void spgw_init_phase2(uint32_t max_vrec_len, st_pgen_writer_t* spgwp, unsigned char* spgw_alloc)

    pglerr_t spgw_append_biallelic_genovec(const uintptr_t* genovec, st_pgen_writer_t* spgwp)
    
    pglerr_t spgw_append_biallelic_genovec_hphase(const uintptr_t* genovec, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, st_pgen_writer_t* spgwp)

    pglerr_t spgw_append_biallelic_genovec_dosage16(const uintptr_t* genovec, const uintptr_t* dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, st_pgen_writer_t* spgwp)
    
    pglerr_t spgw_finish(st_pgen_writer_t* spgwp)
    
    boolerr_t spgw_cleanup(st_pgen_writer_t* spgwp)

    
cdef class PgenReader:
    # todo: nonref_flags, multiallelic variant support
    cdef pgen_file_info_t* _info_ptr
    cdef pgen_reader_t* _state_ptr
    cdef uintptr_t* _subset_include_vec
    cdef uintptr_t* _subset_include_interleaved_vec
    cdef uint32_t* _subset_cumulative_popcounts
    cdef uint32_t _subset_size
    # preallocate buffers we'll use repeatedly
    cdef uintptr_t* _genovec
    cdef uintptr_t* _phasepresent
    cdef uintptr_t* _phaseinfo
    cdef uintptr_t* _dosage_present
    cdef uint16_t* _dosage_vals
    cdef vul_t* _transpose_batch_buf
    # for multi-variant load-and-transpose, we load up to
    # kPglQuaterTransposeBatch (= 256) variants at a time, and then transpose
    cdef uintptr_t* _multivar_vmaj_geno_buf
    cdef uintptr_t* _multivar_vmaj_phasepresent_buf
    cdef uintptr_t* _multivar_vmaj_phaseinfo_buf
    cdef uintptr_t* _multivar_smaj_geno_batch_buf
    cdef uintptr_t* _multivar_smaj_phaseinfo_batch_buf
    cdef uintptr_t* _multivar_smaj_phasepresent_batch_buf

    cdef set_sample_subset_internal(self, np.ndarray[np.uint32_t,mode="c",ndim=1] sample_subset):
        cdef uint32_t raw_sample_ct = self._info_ptr[0].raw_sample_ct
        cdef uint32_t raw_sample_ctv = DIV_UP(raw_sample_ct, kBitsPerVec)
        cdef uint32_t raw_sample_ctaw = raw_sample_ctv * kWordsPerVec
        cdef uintptr_t* sample_include = self._subset_include_vec
        fill_ulong_zero(raw_sample_ctaw, sample_include)
        cdef uint32_t subset_size = sample_subset.size
        if subset_size == 0:
            raise RuntimeError("Empty sample_subset is not currently permitted.")
        cdef sample_uidx = sample_subset[0]
        cdef uint32_t idx = 0
        cdef next_uidx
        while True:
            if sample_uidx >= raw_sample_ct:
                raise RuntimeError("0-based sample idx too large (" + str(sample_uidx) + "; only " + str(raw_sample_ct) + " in file).")
            sample_include[sample_uidx / kBitsPerWord] |= k1LU << (sample_uidx % kBitsPerWord)
            idx += 1
            if idx == subset_size:
                break
            next_uidx = sample_subset[idx]

            # prohibit this since it implies that the caller expects genotypes
            # to be returned in a different order
            if next_uidx <= sample_uidx:
                raise RuntimeError("sample_subset is not in strictly increasing order.")
            
            sample_uidx = next_uidx

        fill_interleaved_mask_vec(sample_include, raw_sample_ctv, self._subset_include_interleaved_vec)
        
        cdef uint32_t raw_sample_ctl = DIV_UP(raw_sample_ct, kBitsPerWord)
        fill_cumulative_popcounts(sample_include, raw_sample_ctl, self._subset_cumulative_popcounts)

        self._subset_size = subset_size
        return

    
    def __cinit__(self, bytes filename, object raw_sample_ct = None,
                  object variant_ct = None, object sample_subset = None):
        self._info_ptr = <pgen_file_info_t*>PyMem_Malloc(sizeof(pgen_file_info_t))
        if not self._info_ptr:
            raise MemoryError()
        pgfi_preinit(self._info_ptr)
        # this depends on pgenlib_internal implementation.  could save
        # pgfi_alloc and pgr_alloc instead.
        self._info_ptr[0].vrtypes = NULL
        cdef uint32_t cur_sample_ct = 0xffffffffU
        if raw_sample_ct is not None:
            cur_sample_ct = raw_sample_ct
        cdef uint32_t cur_variant_ct = 0xffffffffU
        if variant_ct is not None:
            cur_variant_ct = variant_ct
        cdef const char* fname = <const char*>filename
        cdef pgen_header_ctrl_t header_ctrl
        cdef uintptr_t pgfi_alloc_cacheline_ct
        cdef char errstr_buf[kPglErrstrBufBlen]
        if pgfi_init_phase1(fname, cur_variant_ct, cur_sample_ct, 0, &header_ctrl, self._info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != kPglRetSuccess:
            raise RuntimeError(errstr_buf[7:])
        assert (header_ctrl & 0x30) == 0 # no alt allele counts
        assert (header_ctrl & 0xc0) != 0xc0 # no explicit nonref_flags
        cdef uint32_t file_sample_ct = self._info_ptr[0].raw_sample_ct
        assert file_sample_ct != 0
        cdef unsigned char* pgfi_alloc = NULL
        if pgfi_alloc_cacheline_ct != 0:
            if cachealigned_malloc(pgfi_alloc_cacheline_ct * kCacheline, &pgfi_alloc):
                raise MemoryError()
        cdef uint32_t max_vrec_width
        cdef uintptr_t pgr_alloc_cacheline_ct
        if pgfi_init_phase2(header_ctrl, 1, 1, 0, 0, self._info_ptr[0].raw_variant_ct, &max_vrec_width, self._info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf):
            if pgfi_alloc and not self._info_ptr[0].vrtypes:
                aligned_free(pgfi_alloc)
            raise RuntimeError(errstr_buf[7:])

        self._state_ptr = <pgen_reader_t*>PyMem_Malloc(sizeof(pgen_reader_t))
        if not self._state_ptr:
            raise MemoryError()
        pgr_preinit(self._state_ptr)
        self._state_ptr[0].fread_buf = NULL
        cdef uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * kCacheline
        cdef uintptr_t sample_subset_byte_ct = DIV_UP(file_sample_ct, kBitsPerVec) * kBytesPerVec
        cdef uintptr_t cumulative_popcounts_byte_ct = DIV_UP(file_sample_ct, kBitsPerWord * kInt32PerVec) * kBytesPerVec
        cdef uintptr_t genovec_byte_ct = DIV_UP(file_sample_ct, kQuatersPerVec) * kBytesPerVec
        cdef uintptr_t dosage_vals_byte_ct = DIV_UP(file_sample_ct, (2 * kInt32PerVec)) * kBytesPerVec
        cdef unsigned char* pgr_alloc
        if cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * kPglQuaterTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + kPglQuaterTransposeBatch) * genovec_byte_ct + dosage_vals_byte_ct + kPglBitTransposeBufbytes + 4 * (kPglQuaterTransposeBatch * kPglQuaterTransposeBatch / 8), &pgr_alloc):
            raise MemoryError()
        cdef pglerr_t reterr = pgr_init(fname, max_vrec_width, self._info_ptr, self._state_ptr, pgr_alloc)
        if reterr != kPglRetSuccess:
            if not self._state_ptr[0].fread_buf:
                aligned_free(pgr_alloc)
            raise RuntimeError("pgl_init() error " + str(reterr))
        cdef unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct])
        self._subset_include_vec = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._subset_include_interleaved_vec = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        
        # assumes kWordsPerVec <= 2
        self._subset_include_interleaved_vec[-1] = 0

        self._subset_cumulative_popcounts = <uint32_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct])
        self._genovec = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct])
        self._phasepresent = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._phaseinfo = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._dosage_present = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._dosage_vals = <uint16_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[dosage_vals_byte_ct])
        if sample_subset is not None:
            self.set_sample_subset_internal(sample_subset)
        else:
            self._subset_size = file_sample_ct
        self._transpose_batch_buf = <vul_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglBitTransposeBufbytes])
        self._multivar_vmaj_geno_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglQuaterTransposeBatch * genovec_byte_ct])
        self._multivar_vmaj_phasepresent_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglQuaterTransposeBatch * sample_subset_byte_ct])
        self._multivar_vmaj_phaseinfo_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglQuaterTransposeBatch * sample_subset_byte_ct])
        self._multivar_smaj_geno_batch_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglQuaterTransposeBatch * kPglQuaterTransposeBatch / 4])
        self._multivar_smaj_phaseinfo_batch_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglQuaterTransposeBatch * kPglQuaterTransposeBatch / 8])
        self._multivar_smaj_phasepresent_batch_buf = <uintptr_t*>pgr_alloc_iter
        # pgr_alloc_iter = &(pgr_alloc_iter[kPglQuaterTransposeBatch * kPglQuaterTransposeBatch / 8])
        return

    
    cpdef __enter__(self):
        return self


    cpdef get_raw_sample_ct(self):
        return self._info_ptr[0].raw_sample_ct


    cpdef get_variant_ct(self):
        return self._info_ptr[0].raw_variant_ct


    cpdef hardcall_phase_present(self):
        return ((self._info_ptr[0].gflags & kfPgenGlobalHardcallPhasePresent) != 0)
    
    
    cpdef read(self, uint32_t variant_idx, np.ndarray geno_int_out, uint32_t allele_idx = 1):
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            # could have an unsafe mode which doesn't perform this check, but
            # let's default to at least this much bounds-checking
            raise RuntimeError("read() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        if not geno_int_out.flags["C_CONTIGUOUS"]:
            raise RuntimeError("read() requires geno_int_out to be C-contiguous.")
        # for full genotype info for multiallelic variants, use read_phased()
        # instead
        cdef pglerr_t reterr = pgr_read_allele_countvec_subset_unsafe(self._subset_include_vec, self._subset_cumulative_popcounts, self._subset_size, variant_idx, allele_idx, self._state_ptr, self._genovec)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read() error " + str(reterr))
        cdef int8_t* data8_ptr
        cdef int32_t* data32_ptr
        cdef int64_t* data64_ptr
        if geno_int_out.dtype == np.int8:
            data8_ptr = <int8_t*>geno_int_out.data
            genoarr_to_bytes_minus9(self._genovec, self._subset_size, data8_ptr)
        elif geno_int_out.dtype == np.int32:
            data32_ptr = <int32_t*>geno_int_out.data
            genoarr_to_int32s_minus9(self._genovec, self._subset_size, data32_ptr)
        elif geno_int_out.dtype == np.int64:
            data64_ptr = <int64_t*>geno_int_out.data
            genoarr_to_int64s_minus9(self._genovec, self._subset_size, data64_ptr)
        else:
            raise RuntimeError("Invalid read() geno_int_out array element type (int8, int32, or int64 expected).")
        return


    cpdef read_dosages(self, uint32_t variant_idx, np.ndarray floatarr_out, uint32_t allele_idx = 1):
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_dosages() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        if not floatarr_out.flags["C_CONTIGUOUS"]:
            raise RuntimeError("read_dosages() requires floatarr_out to be C-contiguous.")
        # todo: change this when pgenlib_internal supports multiallelic
        # variants
        cdef uint32_t dosage_ct
        cdef uint32_t is_explicit_alt1
        cdef pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(self._subset_include_vec, self._subset_cumulative_popcounts, self._subset_size, variant_idx, self._state_ptr, self._genovec, self._dosage_present, self._dosage_vals, &dosage_ct, &is_explicit_alt1)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read_dosages() error " + str(reterr))
        if allele_idx == 0:
            genovec_invert_unsafe(self._subset_size, self._genovec)
            biallelic_dosage16_invert(dosage_ct, self._dosage_vals)
        # todo: flip on allele_idx == 0
        cdef float* data32_ptr
        cdef double* data64_ptr
        if floatarr_out.dtype == np.float32:
            data32_ptr = <float*>floatarr_out.data
            dosage16_to_floats_minus9(self._genovec, self._dosage_present, self._dosage_vals, self._subset_size, dosage_ct, data32_ptr)
        elif floatarr_out.dtype == np.float64:
            data64_ptr = <double*>floatarr_out.data
            dosage16_to_doubles_minus9(self._genovec, self._dosage_present, self._dosage_vals, self._subset_size, dosage_ct, data64_ptr)
        else:
            raise RuntimeError("Invalid read_dosages() floatarr_out array element type (float32 or float64 expected).")
        return


    cpdef read_alleles(self, uint32_t variant_idx, np.ndarray[np.int32_t,mode="c",ndim=1] allele_int32_out):
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            # could have an unsafe mode which doesn't perform this check, but
            # let's default to at least this much bounds-checking
            raise RuntimeError("read_alleles() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef uint32_t phasepresent_ct
        # upgrade to multiallelic version of this function in the future
        cdef pglerr_t reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(self._subset_include_vec, self._subset_cumulative_popcounts, self._subset_size, variant_idx, self._state_ptr, self._genovec, self._phasepresent, self._phaseinfo, &phasepresent_ct)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read_alleles() error " + str(reterr))
        cdef int32_t* main_data_ptr = <int32_t*>(&(allele_int32_out[0]))
        genoarr_phased_to_allele_codes(self._genovec, self._phasepresent, self._phaseinfo, self._subset_size, phasepresent_ct, NULL, main_data_ptr)
        return


    cpdef read_alleles_and_phasepresent(self, uint32_t variant_idx, np.ndarray[np.int32_t,mode="c",ndim=1] allele_int32_out, np.ndarray[np.uint8_t,mode="c",cast=True] phasepresent_out):
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            # could have an unsafe mode which doesn't perform this check, but
            # let's default to at least this much bounds-checking
            raise RuntimeError("read_alleles_and_phasepresent() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef uint32_t phasepresent_ct
        # upgrade to multiallelic version of this function in the future
        cdef pglerr_t reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(self._subset_include_vec, self._subset_cumulative_popcounts, self._subset_size, variant_idx, self._state_ptr, self._genovec, self._phasepresent, self._phaseinfo, &phasepresent_ct)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read_alleles_and_phasepresent() error " + str(reterr))
        cdef int32_t* main_data_ptr = <int32_t*>(&(allele_int32_out[0]))
        cdef unsigned char* phasepresent_data_ptr = <unsigned char*>(&(phasepresent_out[0]))
        genoarr_phased_to_allele_codes(self._genovec, self._phasepresent, self._phaseinfo, self._subset_size, phasepresent_ct, phasepresent_data_ptr, main_data_ptr)
        return


    cdef read_range_internal8(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int8_t,mode="c",ndim=2] geno_int8_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        if variant_idx_end > self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_range() variant_idx_end too large (" + str(variant_idx_end) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int8_t* data_ptr
        cdef uint32_t variant_idx
        cdef pglerr_t reterr
        if sample_maj == 0:
            if geno_int8_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if geno_int8_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int8_out[(variant_idx - variant_idx_start), 0])
                genoarr_to_bytes_minus9(genovec, subset_size, data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if geno_int8_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int8_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* vmaj_iter
        cdef uintptr_t* smaj_iter
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, uii + variant_idx_offset, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                transpose_quaterblock(vmaj_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int8_out[uii + sample_batch_idx * kPglQuaterTransposeBatch, variant_batch_idx * kPglQuaterTransposeBatch])
                    genoarr_to_bytes_minus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglQuaterTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglQuaterTransposeWords])
            variant_idx_offset += kPglQuaterTransposeBatch
        return

    cdef read_range_internal32(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int32_t,mode="c",ndim=2] geno_int32_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        if variant_idx_end > self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_range() variant_idx_end too large (" + str(variant_idx_end) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* data_ptr
        cdef uint32_t variant_idx
        cdef pglerr_t reterr
        if sample_maj == 0:
            if geno_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if geno_int32_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = <int32_t*>(&(geno_int32_out[(variant_idx - variant_idx_start), 0]))
                genoarr_to_int32s_minus9(genovec, subset_size, data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if geno_int32_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* vmaj_iter
        cdef uintptr_t* smaj_iter
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, uii + variant_idx_offset, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                transpose_quaterblock(vmaj_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = <int32_t*>(&(geno_int32_out[uii + sample_batch_idx * kPglQuaterTransposeBatch, variant_batch_idx * kPglQuaterTransposeBatch]))
                    genoarr_to_int32s_minus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglQuaterTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglQuaterTransposeWords])
            variant_idx_offset += kPglQuaterTransposeBatch
        return

    cdef read_range_internal64(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int64_t,mode="c",ndim=2] geno_int64_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        if variant_idx_end > self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_range() variant_idx_end too large (" + str(variant_idx_end) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int64_t* data_ptr
        cdef uint32_t variant_idx
        cdef pglerr_t reterr
        if sample_maj == 0:
            if geno_int64_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if geno_int64_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int64_out[(variant_idx - variant_idx_start), 0])
                genoarr_to_int64s_minus9(genovec, subset_size, data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if geno_int64_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int64_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* vmaj_iter
        cdef uintptr_t* smaj_iter
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, uii + variant_idx_offset, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                transpose_quaterblock(vmaj_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int64_out[uii + sample_batch_idx * kPglQuaterTransposeBatch, variant_batch_idx * kPglQuaterTransposeBatch])
                    genoarr_to_int64s_minus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglQuaterTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglQuaterTransposeWords])
            variant_idx_offset += kPglQuaterTransposeBatch
        return

    cpdef read_range(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray geno_int_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        # C-contiguity checked by read_range_internal8(), etc.
        if geno_int_out.dtype == np.int8:
            self.read_range_internal8(variant_idx_start, variant_idx_end, geno_int_out, allele_idx, sample_maj)
        elif geno_int_out.dtype == np.int32:
            self.read_range_internal32(variant_idx_start, variant_idx_end, geno_int_out, allele_idx, sample_maj)
        elif geno_int_out.dtype == np.int64:
            self.read_range_internal64(variant_idx_start, variant_idx_end, geno_int_out, allele_idx, sample_maj)
        else:
            raise RuntimeError("Invalid read_range() geno_int_out array element type (int8, int32, or int64 expected).")
        return


    cdef read_list_internal8(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int8_t,mode="c",ndim=2] geno_int8_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int8_t* data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef pglerr_t reterr
        if sample_maj == 0:
            if geno_int8_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if geno_int8_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int8_out[variant_list_idx, 0])
                genoarr_to_bytes_minus9(genovec, subset_size, data_ptr)
            return
        if geno_int8_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int8_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* vmaj_iter
        cdef uintptr_t* smaj_iter
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        variant_list_idx = 0
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                variant_idx = variant_idxs[uii + variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_list() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                transpose_quaterblock(vmaj_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int8_out[uii + sample_batch_idx * kPglQuaterTransposeBatch, variant_batch_idx * kPglQuaterTransposeBatch])
                    genoarr_to_bytes_minus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglQuaterTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglQuaterTransposeWords])
            variant_list_idx += kPglQuaterTransposeBatch
        return

    cdef read_list_internal32(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int32_t,mode="c",ndim=2] geno_int32_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef pglerr_t reterr
        if sample_maj == 0:
            if geno_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if geno_int32_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = <int32_t*>(&(geno_int32_out[variant_list_idx, 0]))
                genoarr_to_int32s_minus9(genovec, subset_size, data_ptr)
            return
        if geno_int32_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* vmaj_iter
        cdef uintptr_t* smaj_iter
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        variant_list_idx = 0
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                variant_idx = variant_idxs[uii + variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_list() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                transpose_quaterblock(vmaj_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = <int32_t*>(&(geno_int32_out[uii + sample_batch_idx * kPglQuaterTransposeBatch, variant_batch_idx * kPglQuaterTransposeBatch]))
                    genoarr_to_int32s_minus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglQuaterTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglQuaterTransposeWords])
            variant_list_idx += kPglQuaterTransposeBatch
        return

    cdef read_list_internal64(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int64_t,mode="c",ndim=2] geno_int64_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int64_t* data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef pglerr_t reterr
        if sample_maj == 0:
            if geno_int64_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if geno_int64_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int64_out[variant_list_idx, 0])
                genoarr_to_int64s_minus9(genovec, subset_size, data_ptr)
            return
        if geno_int64_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int64_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* vmaj_iter
        cdef uintptr_t* smaj_iter
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        variant_list_idx = 0
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                variant_idx = variant_idxs[uii + variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = pgr_read_allele_countvec_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_list() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                transpose_quaterblock(vmaj_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int64_out[uii + sample_batch_idx * kPglQuaterTransposeBatch, variant_batch_idx * kPglQuaterTransposeBatch])
                    genoarr_to_int64s_minus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglQuaterTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglQuaterTransposeWords])
            variant_list_idx += kPglQuaterTransposeBatch
        return

    cpdef read_list(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray geno_int_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        if geno_int_out.dtype == np.int8:
            self.read_list_internal8(variant_idxs, geno_int_out, allele_idx, sample_maj)
        elif geno_int_out.dtype == np.int32:
            self.read_list_internal32(variant_idxs, geno_int_out, allele_idx, sample_maj)
        elif geno_int_out.dtype == np.int64:
            self.read_list_internal64(variant_idxs, geno_int_out, allele_idx, sample_maj)
        else:
            raise RuntimeError("Invalid read_list() geno_int_out array element type (int8, int32, or int64 expected).")
        return


    cpdef read_alleles_range(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out, bint hap_maj = 0):
        # if hap_maj == False, allele_int32_out must have at least
        #   variant_idx_ct rows, 2 * sample_ct columns
        # if hap_maj == True, allele_int32_out must have at least 2 * sample_ct
        #   rows, variant_idx_ct columns
        if variant_idx_end > self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_alleles_range() variant_idx_end too large (" + str(variant_idx_end) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* phasepresent = self._phasepresent
        cdef uintptr_t* phaseinfo = self._phaseinfo
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* main_data_ptr
        cdef uint32_t variant_idx
        cdef uint32_t phasepresent_ct
        cdef pglerr_t reterr
        if hap_maj == 0:
            if allele_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_range() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if allele_int32_out.shape[1] < 2 * subset_size:
                raise RuntimeError("Variant-major read_alleles_range() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ", and column count should be twice that)")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                # upgrade to multiallelic version of this function later
                reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, pgrp, genovec, phasepresent, phaseinfo, &phasepresent_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_alleles_range() error " + str(reterr))
                main_data_ptr = <int32_t*>(&(allele_int32_out[(variant_idx - variant_idx_start), 0]))
                genoarr_phased_to_allele_codes(genovec, phasepresent, phaseinfo, subset_size, phasepresent_ct, NULL, main_data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_alleles_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if allele_int32_out.shape[0] < 2 * subset_size:
            raise RuntimeError("Haplotype-major read_alleles_range() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ", and row count should be twice that)")
        if allele_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Haplotype-major read_alleles_range() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t variant_batch_sizel = DIV_UP(variant_batch_size, kBitsPerWord)
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_ctaw = kWordsPerVec * DIV_UP(subset_size, kBitsPerWord)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_vmaj_phaseinfo_buf = self._multivar_vmaj_phaseinfo_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* multivar_smaj_phaseinfo_batch_buf = self._multivar_smaj_phaseinfo_batch_buf
        cdef uintptr_t* vmaj_geno_iter
        cdef uintptr_t* vmaj_phaseinfo_iter
        cdef uintptr_t* smaj_geno_iter
        cdef uintptr_t* smaj_phaseinfo_iter
        cdef int32_t* main_data1_ptr
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
                variant_batch_sizel = DIV_UP(variant_batch_size, kBitsPerWord)
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for uii in range(variant_batch_size):
                reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, uii + variant_idx_offset, pgrp, vmaj_geno_iter, phasepresent, vmaj_phaseinfo_iter, &phasepresent_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_alleles_range() error " + str(reterr))
                if phasepresent_ct == 0:
                    fill_ulong_zero(sample_ctaw, vmaj_phaseinfo_iter)
                # else:
                    # bitvec_and(phasepresent, sample_ctaw, vmaj_phaseinfo_iter)
                vmaj_geno_iter = &(vmaj_geno_iter[sample_ctaw2])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[sample_ctaw])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_geno_iter = multivar_smaj_geno_batch_buf
                smaj_phaseinfo_iter = multivar_smaj_phaseinfo_batch_buf
                transpose_quaterblock(vmaj_geno_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_geno_iter, transpose_batch_buf)
                # todo: skip bitblock transpose when all phasepresent_ct values
                #       are zero, etc.
                transpose_bitblock(vmaj_phaseinfo_iter, sample_ctaw, <uint32_t>(kPglQuaterTransposeWords / 2), variant_batch_size, sample_batch_size, smaj_phaseinfo_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    main_data_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglQuaterTransposeWords), variant_batch_idx * kPglQuaterTransposeBatch]))
                    main_data1_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglQuaterTransposeWords) + 1, variant_batch_idx * kPglQuaterTransposeBatch]))
                    genoarr_phased_to_hap_codes(smaj_geno_iter, smaj_phaseinfo_iter, variant_batch_size, main_data_ptr, main_data1_ptr)
                    smaj_geno_iter = &(smaj_geno_iter[kPglQuaterTransposeWords])
                    smaj_phaseinfo_iter = &(smaj_phaseinfo_iter[kPglQuaterTransposeWords / 2])
                vmaj_geno_iter = &(vmaj_geno_iter[kPglQuaterTransposeWords])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[kPglQuaterTransposeWords / 2])
            variant_idx_offset += kPglQuaterTransposeBatch
        return


    cpdef read_alleles_list(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out, bint hap_maj = 0):
        # if hap_maj == False, allele_int32_out must have at least
        #   variant_idx_ct rows, 2 * sample_ct columns
        # if hap_maj == True, allele_int32_out must have at least 2 * sample_ct
        #   rows, variant_idx_ct columns
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef const uint32_t* subset_cumulative_popcounts = self._subset_cumulative_popcounts
        cdef pgen_reader_t* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* phasepresent = self._phasepresent
        cdef uintptr_t* phaseinfo = self._phaseinfo
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* main_data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef uint32_t phasepresent_ct
        cdef pglerr_t reterr
        if hap_maj == 0:
            if allele_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_list() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if allele_int32_out.shape[1] < 2 * subset_size:
                raise RuntimeError("Variant-major read_alleles_list() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ", and column count should be twice that)")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_alleles_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                # upgrade to multiallelic version of this function later
                reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, pgrp, genovec, phasepresent, phaseinfo, &phasepresent_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_alleles_list() error " + str(reterr))
                main_data_ptr = <int32_t*>(&(allele_int32_out[variant_list_idx, 0]))
                genoarr_phased_to_allele_codes(genovec, phasepresent, phaseinfo, subset_size, phasepresent_ct, NULL, main_data_ptr)
            return
        if allele_int32_out.shape[0] < 2 * subset_size:
            raise RuntimeError("Haplotype-major read_alleles_list() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ", and row count should be twice that)")
        if allele_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Haplotype-major read_alleles_list() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DIV_UP(variant_idx_ct, kPglQuaterTransposeBatch)
        cdef uint32_t variant_batch_size = kPglQuaterTransposeBatch
        cdef uint32_t variant_batch_sizel = DIV_UP(variant_batch_size, kBitsPerWord)
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DIV_UP(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_ctaw = kWordsPerVec * DIV_UP(subset_size, kBitsPerWord)
        cdef uint32_t sample_batch_ct = DIV_UP(subset_size, kPglQuaterTransposeBatch)
        cdef vul_t* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_vmaj_phaseinfo_buf = self._multivar_vmaj_phaseinfo_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* multivar_smaj_phaseinfo_batch_buf = self._multivar_smaj_phaseinfo_batch_buf
        cdef uintptr_t* vmaj_geno_iter
        cdef uintptr_t* vmaj_phaseinfo_iter
        cdef uintptr_t* smaj_geno_iter
        cdef uintptr_t* smaj_phaseinfo_iter
        cdef int32_t* main_data1_ptr
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t uii
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglQuaterTransposeBatch)
                variant_batch_sizel = DIV_UP(variant_batch_size, kBitsPerWord)
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for variant_list_idx in range(variant_batch_idx * kPglQuaterTransposeBatch, variant_batch_idx * kPglQuaterTransposeBatch + variant_batch_size):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_alleles_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(subset_include_vec, subset_cumulative_popcounts, subset_size, variant_idx, pgrp, vmaj_geno_iter, phasepresent, vmaj_phaseinfo_iter, &phasepresent_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_alleles_list() error " + str(reterr))
                if phasepresent_ct == 0:
                    fill_ulong_zero(sample_ctaw, vmaj_phaseinfo_iter)
                # else:
                    # bitvec_and(phasepresent, sample_ctaw, vmaj_phaseinfo_iter)
                vmaj_geno_iter = &(vmaj_geno_iter[sample_ctaw2])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[sample_ctaw])
            sample_batch_size = kPglQuaterTransposeBatch
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglQuaterTransposeBatch)
                smaj_geno_iter = multivar_smaj_geno_batch_buf
                smaj_phaseinfo_iter = multivar_smaj_phaseinfo_batch_buf
                transpose_quaterblock(vmaj_geno_iter, sample_ctaw2, kPglQuaterTransposeWords, variant_batch_size, sample_batch_size, smaj_geno_iter, transpose_batch_buf)
                # todo: skip bitblock transpose when all phasepresent_ct values
                #       are zero, etc.
                transpose_bitblock(vmaj_phaseinfo_iter, sample_ctaw, <uint32_t>(kPglQuaterTransposeWords / 2), variant_batch_size, sample_batch_size, smaj_phaseinfo_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    main_data_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglQuaterTransposeWords), variant_batch_idx * kPglQuaterTransposeBatch]))
                    main_data1_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglQuaterTransposeWords) + 1, variant_batch_idx * kPglQuaterTransposeBatch]))
                    genoarr_phased_to_hap_codes(smaj_geno_iter, smaj_phaseinfo_iter, variant_batch_size, main_data_ptr, main_data1_ptr)
                    smaj_geno_iter = &(smaj_geno_iter[kPglQuaterTransposeWords])
                    smaj_phaseinfo_iter = &(smaj_phaseinfo_iter[kPglQuaterTransposeWords / 2])
                vmaj_geno_iter = &(vmaj_geno_iter[kPglQuaterTransposeWords])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[kPglQuaterTransposeWords / 2])
        return


    cpdef read_alleles_and_phasepresent_range(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out, np.ndarray[np.uint8_t,mode="c",cast=True,ndim=2] phasepresent_out, bint hap_maj = 0):
        pass


    cpdef read_alleles_and_phasepresent_list(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out, np.ndarray[np.uint8_t,cast=True,mode="c",ndim=2] phasepresent_out, bint hap_maj = 0):
        pass

    
    cpdef count(self, uint32_t variant_idx, np.ndarray[np.uint32_t,mode="c"] genocount_uint32_out, object allele_idx = 1):
        # todo: multiallelic variants
        if allele_idx is None:
            allele_idx = 1
        cdef uint32_t* data_ptr = <uint32_t*>(&(genocount_uint32_out[0]))
        cdef pglerr_t reterr = pgr_get_refalt1_genotype_counts(self._subset_include_vec, self._subset_include_interleaved_vec, self._subset_cumulative_popcounts, self._subset_size, variant_idx, self._state_ptr, data_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("count() error " + str(reterr))
        if allele_idx != 0:
            return
        cdef uint32_t tmp = data_ptr[0]
        data_ptr[0] = data_ptr[2]
        data_ptr[2] = tmp
        return

    
    cpdef change_sample_subset(self, object sample_subset = None):
        if sample_subset is not None:
            self.set_sample_subset_internal(sample_subset)
        else:
            self._subset_size = self._info_ptr[0].raw_sample_ct
        return

    
    cpdef close(self):
        # don't bother propagating file close errors for now
        if self._info_ptr:
            pgfi_cleanup(self._info_ptr)
            if self._info_ptr[0].vrtypes:
                aligned_free(self._info_ptr[0].vrtypes)
                if self._state_ptr:
                    pgr_cleanup(self._state_ptr)
                    if self._state_ptr[0].fread_buf:
                        aligned_free(self._state_ptr[0].fread_buf)
                    PyMem_Free(self._state_ptr)
                    self._state_ptr = NULL
            PyMem_Free(self._info_ptr)
            self._info_ptr = NULL
        return

    
    cpdef __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    
    def __dealloc__(self):
        if self._info_ptr:
            pgfi_cleanup(self._info_ptr)
            if self._info_ptr[0].vrtypes:
                aligned_free(self._info_ptr[0].vrtypes)
                if self._state_ptr:
                    pgr_cleanup(self._state_ptr)
                    if self._state_ptr[0].fread_buf:
                        aligned_free(self._state_ptr[0].fread_buf)
                    PyMem_Free(self._state_ptr)
            PyMem_Free(self._info_ptr)
        return



cdef bytes_to_bits_internal(np.ndarray[np.uint8_t,mode="c",cast=True] boolbytes, uint32_t sample_ct, uintptr_t* bitarr):
    bytes_to_bits_unsafe(boolbytes, sample_ct, bitarr)

cdef class PgenWriter:
    cdef st_pgen_writer_t* _state_ptr
    cdef uintptr_t* _nonref_flags
    # preallocate buffers we'll use repeatedly
    cdef uintptr_t* _genovec
    cdef uintptr_t* _phasepresent
    cdef uintptr_t* _phaseinfo
    cdef uintptr_t* _dosage_present
    cdef uint16_t* _dosage_vals
    
    
    def __cinit__(self, bytes filename, uint32_t sample_ct,
                  uint32_t variant_ct, object nonref_flags,
                  object allele_idx_offsets = None,
                  bint hardcall_phase_present = False,
                  bint dosage_present = False,
                  bint dosage_phase_present = False):
        if dosage_phase_present and not dosage_present:
            raise RuntimeError("Invalid arguments for PgenWriter constructor (dosage_phase_present true but dosage_present false).")
        if allele_idx_offsets is not None:
            for uii in range(variant_ct + 1):
                if allele_idx_offsets[uii] != uii * 2:
                    raise RuntimeError("Multiallelic variants aren't supported by PgenWriter yet.")
                
        self._state_ptr = <st_pgen_writer_t*>PyMem_Malloc(sizeof(st_pgen_writer_t))
        if not self._state_ptr:
            raise MemoryError()
        self._nonref_flags = NULL
        cdef uint32_t nonref_flags_storage = 0
        cdef uint32_t bitvec_cacheline_ct = DIV_UP(sample_ct, kBitsPerCacheline)
        if nonref_flags is not None:
            if type(nonref_flags) == type(True):
                if nonref_flags:
                    nonref_flags_storage = 2
                else:
                    nonref_flags_storage = 1
            else:
                nonref_flags_storage = 3
                if cachealigned_malloc(bitvec_cacheline_ct * kCacheline, &(self._nonref_flags)):
                    raise MemoryError()
                bytes_to_bits_internal(nonref_flags, sample_ct, self._nonref_flags)
        cdef const char* fname = <const char*>filename
        cdef pgen_global_flags_t phase_dosage_gflags = kfPgenGlobal0
        if hardcall_phase_present:
            phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent
        if dosage_present:
            phase_dosage_gflags |= kfPgenGlobalDosagePresent
        assert not dosage_phase_present
        cdef uintptr_t alloc_cacheline_ct
        cdef uint32_t max_vrec_len
        cdef pglerr_t reterr = spgw_init_phase1(fname, NULL, self._nonref_flags, variant_ct, sample_ct, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)
        if reterr != kPglRetSuccess:
            raise RuntimeError("spgw_init_phase1() error " + str(reterr))
        cdef uint32_t genovec_cacheline_ct = DIV_UP(sample_ct, kQuatersPerCacheline)
        cdef uint32_t dosage_vals_cacheline_ct = DIV_UP(sample_ct, (2 * kInt32PerCacheline))
        cdef unsigned char* spgw_alloc
        if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct + dosage_vals_cacheline_ct) * kCacheline, &spgw_alloc):
            raise MemoryError()
        spgw_init_phase2(max_vrec_len, self._state_ptr, spgw_alloc)  
        self._genovec = <uintptr_t*>(&(spgw_alloc[alloc_cacheline_ct * kCacheline]))
        self._phasepresent = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct) * kCacheline]))
        self._phaseinfo = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + bitvec_cacheline_ct) * kCacheline]))
        self._dosage_present = <uintptr_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 2 * bitvec_cacheline_ct) * kCacheline]))
        self._dosage_vals = <uint16_t*>(&(spgw_alloc[(alloc_cacheline_ct + genovec_cacheline_ct + 3 * bitvec_cacheline_ct) * kCacheline]))
        return

    
    cpdef __enter__(self):
        return self


    cpdef append_biallelic(self, np.ndarray[np.int8_t,mode="c"] geno_int8):
        cdef int8_t* genobytes = &(geno_int8[0])
        bytes_to_genoarr_unsafe(genobytes, self._state_ptr[0].pwc.sample_ct, self._genovec)
        cdef pglerr_t reterr = spgw_append_biallelic_genovec(self._genovec, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_biallelic() error " + str(reterr))
        return
    

    cpdef append_alleles(self, np.ndarray[np.int32_t,mode="c"] allele_int32, bint all_phased = False):
        cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
        cdef uintptr_t* genovec = self._genovec
        cdef pglerr_t reterr
        if not all_phased:
            allele_codes_to_genoarr_unsafe(allele_codes, NULL, self._state_ptr[0].pwc.sample_ct, genovec, NULL, NULL)
            reterr = spgw_append_biallelic_genovec(genovec, self._state_ptr)
        else:
            allele_codes_to_genoarr_unsafe(allele_codes, NULL, self._state_ptr[0].pwc.sample_ct, genovec, self._phasepresent, self._phaseinfo)
            reterr = spgw_append_biallelic_genovec_hphase(genovec, self._phasepresent, self._phaseinfo, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_alleles() error " + str(reterr))
        return

    
    cpdef append_partially_phased(self, np.ndarray[np.int32_t,mode="c"] allele_int32, np.ndarray[np.uint8_t,cast=True] phasepresent):
        cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
        cdef unsigned char* phasepresent_bytes = <unsigned char*>(&(phasepresent[0]))
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* phasepresent_buf = self._phasepresent
        cdef uintptr_t* phaseinfo = self._phaseinfo
        allele_codes_to_genoarr_unsafe(allele_codes, phasepresent_bytes, self._state_ptr[0].pwc.sample_ct, genovec, phasepresent_buf, phaseinfo)
        cdef pglerr_t reterr = spgw_append_biallelic_genovec_hphase(genovec, phasepresent_buf, phaseinfo, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_partially_phased() error " + str(reterr))
        return


    cdef append_dosages_internal32(self, np.ndarray[np.float32_t,mode="c"] floatarr):
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_vals = self._dosage_vals
        cdef uint32_t dosage_ct
        floats_to_dosage16(<float*>(&(floatarr[0])), self._state_ptr[0].pwc.sample_ct, 6554, genovec, dosage_present, dosage_vals, &dosage_ct)
        cdef pglerr_t reterr = spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_dosages() error " + str(reterr))
        return

    cdef append_dosages_internal64(self, np.ndarray[np.float64_t,mode="c"] doublearr):
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_vals = self._dosage_vals
        cdef uint32_t dosage_ct
        doubles_to_dosage16(<double*>(&(doublearr[0])), self._state_ptr[0].pwc.sample_ct, 6554, genovec, dosage_present, dosage_vals, &dosage_ct)
        cdef pglerr_t reterr = spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_dosages() error " + str(reterr))
        return

    cpdef append_dosages(self, np.ndarray floatarr):
        if floatarr.dtype == np.float32:
            self.append_dosages_internal32(floatarr)
        elif floatarr.dtype == np.float64:
            self.append_dosages_internal64(floatarr)
        else:
            raise RuntimeError("Invalid append_dosages() dosage array element type (float32 or float64 expected).")
        return


    cpdef append_biallelic_batch(self, np.ndarray[np.int8_t,mode="c",ndim=2] geno_int8_batch):
        cdef uint32_t batch_size = <uint32_t>geno_int8_batch.shape[0]
        cdef int8_t* genobytes
        cdef uint32_t uii
        cdef pglerr_t reterr
        for uii in range(batch_size):
            genobytes = &(geno_int8_batch[uii, 0])
            bytes_to_genoarr_unsafe(genobytes, self._state_ptr[0].pwc.sample_ct, self._genovec)
            reterr = spgw_append_biallelic_genovec(self._genovec, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_biallelic_batch() error " + str(reterr))
        return
    

    cpdef append_alleles_batch(self, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_batch, bint all_phased = False):
        cdef uint32_t batch_size = <uint32_t>allele_int32_batch.shape[0]
        cdef uintptr_t* genovec = self._genovec
        cdef int32_t* allele_codes
        cdef uint32_t uii
        cdef pglerr_t reterr
        if not all_phased:
            for uii in range(batch_size):
                allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
                allele_codes_to_genoarr_unsafe(allele_codes, NULL, self._state_ptr[0].pwc.sample_ct, genovec, NULL, NULL)
                reterr = spgw_append_biallelic_genovec(genovec, self._state_ptr)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("append_alleles_batch() error " + str(reterr))
        else:
            for uii in range(batch_size):
                allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
                allele_codes_to_genoarr_unsafe(allele_codes, NULL, self._state_ptr[0].pwc.sample_ct, genovec, self._phasepresent, self._phaseinfo)
                reterr = spgw_append_biallelic_genovec_hphase(genovec, self._phasepresent, self._phaseinfo, self._state_ptr)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("append_alleles_batch() error " + str(reterr))
        return

    
    cpdef append_partially_phased_batch(self, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_batch, np.ndarray[np.uint8_t,mode="c",cast=True,ndim=2] phasepresent_batch):
        cdef uint32_t batch_size = <uint32_t>allele_int32_batch.shape[0]
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* phasepresent_buf = self._phasepresent
        cdef uintptr_t* phaseinfo = self._phaseinfo
        cdef int32_t* allele_codes
        cdef unsigned char* phasepresent_bytes
        cdef uint32_t uii
        cdef pglerr_t reterr
        for uii in range(batch_size):
            allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
            phasepresent_bytes = <unsigned char*>(&(phasepresent_batch[uii, 0]))
            allele_codes_to_genoarr_unsafe(allele_codes, phasepresent_bytes, self._state_ptr[0].pwc.sample_ct, genovec, phasepresent_buf, phaseinfo)
            reterr = spgw_append_biallelic_genovec_hphase(genovec, phasepresent_buf, phaseinfo, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_partially_phased_batch() error " + str(reterr))
        return


    cdef append_dosages_batch_internal32(self, np.ndarray[np.float32_t,mode="c",ndim=2] floatarr_batch):
        cdef uint32_t batch_size = <uint32_t>floatarr_batch.shape[0]
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_vals = self._dosage_vals
        cdef uint32_t dosage_ct
        cdef uint32_t uii
        cdef pglerr_t reterr
        for uii in range(batch_size):
            floats_to_dosage16(<float*>(&(floatarr_batch[uii, 0])), self._state_ptr[0].pwc.sample_ct, 6554, genovec, dosage_present, dosage_vals, &dosage_ct)
            reterr = spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_dosages_batch() error " + str(reterr))
        return

    cdef append_dosages_batch_internal64(self, np.ndarray[np.float64_t,mode="c",ndim=2] doublearr_batch):
        cdef uint32_t batch_size = <uint32_t>doublearr_batch.shape[0]
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_vals = self._dosage_vals
        cdef uint32_t dosage_ct
        cdef uint32_t uii
        cdef pglerr_t reterr
        for uii in range(batch_size):
            doubles_to_dosage16(<double*>(&(doublearr_batch[uii, 0])), self._state_ptr[0].pwc.sample_ct, 6554, genovec, dosage_present, dosage_vals, &dosage_ct)
            reterr = spgw_append_biallelic_genovec_dosage16(genovec, dosage_present, dosage_vals, dosage_ct, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_dosages_batch() error " + str(reterr))
        return

    cpdef append_dosages_batch(self, np.ndarray floatarr_batch):
        if floatarr_batch.dtype == np.float32:
            self.append_dosages_batch_internal32(floatarr_batch)
        elif floatarr_batch.dtype == np.float64:
            self.append_dosages_batch_internal64(floatarr_batch)
        else:
            raise RuntimeError("Invalid append_dosages_batch() dosage array element type (float32 or float64 expected).")
        return
    
    
    cpdef close(self):
        if self._state_ptr:
            if self._state_ptr[0].pwc.vidx != self._state_ptr[0].pwc.variant_ct:
                raise RuntimeError("PgenWriter.close() called when number of written variants (" + str(self._state_ptr[0].pwc.vidx) + ") unequal to initially declared value (" + str(self._state_ptr[0].pwc.variant_ct) + ").")
            spgw_finish(self._state_ptr)
            if self._nonref_flags:
                aligned_free(self._nonref_flags)
            PyMem_Free(self._state_ptr)
            self._state_ptr = NULL
        return


    cpdef __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return


    def __dealloc__(self):
        if self._state_ptr:
            if self._state_ptr[0].pwc.vidx == self._state_ptr[0].pwc.variant_ct:
                spgw_finish(self._state_ptr)
            else:
                spgw_cleanup(self._state_ptr)
            if self._nonref_flags:
                aligned_free(self._nonref_flags)
            PyMem_Free(self._state_ptr)
        return
