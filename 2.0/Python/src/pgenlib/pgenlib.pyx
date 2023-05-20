# cython: language_level=3
# from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t, uintptr_t, uint32_t, int32_t, uint16_t, uint8_t, int8_t
from libc.string cimport memcpy
from cpython.mem cimport PyMem_Malloc, PyMem_Free
# from cpython.view cimport array as cvarray
import numpy as np
cimport numpy as np
import sys

cdef extern from "../plink2/include/pgenlib_misc.h" namespace "plink2":
    ctypedef uint32_t BoolErr
    ctypedef enum PglErr:
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
        kNypsPerVec
    cdef enum:
        kNypsPerCacheline
    cdef enum:
        kBytesPerWord
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
        kPglNypTransposeBatch
    cdef enum:
        kPglNypTransposeWords
    cdef enum:
        kPglNypTransposeBufbytes
    cdef enum:
        kPglBitTransposeBufbytes

    ctypedef uint32_t PgenGlobalFlags
    cdef enum:
        kfPgenGlobal0
    cdef enum:
        kfPgenGlobalMultiallelicHardcallFound
    cdef enum:
        kfPgenGlobalHardcallPhasePresent
    cdef enum:
        kfPgenGlobalDosagePresent

    BoolErr cachealigned_malloc(uintptr_t size, void* aligned_pp)
    void aligned_free(void* aligned_ptr)

    uintptr_t RoundUpPow2(uintptr_t val, uintptr_t alignment)
    uintptr_t DivUp(uintptr_t val, uintptr_t divisor)
    void ZeroTrailingNyps(uintptr_t nyp_ct, uintptr_t* bitarr)
    void ZeroWArr(uintptr_t entry_ct, uintptr_t* ularr)
    # void BitvecAnd(const uintptr_t* arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec)
    void FillInterleavedMaskVec(const uintptr_t* subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec)
    void FillCumulativePopcounts(const uintptr_t* subset_mask, uint32_t word_ct, uint32_t* cumulative_popcounts)

    ctypedef uintptr_t VecW
    void TransposeNypblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* vecaligned_buf)
    void TransposeBitblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* vecaligned_buf)

    void GenoarrCountFreqsUnsafe(const uintptr_t* genoarr, uint32_t sample_ct, uint32_t* genocounts)
    void GenovecInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec)
    void BiallelicDosage16Invert(uint32_t dosage_ct, uint16_t* dosage_main)

    ctypedef unsigned char AlleleCode

    cdef struct PgenVariantStruct:
        uintptr_t* genovec
        uintptr_t* patch_01_set
        AlleleCode* patch_01_vals
        uintptr_t* patch_10_set
        AlleleCode* patch_10_vals
        uintptr_t* phasepresent
        uintptr_t* phaseinfo
        uintptr_t* dosage_present
        uint16_t* dosage_main
        uint32_t patch_01_ct
        uint32_t patch_10_ct
        uint32_t phasepresent_ct
    cdef enum:
        kPglMaxAlleleCt
    cdef enum:
        kPglMaxVariantCt

    uintptr_t PglComputeMaxAlleleCt(const uintptr_t* allele_idx_offsets, uint32_t variant_ct)


cdef extern from "../plink2/pvar_ffi_support.h" namespace "plink2":

    cdef struct RefcountedWptrStruct:
        uintptr_t ref_ct
        uintptr_t p[0]

    cdef struct MinimalPvarStruct:
        const char** variant_ids
        const char** allele_storage
        RefcountedWptrStruct* allele_idx_offsetsp
        uint32_t variant_ct
        uint32_t max_allele_ct

    void PreinitMinimalPvar(MinimalPvarStruct* mpp)
    PglErr LoadMinimalPvar(const char* fname, MinimalPvarStruct* mpp, char* errstr_buf)
    void CleanupMinimalPvar(MinimalPvarStruct* mpp)


cdef extern from "../plink2/pgenlib_ffi_support.h" namespace "plink2":

    void GenoarrToBytesMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int8_t* genobytes)
    void GenoarrToInt32sMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* geno_int32)
    void GenoarrToInt64sMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int64_t* geno_int64)
    void GenoarrMPToAlleleCodesMinus9(const PgenVariantStruct* pgv, uint32_t sample_ct, unsigned char* phasebytes, int32_t* allele_codes)
    void GenoarrPhasedToHapCodes(const uintptr_t* genoarr, const uintptr_t* phaseinfo, uint32_t variant_batch_size, int32_t* hap0_codes_iter, int32_t* hap1_codes_iter)
    void Dosage16ToFloatsMinus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, float* geno_float)
    void Dosage16ToDoublesMinus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double)
    void BytesToBitsUnsafe(const uint8_t* boolbytes, uint32_t sample_ct, uintptr_t* bitarr)
    void BytesToGenoarrUnsafe(const int8_t* genobytes, uint32_t sample_ct, uintptr_t* genoarr)

    int32_t ConvertMultiAlleleCodesUnsafe(const int32_t* allele_codes, const unsigned char* phasepresent_bytes, uint32_t sample_ct, uintptr_t* genoarr, uintptr_t* patch_01_set, AlleleCode* patch_01_vals, uintptr_t* patch_10_set, AlleleCode* patch_10_vals, uint32_t* patch_01_ctp, uint32_t* patch_10_ctp, uintptr_t* phasepresent, uintptr_t* phaseinfo)
    void FloatsToDosage16(const float* floatarr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr)
    void DoublesToDosage16(const double* doublearr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr)


cdef extern from "../plink2/include/pgenlib_read.h" namespace "plink2":
    cdef cppclass PgenFileInfo:
        uint32_t raw_variant_ct
        uint32_t raw_sample_ct
        unsigned char* vrtypes
        uintptr_t* allele_idx_offsets
        uintptr_t* nonref_flags
        uint32_t gflags
        uint32_t max_allele_ct

    ctypedef uint32_t PgenHeaderCtrl

    void PreinitPgfi(PgenFileInfo* pgfip)

    PglErr PgfiInitPhase1(const char* fname, const char* pgi_fname, uint32_t raw_variant_ct, uint32_t raw_sample_ct, PgenHeaderCtrl* header_ctrl_ptr, PgenFileInfo* pgfip, uintptr_t* pgfi_alloc_cacheline_ct_ptr, char* errstr_buf)

    PglErr PgfiInitPhase2(PgenHeaderCtrl header_ctrl, uint32_t allele_cts_already_loaded, uint32_t nonref_flags_already_loaded, uint32_t use_blockload, uint32_t vblock_idx_start, uint32_t vidx_end, uint32_t* max_vrec_width_ptr, PgenFileInfo* pgfip, unsigned char* pgfi_alloc, uintptr_t* pgr_alloc_cacheline_ct_ptr, char* errstr_buf)

    cdef struct PgenReaderStruct:
        pass

    unsigned char* PgrGetFreadBuf(PgenReaderStruct* pgr_ptr)

    void PgrSetFreadBuf(unsigned char* fread_buf, PgenReaderStruct* pgr_ptr)

    cdef struct PgrSampleSubsetIndexStruct:
        pass

    void PgrSetSampleSubsetIndex(const uint32_t* sample_include_cumulative_popcounts, PgenReaderStruct* pgr_ptr, PgrSampleSubsetIndexStruct* pssi_ptr)

    void PgrClearSampleSubsetIndex(PgenReaderStruct* pgr_ptr, PgrSampleSubsetIndexStruct* pssi_ptr)

    void PreinitPgr(PgenReaderStruct* pgr_ptr)

    PglErr PgrInit(const char* fname, uint32_t max_vrec_width, PgenFileInfo* pgfip, PgenReaderStruct* pgr_ptr, unsigned char* pgr_alloc)

    PglErr PgrGet1(const uintptr_t* sample_include, PgrSampleSubsetIndexStruct pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReaderStruct* pgr_ptr, uintptr_t* allele_countvec)

    PglErr PgrGetM(const uintptr_t* sample_include, PgrSampleSubsetIndexStruct pssi, uint32_t sample_ct, uint32_t vidx, PgenReaderStruct* pgr_ptr, PgenVariantStruct* pgvp)

    PglErr PgrGetP(const uintptr_t* sample_include, PgrSampleSubsetIndexStruct pssi, uint32_t sample_ct, uint32_t vidx, PgenReaderStruct* pgr_ptr, uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo, uint32_t* phasepresent_ct_ptr)

    PglErr PgrGetMP(const uintptr_t* sample_include, PgrSampleSubsetIndexStruct pssi, uint32_t sample_ct, uint32_t vidx, PgenReaderStruct* pgr_ptr, PgenVariantStruct* pgvp)

    PglErr PgrGet1D(const uintptr_t* sample_include, PgrSampleSubsetIndexStruct pssi, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReaderStruct* pgr_ptr, uintptr_t* allele_countvec, uintptr_t* dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr)

    PglErr PgrGetCounts(const uintptr_t* sample_include, const uintptr_t* sample_include_interleaved_vec, PgrSampleSubsetIndexStruct pssi, uint32_t sample_ct, uint32_t vidx, PgenReaderStruct* pgr_ptr, uint32_t* genocounts)

    BoolErr CleanupPgfi(PgenFileInfo* pgfip, PglErr* reterrp)
    BoolErr CleanupPgr(PgenReaderStruct* pgr_ptr, PglErr* reterrp)


cdef extern from "../plink2/include/pgenlib_write.h" namespace "plink2":
    cdef cppclass PgenWriterCommon:
        pass

    ctypedef enum PgenWriteMode:
        kPgenWriteBackwardSeek
        kPgenWriteAndCopy

    cdef cppclass STPgenWriter:
        pass

    uint32_t SpgwGetVariantCt(STPgenWriter* spgwp)

    uint32_t SpgwGetSampleCt(STPgenWriter* spgwp)

    uint32_t SpgwGetVidx(STPgenWriter* spgwp)

    PglErr SpgwInitPhase1(const char* fname, uintptr_t* allele_idx_offsets, uintptr_t* explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, uint32_t optional_max_allele_ct, uint32_t separate_index, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, STPgenWriter* spgwp, uintptr_t* alloc_cacheline_ct_ptr, uint32_t* max_vrec_len_ptr)

    void SpgwInitPhase2(uint32_t max_vrec_len, STPgenWriter* spgwp, unsigned char* spgw_alloc)

    PglErr SpgwAppendBiallelicGenovec(const uintptr_t* genovec, STPgenWriter* spgwp)

    PglErr SpgwAppendMultiallelicSparse(const uintptr_t* genovec, const uintptr_t* patch_01_set, const AlleleCode* patch_01_vals, const uintptr_t* patch_10_set, const AlleleCode* patch_10_vals, uint32_t allele_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, STPgenWriter* spgwp)

    PglErr SpgwAppendBiallelicGenovecHphase(const uintptr_t* genovec, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, STPgenWriter* spgwp)

    PglErr SpgwAppendMultiallelicGenovecHphase(const uintptr_t* genovec, const uintptr_t* patch_01_set, const AlleleCode* patch_01_vals, const uintptr_t* patch_10_set, const AlleleCode* patch_10_vals, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t allele_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, STPgenWriter* spgwp)

    PglErr SpgwAppendBiallelicGenovecDosage16(const uintptr_t* genovec, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, STPgenWriter* spgwp)

    PglErr SpgwFinish(STPgenWriter* spgwp)

    BoolErr CleanupSpgw(STPgenWriter* spgwp, PglErr* reterrp)


cdef class PvarReader:
    cdef MinimalPvarStruct _mp

    def __cinit__(self, bytes filename):
        PreinitMinimalPvar(&self._mp)
        cdef const char* fname = <const char*>filename
        cdef char errstr_buf[kPglErrstrBufBlen]
        if LoadMinimalPvar(fname, &self._mp, errstr_buf) != kPglRetSuccess:
            raise RuntimeError(errstr_buf[7:])
        return


    cpdef __enter__(self):
        return self


    cpdef get_variant_ct(self):
        return self._mp.variant_ct


    cpdef get_variant_id(self, uint32_t variant_idx):
        cdef uint32_t variant_ct = self._mp.variant_ct
        if variant_idx >= variant_ct:
            raise RuntimeError("0-based variant idx too large (" + str(variant_idx) + "; only " + str(variant_ct) + " in file).")
        cdef py_string = self._mp.variant_ids[variant_idx]
        return py_string


    cpdef get_allele_ct(self, uint32_t variant_idx):
        cdef uint32_t variant_ct = self._mp.variant_ct
        if variant_idx >= variant_ct:
            raise RuntimeError("0-based variant idx too large (" + str(variant_idx) + "; only " + str(variant_ct) + " in file).")
        if self._mp.allele_idx_offsetsp == NULL:
            return 2
        cdef const uintptr_t* allele_idx_offsets = self._mp.allele_idx_offsetsp.p
        return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx]


    cpdef get_allele_code(self, uint32_t variant_idx, uint32_t allele_idx):
        cdef uint32_t variant_ct = self._mp.variant_ct
        if variant_idx >= variant_ct:
            raise RuntimeError("0-based variant idx too large (" + str(variant_idx) + "; only " + str(variant_ct) + " in file).")
        cdef uintptr_t allele_idx_offset_base = 2 * variant_idx
        cdef uint32_t allele_ct = 2
        cdef const uintptr_t* allele_idx_offsets
        if self._mp.allele_idx_offsetsp != NULL:
            allele_idx_offsets = self._mp.allele_idx_offsetsp.p
            allele_idx_offset_base = allele_idx_offsets[variant_idx]
            allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offset_base
        if allele_idx >= allele_ct:
            raise RuntimeError("0-based allele idx too large (" + str(allele_idx) + "; only " + str(allele_ct) + " for this variant).")
        cdef bytes py_string = self._mp.allele_storage[allele_idx_offset_base + allele_idx]
        return py_string


    cpdef get_allele_idx_offsets(self):
        cdef uint32_t variant_ct_p1 = 1 + self._mp.variant_ct
        cdef np.ndarray[np.uintp_t,mode="c",ndim=1] allele_idx_offsets = np.zeros([variant_ct_p1], dtype=np.uintp)
        if self._mp.allele_idx_offsetsp == NULL:
            for variant_idx in range(variant_ct_p1):
                allele_idx_offsets[variant_idx] = 2 * variant_idx
        else:
            memcpy(&(allele_idx_offsets[0]), self._mp.allele_idx_offsetsp.p, variant_ct_p1 * sizeof(uintptr_t))
        return allele_idx_offsets


    cpdef get_max_allele_ct(self):
        return self._mp.max_allele_ct


    cpdef close(self):
        CleanupMinimalPvar(&self._mp)
        return


    cpdef __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return


    def __dealloc__(self):
        CleanupMinimalPvar(&self._mp)
        return


cdef class PgenReader:
    cdef PgenFileInfo* _info_ptr
    cdef PgenReaderStruct* _state_ptr
    cdef uintptr_t* _subset_include_vec
    cdef uintptr_t* _subset_include_interleaved_vec

    cdef uint32_t* _subset_cumulative_popcounts
    cdef PgrSampleSubsetIndexStruct _subset_index

    cdef uint32_t _subset_size
    # preallocate buffers we'll use repeatedly
    cdef PgenVariantStruct _pgv
    cdef VecW* _transpose_batch_buf
    # for multi-variant load-and-transpose, we load up to
    # kPglNypTransposeBatch (= 256) variants at a time, and then transpose
    cdef uintptr_t* _multivar_vmaj_geno_buf
    cdef uintptr_t* _multivar_vmaj_phasepresent_buf
    cdef uintptr_t* _multivar_vmaj_phaseinfo_buf
    cdef uintptr_t* _multivar_smaj_geno_batch_buf
    cdef uintptr_t* _multivar_smaj_phaseinfo_batch_buf
    cdef uintptr_t* _multivar_smaj_phasepresent_batch_buf

    cdef set_allele_idx_offsets_internal(self, np.ndarray[np.uintp_t,mode="c",ndim=1] allele_idx_offsets):
        # Make a copy instead of trying to share this with the caller.
        cdef uint32_t nvariant = self._info_ptr[0].raw_variant_ct
        cdef uintptr_t nbytes = (1 + nvariant) * sizeof(uintptr_t)
        if cachealigned_malloc(nbytes, &self._info_ptr[0].allele_idx_offsets):
            raise MemoryError()
        memcpy(self._info_ptr[0].allele_idx_offsets, &(allele_idx_offsets[0]), nbytes)
        self._info_ptr[0].max_allele_ct = PglComputeMaxAlleleCt(self._info_ptr[0].allele_idx_offsets, nvariant)

    cdef set_sample_subset_internal(self, np.ndarray[np.uint32_t,mode="c",ndim=1] sample_subset):
        cdef uint32_t raw_sample_ct = self._info_ptr[0].raw_sample_ct
        cdef uint32_t raw_sample_ctv = DivUp(raw_sample_ct, kBitsPerVec)
        cdef uint32_t raw_sample_ctaw = raw_sample_ctv * kWordsPerVec
        cdef uintptr_t* sample_include = self._subset_include_vec
        ZeroWArr(raw_sample_ctaw, sample_include)
        cdef uint32_t subset_size = sample_subset.size
        if subset_size == 0:
            raise RuntimeError("Empty sample_subset is not currently permitted.")
        cdef sample_uidx = sample_subset[0]
        cdef uint32_t idx = 0
        cdef next_uidx
        while True:
            if sample_uidx >= raw_sample_ct:
                raise RuntimeError("0-based sample idx too large (" + str(sample_uidx) + "; only " + str(raw_sample_ct) + " in file).")
            sample_include[sample_uidx // kBitsPerWord] |= k1LU << (sample_uidx % kBitsPerWord)
            idx += 1
            if idx == subset_size:
                break
            next_uidx = sample_subset[idx]

            # prohibit this since it implies that the caller expects genotypes
            # to be returned in a different order
            if next_uidx <= sample_uidx:
                raise RuntimeError("sample_subset is not in strictly increasing order.")

            sample_uidx = next_uidx

        FillInterleavedMaskVec(sample_include, raw_sample_ctv, self._subset_include_interleaved_vec)

        cdef uint32_t raw_sample_ctl = DivUp(raw_sample_ct, kBitsPerWord)
        # er, this isn't a safe usage pattern.
        FillCumulativePopcounts(sample_include, raw_sample_ctl, self._subset_cumulative_popcounts)
        PgrSetSampleSubsetIndex(self._subset_cumulative_popcounts, self._state_ptr, &(self._subset_index))

        self._subset_size = subset_size
        return


    def __cinit__(self, bytes filename, object raw_sample_ct = None,
                  object variant_ct = None, object sample_subset = None,
                  object allele_idx_offsets = None, object pvar = None):
        self._info_ptr = <PgenFileInfo*>PyMem_Malloc(sizeof(PgenFileInfo))
        if not self._info_ptr:
            raise MemoryError()
        PreinitPgfi(self._info_ptr)
        # this depends on pgenlib_internal implementation.  could save
        # pgfi_alloc and pgr_alloc instead.
        self._info_ptr[0].vrtypes = NULL
        cdef uint32_t cur_sample_ct = 0xffffffffU
        if raw_sample_ct is not None:
            cur_sample_ct = raw_sample_ct
        cdef uint32_t cur_variant_ct = 0xffffffffU
        cdef PvarReader pr
        if pvar is not None:
            if (variant_ct is not None) or (allele_idx_offsets is not None):
                raise RuntimeError("pvar cannot be specified at the same time as variant_ct and/or allele_idx_offsets")
            pr = pvar
            variant_ct = pr.get_variant_ct()
            if pr.get_max_allele_ct() == 2:
                allele_idx_offsets = None
            else:
                allele_idx_offsets = pr.get_allele_idx_offsets()
        if variant_ct is not None:
            cur_variant_ct = variant_ct
        cdef const char* fname = <const char*>filename
        cdef PgenHeaderCtrl header_ctrl
        cdef uintptr_t pgfi_alloc_cacheline_ct
        cdef char errstr_buf[kPglErrstrBufBlen]
        if PgfiInitPhase1(fname, NULL, cur_variant_ct, cur_sample_ct, &header_ctrl, self._info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != kPglRetSuccess:
            raise RuntimeError(errstr_buf[7:])
        cdef uint32_t file_variant_ct = self._info_ptr[0].raw_variant_ct
        if allele_idx_offsets is not None:
            self.set_allele_idx_offsets_internal(allele_idx_offsets)
        else:
            if (header_ctrl & 0x30) != 0:
                # Not tested yet, since C pgenlib doesn't yet have a mode which
                # embeds a copy of allele_idx_offsets in the .pgen.
                if cachealigned_malloc((file_variant_ct + 1) * sizeof(uintptr_t), &self._info_ptr[0].allele_idx_offsets):
                    raise MemoryError()
        cdef uintptr_t nonref_flags_byte_ct = DivUp(file_variant_ct, kBitsPerWord) * kBytesPerWord
        if (header_ctrl & 0xc0) == 0xc0:
            # explicit nonref_flags
            # might want to modify PgfiInitPhase2 to let us not load this.
            if cachealigned_malloc(nonref_flags_byte_ct, &self._info_ptr[0].nonref_flags):
                raise MemoryError()
        cdef uint32_t file_sample_ct = self._info_ptr[0].raw_sample_ct
        assert file_sample_ct != 0
        cdef unsigned char* pgfi_alloc = NULL
        if pgfi_alloc_cacheline_ct != 0:
            if cachealigned_malloc(pgfi_alloc_cacheline_ct * kCacheline, &pgfi_alloc):
                raise MemoryError()
        cdef uint32_t max_vrec_width
        cdef uintptr_t pgr_alloc_cacheline_ct
        if PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, file_variant_ct, &max_vrec_width, self._info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf):
            if pgfi_alloc and not self._info_ptr[0].vrtypes:
                aligned_free(pgfi_alloc)
            raise RuntimeError(errstr_buf[7:])
        if self._info_ptr[0].gflags & kfPgenGlobalMultiallelicHardcallFound:
            if self._info_ptr[0].allele_idx_offsets == NULL:
                raise RuntimeError("PgenReader: multiallelic variants present, but allele_idx_offsets not provided")

        self._state_ptr = <PgenReaderStruct*>PyMem_Malloc(sizeof(PgenReaderStruct))
        if not self._state_ptr:
            raise MemoryError()
        PreinitPgr(self._state_ptr)
        PgrSetFreadBuf(NULL, self._state_ptr)
        cdef uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * kCacheline
        cdef uintptr_t sample_subset_byte_ct = DivUp(file_sample_ct, kBitsPerVec) * kBytesPerVec
        cdef uintptr_t cumulative_popcounts_byte_ct = DivUp(file_sample_ct, kBitsPerWord * kInt32PerVec) * kBytesPerVec
        cdef uintptr_t genovec_byte_ct = DivUp(file_sample_ct, kNypsPerVec) * kBytesPerVec
        cdef uintptr_t patch_01_vals_byte_ct = RoundUpPow2(file_sample_ct * sizeof(AlleleCode), kBytesPerVec)
        cdef uintptr_t patch_10_vals_byte_ct = RoundUpPow2(file_sample_ct * 2 * sizeof(AlleleCode), kBytesPerVec)
        cdef uintptr_t dosage_main_byte_ct = DivUp(file_sample_ct, (2 * kInt32PerVec)) * kBytesPerVec
        cdef unsigned char* pgr_alloc
        if cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * kPglNypTransposeBatch + 7) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + kPglNypTransposeBatch) * genovec_byte_ct + patch_01_vals_byte_ct + patch_10_vals_byte_ct + dosage_main_byte_ct + kPglBitTransposeBufbytes + 4 * (kPglNypTransposeBatch * kPglNypTransposeBatch // 8), &pgr_alloc):
            raise MemoryError()
        cdef PglErr reterr = PgrInit(fname, max_vrec_width, self._info_ptr, self._state_ptr, pgr_alloc)
        if reterr != kPglRetSuccess:
            if not PgrGetFreadBuf(self._state_ptr):
                aligned_free(pgr_alloc)
            raise RuntimeError("PgrInit() error " + str(reterr))
        cdef unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct])
        self._subset_include_vec = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._subset_include_interleaved_vec = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])

        # assumes kWordsPerVec <= 2
        self._subset_include_interleaved_vec[-1] = 0

        self._subset_cumulative_popcounts = <uint32_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct])
        self._pgv.genovec = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct])
        self._pgv.patch_01_set = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._pgv.patch_01_vals = <AlleleCode*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[patch_01_vals_byte_ct])
        self._pgv.patch_10_set = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._pgv.patch_10_vals = <AlleleCode*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[patch_10_vals_byte_ct])
        self._pgv.phasepresent = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._pgv.phaseinfo = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._pgv.dosage_present = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct])
        self._pgv.dosage_main = <uint16_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct])
        if sample_subset is not None:
            self.set_sample_subset_internal(sample_subset)
        else:
            self._subset_size = file_sample_ct
        self._transpose_batch_buf = <VecW*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglBitTransposeBufbytes])
        self._multivar_vmaj_geno_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglNypTransposeBatch * genovec_byte_ct])
        self._multivar_vmaj_phasepresent_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglNypTransposeBatch * sample_subset_byte_ct])
        self._multivar_vmaj_phaseinfo_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglNypTransposeBatch * sample_subset_byte_ct])
        self._multivar_smaj_geno_batch_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglNypTransposeBatch * kPglNypTransposeBatch // 4])
        self._multivar_smaj_phaseinfo_batch_buf = <uintptr_t*>pgr_alloc_iter
        pgr_alloc_iter = &(pgr_alloc_iter[kPglNypTransposeBatch * kPglNypTransposeBatch // 8])
        self._multivar_smaj_phasepresent_batch_buf = <uintptr_t*>pgr_alloc_iter
        # pgr_alloc_iter = &(pgr_alloc_iter[kPglNypTransposeBatch * kPglNypTransposeBatch // 8])
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
        # for full genotype info for multiallelic variants, use read_alleles()
        # instead

        # when this is too much bounds-checking, caller should be using
        # read_range() or read_list()
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        if geno_int_out.ndim != 1:
            raise RuntimeError("read() requires geno_int_out to be one-dimensional.")
        cdef uint32_t subset_size = self._subset_size
        if geno_int_out.shape[0] < subset_size:
            raise RuntimeError("read() geno_int_out is too small (" + str(geno_int_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ").")

        cdef uintptr_t* genovec = self._pgv.genovec
        cdef PglErr reterr = PgrGet1(self._subset_include_vec, self._subset_index, subset_size, variant_idx, allele_idx, self._state_ptr, genovec)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read() error " + str(reterr))
        cdef int8_t* data8_ptr
        cdef int32_t* data32_ptr
        cdef int64_t* data64_ptr
        if geno_int_out.dtype == np.int8:
            data8_ptr = <int8_t*>geno_int_out.data
            GenoarrToBytesMinus9(genovec, subset_size, data8_ptr)
        elif geno_int_out.dtype == np.int32:
            data32_ptr = <int32_t*>geno_int_out.data
            GenoarrToInt32sMinus9(genovec, subset_size, data32_ptr)
        elif geno_int_out.dtype == np.int64:
            data64_ptr = <int64_t*>geno_int_out.data
            GenoarrToInt64sMinus9(genovec, subset_size, data64_ptr)
        else:
            raise RuntimeError("Invalid read() geno_int_out array element type (int8, int32, or int64 expected).")
        return


    cpdef read_dosages(self, uint32_t variant_idx, np.ndarray floatarr_out, uint32_t allele_idx = 1):
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_dosages() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        if floatarr_out.ndim != 1:
            raise RuntimeError("read_dosages() requires floatarr_out to be one-dimensional.")
        cdef uint32_t subset_size = self._subset_size
        if floatarr_out.shape[0] < subset_size:
            raise RuntimeError("read_dosages() floatarr_out is too small (" + str(floatarr_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ").")

        cdef uint32_t dosage_ct
        cdef PglErr reterr = PgrGet1D(self._subset_include_vec, self._subset_index, subset_size, variant_idx, allele_idx, self._state_ptr, self._pgv.genovec, self._pgv.dosage_present, self._pgv.dosage_main, &dosage_ct)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read_dosages() error " + str(reterr))
        cdef float* data32_ptr
        cdef double* data64_ptr
        if floatarr_out.dtype == np.float32:
            data32_ptr = <float*>floatarr_out.data
            Dosage16ToFloatsMinus9(self._pgv.genovec, self._pgv.dosage_present, self._pgv.dosage_main, subset_size, dosage_ct, data32_ptr)
        elif floatarr_out.dtype == np.float64:
            data64_ptr = <double*>floatarr_out.data
            Dosage16ToDoublesMinus9(self._pgv.genovec, self._pgv.dosage_present, self._pgv.dosage_main, subset_size, dosage_ct, data64_ptr)
        else:
            raise RuntimeError("Invalid read_dosages() floatarr_out array element type (float32 or float64 expected).")
        return


    cpdef read_alleles(self, uint32_t variant_idx, np.ndarray[np.int32_t,mode="c",ndim=1] allele_int32_out):
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_alleles() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef uint32_t subset_size = self._subset_size
        if allele_int32_out.shape[0] < subset_size:
            raise RuntimeError("read_alleles() allele_int32_out is too small (" + str(allele_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ").")

        cdef PglErr reterr = PgrGetMP(self._subset_include_vec, self._subset_index, subset_size, variant_idx, self._state_ptr, &self._pgv)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read_alleles() error " + str(reterr))
        cdef int32_t* main_data_ptr = <int32_t*>(&(allele_int32_out[0]))
        GenoarrMPToAlleleCodesMinus9(&self._pgv, subset_size, NULL, main_data_ptr)
        return


    cpdef read_alleles_and_phasepresent(self, uint32_t variant_idx, np.ndarray[np.int32_t,mode="c",ndim=1] allele_int32_out, np.ndarray[np.uint8_t,mode="c",cast=True] phasepresent_out):
        if variant_idx >= self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_alleles_and_phasepresent() variant_idx too large (" + str(variant_idx) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef uint32_t subset_size = self._subset_size
        if allele_int32_out.shape[0] < 2 * subset_size:
            raise RuntimeError("read_alleles_and_phasepresent() allele_int32_out is too small (" + str(allele_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ", , and column count should be twice that).")
        if phasepresent_out.shape[0] < subset_size:
            raise RuntimeError("read_alleles_and_phasepresent() phasepresent_out is too small (" + str(phasepresent_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ").")

        cdef PglErr reterr = PgrGetMP(self._subset_include_vec, self._subset_index, subset_size, variant_idx, self._state_ptr, &self._pgv)
        if reterr != kPglRetSuccess:
            raise RuntimeError("read_alleles_and_phasepresent() error " + str(reterr))
        cdef int32_t* main_data_ptr = <int32_t*>(&(allele_int32_out[0]))
        cdef unsigned char* phasepresent_data_ptr = <unsigned char*>(&(phasepresent_out[0]))
        GenoarrMPToAlleleCodesMinus9(&self._pgv, subset_size, phasepresent_data_ptr, main_data_ptr)
        return


    cdef read_range_internal8(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int8_t,mode="c",ndim=2] geno_int8_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        # todo: benchmark effect of adding @cython.boundscheck(False) and
        # @cython.wraparound(False) annotations to these multiple-variant-load
        # functions
        # todo: experiment with using multiple pthreads under the hood for
        # large jobs (this probably belongs under pgenlib_ffi_support, rather
        # than here)
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int8_t* data_ptr
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if sample_maj == 0:
            if geno_int8_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if geno_int8_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int8_out[(variant_idx - variant_idx_start), 0])
                GenoarrToBytesMinus9(genovec, subset_size, data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if geno_int8_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int8_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
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
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, uii + variant_idx_offset, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                TransposeNypblock(vmaj_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int8_out[uii + sample_batch_idx * kPglNypTransposeBatch, variant_batch_idx * kPglNypTransposeBatch])
                    GenoarrToBytesMinus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglNypTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglNypTransposeWords])
            variant_idx_offset += kPglNypTransposeBatch
        return

    cdef read_range_internal32(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int32_t,mode="c",ndim=2] geno_int32_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* data_ptr
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if sample_maj == 0:
            if geno_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if geno_int32_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = <int32_t*>(&(geno_int32_out[(variant_idx - variant_idx_start), 0]))
                GenoarrToInt32sMinus9(genovec, subset_size, data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if geno_int32_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
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
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, uii + variant_idx_offset, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                TransposeNypblock(vmaj_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = <int32_t*>(&(geno_int32_out[uii + sample_batch_idx * kPglNypTransposeBatch, variant_batch_idx * kPglNypTransposeBatch]))
                    GenoarrToInt32sMinus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglNypTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglNypTransposeWords])
            variant_idx_offset += kPglNypTransposeBatch
        return

    cdef read_range_internal64(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int64_t,mode="c",ndim=2] geno_int64_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int64_t* data_ptr
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if sample_maj == 0:
            if geno_int64_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if geno_int64_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_range() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int64_out[(variant_idx - variant_idx_start), 0])
                GenoarrToInt64sMinus9(genovec, subset_size, data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if geno_int64_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int64_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_range() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
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
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, uii + variant_idx_offset, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                TransposeNypblock(vmaj_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int64_out[uii + sample_batch_idx * kPglNypTransposeBatch, variant_batch_idx * kPglNypTransposeBatch])
                    GenoarrToInt64sMinus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglNypTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglNypTransposeWords])
            variant_idx_offset += kPglNypTransposeBatch
        return

    cpdef read_range(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray geno_int_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        if variant_idx_end > self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_range() variant_idx_end too large (" + str(variant_idx_end) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
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
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int8_t* data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if sample_maj == 0:
            if geno_int8_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if geno_int8_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int8_out[variant_list_idx, 0])
                GenoarrToBytesMinus9(genovec, subset_size, data_ptr)
            return
        if geno_int8_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few rows (" + str(geno_int8_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int8_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few columns (" + str(geno_int8_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
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
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                variant_idx = variant_idxs[uii + variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_list() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                TransposeNypblock(vmaj_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int8_out[uii + sample_batch_idx * kPglNypTransposeBatch, variant_batch_idx * kPglNypTransposeBatch])
                    GenoarrToBytesMinus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglNypTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglNypTransposeWords])
            variant_list_idx += kPglNypTransposeBatch
        return

    cdef read_list_internal32(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int32_t,mode="c",ndim=2] geno_int32_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if sample_maj == 0:
            if geno_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if geno_int32_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = <int32_t*>(&(geno_int32_out[variant_list_idx, 0]))
                GenoarrToInt32sMinus9(genovec, subset_size, data_ptr)
            return
        if geno_int32_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few rows (" + str(geno_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few columns (" + str(geno_int32_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
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
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                variant_idx = variant_idxs[uii + variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_list() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                TransposeNypblock(vmaj_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = <int32_t*>(&(geno_int32_out[uii + sample_batch_idx * kPglNypTransposeBatch, variant_batch_idx * kPglNypTransposeBatch]))
                    GenoarrToInt32sMinus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglNypTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglNypTransposeWords])
            variant_list_idx += kPglNypTransposeBatch
        return

    cdef read_list_internal64(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int64_t,mode="c",ndim=2] geno_int64_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int64_t* data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if sample_maj == 0:
            if geno_int64_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if geno_int64_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_list() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data_ptr = &(geno_int64_out[variant_list_idx, 0])
                GenoarrToInt64sMinus9(genovec, subset_size, data_ptr)
            return
        if geno_int64_out.shape[0] < subset_size:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few rows (" + str(geno_int64_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ")")
        if geno_int64_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Sample-major read_list() geno_int_out buffer has too few columns (" + str(geno_int64_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
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
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
            vmaj_iter = multivar_vmaj_geno_buf
            for uii in range(variant_batch_size):
                variant_idx = variant_idxs[uii + variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, vmaj_iter)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_list() error " + str(reterr))
                vmaj_iter = &(vmaj_iter[sample_ctaw2])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_iter = multivar_vmaj_geno_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_iter = multivar_smaj_geno_batch_buf
                TransposeNypblock(vmaj_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    data_ptr = &(geno_int64_out[uii + sample_batch_idx * kPglNypTransposeBatch, variant_batch_idx * kPglNypTransposeBatch])
                    GenoarrToInt64sMinus9(smaj_iter, variant_batch_size, data_ptr)
                    smaj_iter = &(smaj_iter[kPglNypTransposeWords])
                vmaj_iter = &(vmaj_iter[kPglNypTransposeWords])
            variant_list_idx += kPglNypTransposeBatch
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
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uintptr_t* phasepresent = self._pgv.phasepresent
        cdef uintptr_t* phaseinfo = self._pgv.phaseinfo
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* main_data_ptr
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if hap_maj == 0:
            if allele_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_range() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if allele_int32_out.shape[1] < 2 * subset_size:
                raise RuntimeError("Variant-major read_alleles_range() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ", and column count should be twice that)")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = PgrGetMP(subset_include_vec, subset_index, subset_size, variant_idx, pgrp, &self._pgv)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("variant_idx " + str(variant_idx) + " read_alleles_range() error " + str(reterr))
                main_data_ptr = <int32_t*>(&(allele_int32_out[(variant_idx - variant_idx_start), 0]))
                GenoarrMPToAlleleCodesMinus9(&self._pgv, subset_size, NULL, main_data_ptr)
            return
        if variant_idx_start >= variant_idx_end:
            raise RuntimeError("read_alleles_range() variant_idx_start >= variant_idx_end (" + str(variant_idx_start) + ", " + str(variant_idx_end) + ")")
        if allele_int32_out.shape[0] < 2 * subset_size:
            raise RuntimeError("Haplotype-major read_alleles_range() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ", and row count should be twice that)")
        if allele_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Haplotype-major read_alleles_range() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
        cdef uintptr_t* allele_idx_offsets = self._info_ptr[0].allele_idx_offsets
        if allele_idx_offsets != NULL:
            if allele_idx_offsets[variant_idx_end] - allele_idx_offsets[variant_idx_start] != (variant_idx_end - variant_idx_start) * k1LU * 2:
                raise RuntimeError("Haplotype-major read_alleles_range() doesn't support multiallelic variants yet.")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t variant_batch_sizel = DivUp(variant_batch_size, kBitsPerWord)
        cdef uint32_t variant_idx_offset = variant_idx_start
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_ctaw = kWordsPerVec * DivUp(subset_size, kBitsPerWord)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
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
        cdef uint32_t phasepresent_ct
        cdef uint32_t uii
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
                variant_batch_sizel = DivUp(variant_batch_size, kBitsPerWord)
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for uii in range(variant_batch_size):
                reterr = PgrGetP(subset_include_vec, subset_index, subset_size, uii + variant_idx_offset, pgrp, vmaj_geno_iter, phasepresent, vmaj_phaseinfo_iter, &phasepresent_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("variant_idx " + str(uii + variant_idx_offset) + " read_alleles_range() error " + str(reterr))
                if phasepresent_ct == 0:
                    ZeroWArr(sample_ctaw, vmaj_phaseinfo_iter)
                # else:
                    # BitvecAnd(phasepresent, sample_ctaw, vmaj_phaseinfo_iter)
                vmaj_geno_iter = &(vmaj_geno_iter[sample_ctaw2])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[sample_ctaw])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_geno_iter = multivar_smaj_geno_batch_buf
                smaj_phaseinfo_iter = multivar_smaj_phaseinfo_batch_buf
                TransposeNypblock(vmaj_geno_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_geno_iter, transpose_batch_buf)
                # todo: skip bitblock transpose when all phasepresent_ct values
                #       are zero, etc.
                TransposeBitblock(vmaj_phaseinfo_iter, sample_ctaw, <uint32_t>(kPglNypTransposeWords // 2), variant_batch_size, sample_batch_size, smaj_phaseinfo_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    # bugfix (11 Mar 2023): first index was incorrect
                    main_data_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglNypTransposeBatch), variant_batch_idx * kPglNypTransposeBatch]))
                    main_data1_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglNypTransposeBatch) + 1, variant_batch_idx * kPglNypTransposeBatch]))
                    GenoarrPhasedToHapCodes(smaj_geno_iter, smaj_phaseinfo_iter, variant_batch_size, main_data_ptr, main_data1_ptr)
                    smaj_geno_iter = &(smaj_geno_iter[kPglNypTransposeWords])
                    smaj_phaseinfo_iter = &(smaj_phaseinfo_iter[kPglNypTransposeWords // 2])
                vmaj_geno_iter = &(vmaj_geno_iter[kPglNypTransposeWords])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[kPglNypTransposeWords // 2])
            variant_idx_offset += kPglNypTransposeBatch
        return


    cpdef read_alleles_list(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out, bint hap_maj = 0):
        # if hap_maj == False, allele_int32_out must have at least
        #   variant_idx_ct rows, 2 * sample_ct columns
        # if hap_maj == True, allele_int32_out must have at least 2 * sample_ct
        #   rows, variant_idx_ct columns
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uintptr_t* phasepresent = self._pgv.phasepresent
        cdef uintptr_t* phaseinfo = self._pgv.phaseinfo
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* main_data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if hap_maj == 0:
            if allele_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_list() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if allele_int32_out.shape[1] < 2 * subset_size:
                raise RuntimeError("Variant-major read_alleles_list() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ", and column count should be twice that)")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_alleles_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGetMP(subset_include_vec, subset_index, subset_size, variant_idx, pgrp, &self._pgv)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_alleles_list() error " + str(reterr))
                main_data_ptr = <int32_t*>(&(allele_int32_out[variant_list_idx, 0]))
                GenoarrMPToAlleleCodesMinus9(&self._pgv, subset_size, NULL, main_data_ptr)
            return
        if allele_int32_out.shape[0] < 2 * subset_size:
            raise RuntimeError("Haplotype-major read_alleles_list() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; current sample subset has size " + str(subset_size) + ", and row count should be twice that)")
        if allele_int32_out.shape[1] < variant_idx_ct:
            raise RuntimeError("Haplotype-major read_alleles_list() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
        cdef uint32_t variant_batch_ct = DivUp(variant_idx_ct, kPglNypTransposeBatch)
        cdef uint32_t variant_batch_size = kPglNypTransposeBatch
        cdef uint32_t variant_batch_sizel = DivUp(variant_batch_size, kBitsPerWord)
        cdef uint32_t sample_ctaw2 = kWordsPerVec * DivUp(subset_size, kBitsPerWordD2)
        cdef uint32_t sample_ctaw = kWordsPerVec * DivUp(subset_size, kBitsPerWord)
        cdef uint32_t sample_batch_ct = DivUp(subset_size, kPglNypTransposeBatch)
        cdef VecW* transpose_batch_buf = self._transpose_batch_buf
        cdef uintptr_t* multivar_vmaj_geno_buf = self._multivar_vmaj_geno_buf
        cdef uintptr_t* multivar_vmaj_phaseinfo_buf = self._multivar_vmaj_phaseinfo_buf
        cdef uintptr_t* multivar_smaj_geno_batch_buf = self._multivar_smaj_geno_batch_buf
        cdef uintptr_t* multivar_smaj_phaseinfo_batch_buf = self._multivar_smaj_phaseinfo_batch_buf
        cdef uintptr_t* allele_idx_offsets = self._info_ptr[0].allele_idx_offsets
        cdef uintptr_t* vmaj_geno_iter
        cdef uintptr_t* vmaj_phaseinfo_iter
        cdef uintptr_t* smaj_geno_iter
        cdef uintptr_t* smaj_phaseinfo_iter
        cdef int32_t* main_data1_ptr
        cdef uint32_t variant_batch_idx
        cdef uint32_t sample_batch_size
        cdef uint32_t sample_batch_idx
        cdef uint32_t phasepresent_ct
        cdef uint32_t uii
        for variant_batch_idx in range(variant_batch_ct):
            if variant_batch_idx == (variant_batch_ct - 1):
                variant_batch_size = 1 + <uint32_t>((variant_idx_ct - 1) % kPglNypTransposeBatch)
                variant_batch_sizel = DivUp(variant_batch_size, kBitsPerWord)
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for variant_list_idx in range(variant_batch_idx * kPglNypTransposeBatch, variant_batch_idx * kPglNypTransposeBatch + variant_batch_size):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_alleles_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                if allele_idx_offsets != NULL:
                    if allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx] != 2:
                        raise RuntimeError("Haplotype-major read_alleles_list() doesn't support multiallelic variants yet.")
                reterr = PgrGetP(subset_include_vec, subset_index, subset_size, variant_idx, pgrp, vmaj_geno_iter, phasepresent, vmaj_phaseinfo_iter, &phasepresent_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_alleles_list() error " + str(reterr))
                if phasepresent_ct == 0:
                    ZeroWArr(sample_ctaw, vmaj_phaseinfo_iter)
                # else:
                    # BitvecAnd(phasepresent, sample_ctaw, vmaj_phaseinfo_iter)
                vmaj_geno_iter = &(vmaj_geno_iter[sample_ctaw2])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[sample_ctaw])
            sample_batch_size = kPglNypTransposeBatch
            vmaj_geno_iter = multivar_vmaj_geno_buf
            vmaj_phaseinfo_iter = multivar_vmaj_phaseinfo_buf
            for sample_batch_idx in range(sample_batch_ct):
                if sample_batch_idx == sample_batch_ct - 1:
                    sample_batch_size = 1 + <uint32_t>((subset_size - 1) % kPglNypTransposeBatch)
                smaj_geno_iter = multivar_smaj_geno_batch_buf
                smaj_phaseinfo_iter = multivar_smaj_phaseinfo_batch_buf
                TransposeNypblock(vmaj_geno_iter, sample_ctaw2, kPglNypTransposeWords, variant_batch_size, sample_batch_size, smaj_geno_iter, transpose_batch_buf)
                # todo: skip bitblock transpose when all phasepresent_ct values
                #       are zero, etc.
                TransposeBitblock(vmaj_phaseinfo_iter, sample_ctaw, <uint32_t>(kPglNypTransposeWords // 2), variant_batch_size, sample_batch_size, smaj_phaseinfo_iter, transpose_batch_buf)
                for uii in range(sample_batch_size):
                    main_data_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglNypTransposeBatch), variant_batch_idx * kPglNypTransposeBatch]))
                    main_data1_ptr = <int32_t*>(&(allele_int32_out[2 * (uii + sample_batch_idx * kPglNypTransposeBatch) + 1, variant_batch_idx * kPglNypTransposeBatch]))
                    GenoarrPhasedToHapCodes(smaj_geno_iter, smaj_phaseinfo_iter, variant_batch_size, main_data_ptr, main_data1_ptr)
                    smaj_geno_iter = &(smaj_geno_iter[kPglNypTransposeWords])
                    smaj_phaseinfo_iter = &(smaj_phaseinfo_iter[kPglNypTransposeWords // 2])
                vmaj_geno_iter = &(vmaj_geno_iter[kPglNypTransposeWords])
                vmaj_phaseinfo_iter = &(vmaj_phaseinfo_iter[kPglNypTransposeWords // 2])
        return


    cpdef read_alleles_and_phasepresent_range(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out, np.ndarray[np.uint8_t,mode="c",cast=True,ndim=2] phasepresent_out, bint hap_maj = 0):
        # if hap_maj == False, allele_int32_out must have at least
        #   variant_idx_ct rows, 2 * sample_ct columns
        # if hap_maj == True, allele_int32_out must have at least 2 * sample_ct
        #   rows, variant_idx_ct columns
        if variant_idx_end > self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_alleles_range() variant_idx_end too large (" + str(variant_idx_end) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* main_data_ptr
        cdef unsigned char* phasepresent_data_ptr
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if hap_maj == 0:
            if allele_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_range() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if phasepresent_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_range() phasepresent_out buffer has too few rows (" + str(phasepresent_out.shape[0]) + "; (variant_idx_end - variant_idx_start) is " + str(variant_idx_ct) + ")")
            if allele_int32_out.shape[1] < 2 * subset_size:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_range() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ", and column count should be twice that)")
            if phasepresent_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_range() phasepresent_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = PgrGetMP(subset_include_vec, subset_index, subset_size, variant_idx, pgrp, &self._pgv)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("variant_idx " + str(variant_idx) + " read_alleles_and_phasepresent_range() error " + str(reterr))
                main_data_ptr = <int32_t*>(&(allele_int32_out[(variant_idx - variant_idx_start), 0]))
                phasepresent_data_ptr = <unsigned char*>(&(phasepresent_out[(variant_idx - variant_idx_start), 0]))
                GenoarrMPToAlleleCodesMinus9(&self._pgv, subset_size, phasepresent_data_ptr, main_data_ptr)
            return
        raise RuntimeError("read_alleles_and_phasepresent_range() does not support hap_maj == 1 yet.")


    cpdef read_alleles_and_phasepresent_list(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_out, np.ndarray[np.uint8_t,cast=True,mode="c",ndim=2] phasepresent_out, bint hap_maj = 0):
        # if hap_maj == False, allele_int32_out must have at least
        #   variant_idx_ct rows, 2 * sample_ct columns
        # if hap_maj == True, allele_int32_out must have at least 2 * sample_ct
        #   rows, variant_idx_ct columns
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uintptr_t* phasepresent = self._pgv.phasepresent
        cdef uintptr_t* phaseinfo = self._pgv.phaseinfo
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef int32_t* main_data_ptr
        cdef unsigned char* phasepresent_data_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef PglErr reterr
        if hap_maj == 0:
            if allele_int32_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_list() allele_int32_out buffer has too few rows (" + str(allele_int32_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if phasepresent_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_list() phasepresent_out buffer has too few rows (" + str(phasepresent_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if allele_int32_out.shape[1] < 2 * subset_size:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_list() allele_int32_out buffer has too few columns (" + str(allele_int32_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ", and column count should be twice that)")
            if phasepresent_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_alleles_and_phasepresent_list() phasepresent_out buffer has too few columns (" + str(phasepresent_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_alleles_and_phasepresent_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGetMP(subset_include_vec, subset_index, subset_size, variant_idx, pgrp, &self._pgv)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_alleles_and_phasepresent_list() error " + str(reterr))
                main_data_ptr = <int32_t*>(&(allele_int32_out[variant_list_idx, 0]))
                phasepresent_data_ptr = <unsigned char*>(&(phasepresent_out[variant_list_idx, 0]))
                GenoarrMPToAlleleCodesMinus9(&self._pgv, subset_size, phasepresent_data_ptr, main_data_ptr)
            return
        raise RuntimeError("read_alleles_and_phasepresent_list() does not support hap_maj == 1 yet.")


    cdef read_dosages_range_internal32(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.float32_t,mode="c",ndim=2] floatarr_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uintptr_t* dosage_present = self._pgv.dosage_present
        cdef uint16_t* dosage_main = self._pgv.dosage_main
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef float* data32_ptr
        cdef uint32_t variant_idx
        cdef uint32_t dosage_ct
        cdef PglErr reterr
        if sample_maj == 0:
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = PgrGet1D(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_dosages_range() error " + str(reterr))
                data32_ptr = <float*>(&(floatarr_out[(variant_idx - variant_idx_start), 0]))
                Dosage16ToFloatsMinus9(genovec, dosage_present, dosage_main, subset_size, dosage_ct, data32_ptr)
            return
        raise RuntimeError("read_dosages_range() does not support sample_maj == 1 yet.")


    cdef read_dosages_range_internal64(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray[np.float64_t,mode="c",ndim=2] floatarr_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uintptr_t* dosage_present = self._pgv.dosage_present
        cdef uint16_t* dosage_main = self._pgv.dosage_main
        cdef uint32_t variant_idx_ct = variant_idx_end - variant_idx_start
        cdef uint32_t subset_size = self._subset_size
        cdef double* data64_ptr
        cdef uint32_t variant_idx
        cdef uint32_t dosage_ct
        cdef PglErr reterr
        if sample_maj == 0:
            for variant_idx in range(variant_idx_start, variant_idx_end):
                reterr = PgrGet1D(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_dosages_range() error " + str(reterr))
                data64_ptr = <double*>(&(floatarr_out[(variant_idx - variant_idx_start), 0]))
                Dosage16ToDoublesMinus9(genovec, dosage_present, dosage_main, subset_size, dosage_ct, data64_ptr)
            return
        raise RuntimeError("read_dosages_range() does not support sample_maj == 1 yet.")


    cpdef read_dosages_range(self, uint32_t variant_idx_start, uint32_t variant_idx_end, np.ndarray floatarr_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        if variant_idx_end > self._info_ptr[0].raw_variant_ct:
            raise RuntimeError("read_dosages_range() variant_idx_end too large (" + str(variant_idx_end) + "; only " + str(self._info_ptr[0].raw_variant_ct) + " in file)")
        if floatarr_out.dtype == np.float32:
            self.read_dosages_range_internal32(variant_idx_start, variant_idx_end, floatarr_out, allele_idx, sample_maj)
        elif floatarr_out.dtype == np.float64:
            self.read_dosages_range_internal64(variant_idx_start, variant_idx_end, floatarr_out, allele_idx, sample_maj)
        else:
            raise RuntimeError("Invalid read_dosages_range() floatarr_out array element type (float32 or float64 expected).")
        return


    cdef read_dosages_list_internal32(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.float32_t,mode="c",ndim=2] floatarr_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uintptr_t* dosage_present = self._pgv.dosage_present
        cdef uint16_t* dosage_main = self._pgv.dosage_main
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef float* data32_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef uint32_t dosage_ct
        cdef PglErr reterr
        if sample_maj == 0:
            if floatarr_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_dosages_list() floatarr_out buffer has too few rows (" + str(floatarr_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if floatarr_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_dosages_list() floatarr_out buffer has too few columns (" + str(floatarr_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_dosages_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1D(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data32_ptr = <float*>(&(floatarr_out[variant_list_idx, 0]))
                Dosage16ToFloatsMinus9(genovec, dosage_present, dosage_main, subset_size, dosage_ct, data32_ptr)
            return
        raise RuntimeError("read_dosages_list() does not support sample_maj == 1 yet.")


    cdef read_dosages_list_internal64(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray[np.float64_t,mode="c",ndim=2] floatarr_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        cdef uint32_t raw_variant_ct = self._info_ptr[0].raw_variant_ct
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uintptr_t* genovec = self._pgv.genovec
        cdef uintptr_t* dosage_present = self._pgv.dosage_present
        cdef uint16_t* dosage_main = self._pgv.dosage_main
        cdef uint32_t variant_idx_ct = <uint32_t>variant_idxs.shape[0]
        cdef uint32_t subset_size = self._subset_size
        cdef double* data64_ptr
        cdef uint32_t variant_list_idx
        cdef uint32_t variant_idx
        cdef uint32_t dosage_ct
        cdef PglErr reterr
        if sample_maj == 0:
            if floatarr_out.shape[0] < variant_idx_ct:
                raise RuntimeError("Variant-major read_dosages_list() floatarr_out buffer has too few rows (" + str(floatarr_out.shape[0]) + "; variant_idxs length is " + str(variant_idx_ct) + ")")
            if floatarr_out.shape[1] < subset_size:
                raise RuntimeError("Variant-major read_dosages_list() floatarr_out buffer has too few columns (" + str(floatarr_out.shape[1]) + "; current sample subset has size " + str(subset_size) + ")")
            for variant_list_idx in range(variant_idx_ct):
                variant_idx = variant_idxs[variant_list_idx]
                if variant_idx >= raw_variant_ct:
                    raise RuntimeError("read_dosages_list() variant index too large (" + str(variant_idx) + "; only " + str(raw_variant_ct) + " in file)")
                reterr = PgrGet1D(subset_include_vec, subset_index, subset_size, variant_idx, allele_idx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("read_range() error " + str(reterr))
                data64_ptr = <double*>(&(floatarr_out[variant_list_idx, 0]))
                Dosage16ToDoublesMinus9(genovec, dosage_present, dosage_main, subset_size, dosage_ct, data64_ptr)
            return
        raise RuntimeError("read_dosages_list() does not support sample_maj == 1 yet.")


    cpdef read_dosages_list(self, np.ndarray[np.uint32_t] variant_idxs, np.ndarray floatarr_out, uint32_t allele_idx = 1, bint sample_maj = 0):
        if floatarr_out.dtype == np.float32:
            self.read_dosages_list_internal32(variant_idxs, floatarr_out, allele_idx, sample_maj)
        elif floatarr_out.dtype == np.float64:
            self.read_dosages_list_internal64(variant_idxs, floatarr_out, allele_idx, sample_maj)
        else:
            raise RuntimeError("Invalid read_dosages_list() floatarr_out array element type (float32 or float64 expected).")
        return

    cpdef count(self, uint32_t variant_idx, np.ndarray[np.uint32_t,mode="c"] genocount_uint32_out, object allele_idx = 1):
        cdef const uintptr_t* subset_include_vec = self._subset_include_vec
        cdef PgrSampleSubsetIndexStruct subset_index = self._subset_index
        cdef uint32_t subset_size = self._subset_size
        cdef PgenReaderStruct* pgrp = self._state_ptr
        cdef uint32_t* data_ptr = <uint32_t*>(&(genocount_uint32_out[0]))
        cdef uintptr_t* allele_idx_offsets
        cdef uintptr_t* genovec
        cdef AlleleCode* patch_01_vals
        cdef AlleleCode* patch_10_vals
        cdef uintptr_t type_ct
        cdef uintptr_t allele_ct
        cdef uint32_t patch_01_ct
        cdef uint32_t patch_10_ct
        cdef uintptr_t ac0
        cdef uintptr_t ac1
        cdef PglErr reterr
        if allele_idx is None:
            allele_ct = 2
            allele_idx_offsets = self._info_ptr[0].allele_idx_offsets
            if allele_idx_offsets != NULL:
                allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx]
            if allele_ct > 2:
                # Multiallelic special case.
                type_ct = (allele_ct * (allele_ct + 1)) >> 1
                if <uintptr_t>(genocount_uint32_out.shape[0]) <= type_ct:
                    raise RuntimeError("count() genocount_uint32_out buffer is too small for multiallelic variant")
                reterr = PgrGetM(subset_include_vec, subset_index, subset_size, variant_idx, pgrp, &self._pgv)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("count() error " + str(reterr))
                genovec = self._pgv.genovec
                ZeroTrailingNyps(subset_size, genovec)
                GenoarrCountFreqsUnsafe(genovec, subset_size, data_ptr)
                data_ptr[type_ct] = data_ptr[3]
                for i in range(3, type_ct):
                    data_ptr[i] = 0
                # Now correct for patch_01 and patch_10.
                patch_01_ct = self._pgv.patch_01_ct
                if patch_01_ct > 0:
                    data_ptr[1] -= patch_01_ct
                    patch_01_vals = self._pgv.patch_01_vals
                    # could special-case allele_ct == 3, etc.
                    for i in range(patch_01_ct):
                        ac1 = patch_01_vals[i]
                        data_ptr[(ac1 * (ac1 + 1)) >> 1] += 1
                patch_10_ct = self._pgv.patch_10_ct
                if patch_10_ct > 0:
                    data_ptr[2] -= patch_10_ct
                    patch_10_vals = self._pgv.patch_10_vals
                    for i in range(patch_10_ct):
                        ac0 = patch_10_vals[2*i]
                        ac1 = patch_10_vals[2*i+1]
                        data_ptr[ac0 + ((ac1 * (ac1 + 1)) >> 1)] += 1
                return
            allele_idx = 1
        reterr = PgrGetCounts(subset_include_vec, self._subset_include_interleaved_vec, subset_index, subset_size, variant_idx, pgrp, data_ptr)
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
            PgrClearSampleSubsetIndex(self._state_ptr, &(self._subset_index))
            self._subset_size = self._info_ptr[0].raw_sample_ct
        return


    cpdef close(self):
        cdef PglErr reterr = kPglRetSuccess
        if self._info_ptr:
            CleanupPgfi(self._info_ptr, &reterr)
            if self._info_ptr[0].vrtypes:
                aligned_free(self._info_ptr[0].vrtypes)
            if self._info_ptr[0].allele_idx_offsets:
                aligned_free(self._info_ptr[0].allele_idx_offsets)
            if self._info_ptr[0].nonref_flags:
                aligned_free(self._info_ptr[0].nonref_flags)
            if self._state_ptr:
                CleanupPgr(self._state_ptr, &reterr)
                if PgrGetFreadBuf(self._state_ptr):
                    aligned_free(PgrGetFreadBuf(self._state_ptr))
                PyMem_Free(self._state_ptr)
                self._state_ptr = NULL
            PyMem_Free(self._info_ptr)
            self._info_ptr = NULL
            if reterr != kPglRetSuccess:
                raise RuntimeError("close() error " + str(reterr))
        return


    cpdef __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return


    def __dealloc__(self):
        cdef PglErr reterr = kPglRetSuccess
        if self._info_ptr:
            CleanupPgfi(self._info_ptr, &reterr)
            if self._info_ptr[0].vrtypes:
                aligned_free(self._info_ptr[0].vrtypes)
            if self._info_ptr[0].allele_idx_offsets:
                aligned_free(self._info_ptr[0].allele_idx_offsets)
            if self._info_ptr[0].nonref_flags:
                aligned_free(self._info_ptr[0].nonref_flags)
            if self._state_ptr:
                CleanupPgr(self._state_ptr, &reterr)
                if PgrGetFreadBuf(self._state_ptr):
                    aligned_free(PgrGetFreadBuf(self._state_ptr))
                PyMem_Free(self._state_ptr)
            PyMem_Free(self._info_ptr)
        return


cdef bytes_to_bits_internal(np.ndarray[np.uint8_t,mode="c",cast=True] boolbytes, uint32_t sample_ct, uintptr_t* bitarr):
    BytesToBitsUnsafe(boolbytes, sample_ct, bitarr)


cdef class PgenWriter:
    cdef STPgenWriter* _state_ptr
    cdef uintptr_t* _nonref_flags
    cdef PgenGlobalFlags _phase_dosage_gflags
    cdef uint32_t _allele_ct_limit
    # preallocate buffers we'll use repeatedly
    cdef uintptr_t* _genovec
    cdef uintptr_t* _patch_01_set
    cdef AlleleCode* _patch_01_vals
    cdef uintptr_t* _patch_10_set
    cdef AlleleCode* _patch_10_vals
    cdef uintptr_t* _phasepresent
    cdef uintptr_t* _phaseinfo
    cdef uintptr_t* _dosage_present
    cdef uint16_t* _dosage_main

    def __cinit__(self, bytes filename, uint32_t sample_ct,
                  object variant_ct = None, object nonref_flags = True,
                  uint32_t allele_ct_limit = 2,
                  bint hardcall_phase_present = False,
                  bint dosage_present = False,
                  bint dosage_phase_present = False,
                  object variant_ct_limit = None):
        cdef PgenWriteMode write_mode = kPgenWriteBackwardSeek
        cdef uint32_t cur_variant_ct_limit
        if variant_ct is not None:
            # Could also enforce variant_ct >= variant_ct_limit when both are
            # provided?
            cur_variant_ct_limit = variant_ct
        elif variant_ct_limit is not None:
            write_mode = kPgenWriteAndCopy
            cur_variant_ct_limit = variant_ct_limit
        else:
            raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct or variant_ct_upper_bound must be provided).")
        if cur_variant_ct_limit == 0 or cur_variant_ct_limit > kPglMaxVariantCt:
            raise RuntimeError("Invalid arguments for PgenWriter constructor (variant_ct must be positive, and less than ~2^31).")
        if dosage_phase_present and not dosage_present:
            raise RuntimeError("Invalid arguments for PgenWriter constructor (dosage_phase_present true but dosage_present false).")
        if allele_ct_limit != 2:
            if (allele_ct_limit < 2) or (allele_ct_limit > kPglMaxAlleleCt):
                raise RuntimeError("Invalid arguments for PgenWriter constructor (allele_ct_limit must be in [2, 255]).")
            write_mode = kPgenWriteAndCopy

        self._state_ptr = <STPgenWriter*>PyMem_Malloc(sizeof(STPgenWriter))
        if not self._state_ptr:
            raise MemoryError()
        self._nonref_flags = NULL
        cdef uint32_t nonref_flags_storage = 0
        cdef uint32_t bitvec_cacheline_ct = DivUp(sample_ct, kBitsPerCacheline)
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
        cdef PgenGlobalFlags phase_dosage_gflags = kfPgenGlobal0
        if hardcall_phase_present:
            phase_dosage_gflags |= kfPgenGlobalHardcallPhasePresent
        if dosage_present:
            phase_dosage_gflags |= kfPgenGlobalDosagePresent
        self._phase_dosage_gflags = phase_dosage_gflags
        assert not dosage_phase_present

        cdef uintptr_t alloc_cacheline_ct
        cdef uint32_t max_vrec_len
        cdef PglErr reterr = SpgwInitPhase1(fname, NULL, self._nonref_flags, cur_variant_ct_limit, sample_ct, allele_ct_limit, write_mode, phase_dosage_gflags, nonref_flags_storage, self._state_ptr, &alloc_cacheline_ct, &max_vrec_len)
        if reterr != kPglRetSuccess:
            raise RuntimeError("SpgwInitPhase1() error " + str(reterr))
        cdef uint32_t genovec_cacheline_ct = DivUp(sample_ct, kNypsPerCacheline)
        cdef uint32_t patch_01_vals_cacheline_ct = DivUp(sample_ct * sizeof(AlleleCode), kCacheline)
        cdef uint32_t patch_10_vals_cacheline_ct = DivUp(sample_ct * 2 * sizeof(AlleleCode), kCacheline)
        cdef uint32_t dosage_main_cacheline_ct = DivUp(sample_ct, (2 * kInt32PerCacheline))
        cdef unsigned char* spgw_alloc
        if cachealigned_malloc((alloc_cacheline_ct + genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct) * kCacheline, &spgw_alloc):
            raise MemoryError()
        SpgwInitPhase2(max_vrec_len, self._state_ptr, spgw_alloc)
        cdef unsigned char* spgw_alloc_iter = &(spgw_alloc[alloc_cacheline_ct * kCacheline])
        self._allele_ct_limit = allele_ct_limit
        self._genovec = <uintptr_t*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * kCacheline])

        # Can't skimp on patch_{01,10}_{set,vals} allocations even when
        # allele_ct_limit == 2, due to how ConvertMultiAlleleCodesUnsafe()
        # works.
        # Could skimp on dosage/phase, but that doesn't gain us much.
        self._patch_01_set = <uintptr_t*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        self._patch_01_vals = <AlleleCode*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * kCacheline])
        self._patch_10_set = <uintptr_t*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        self._patch_10_vals = <AlleleCode*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[patch_10_vals_cacheline_ct * kCacheline])
        self._phasepresent = <uintptr_t*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        self._phaseinfo = <uintptr_t*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        self._dosage_present = <uintptr_t*>(spgw_alloc_iter)
        spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * kCacheline])
        self._dosage_main = <uint16_t*>(spgw_alloc_iter)
        # bugfix (16 Apr 2023): SpgwAppendBiallelicGenovec[Hphase] assumes
        # trailing bits are clear
        self._genovec[(sample_ct - 1) // kBitsPerWordD2] = 0
        self._phasepresent[(sample_ct - 1) // kBitsPerWord] = 0
        return


    cpdef __enter__(self):
        return self


    cpdef append_biallelic(self, np.ndarray[np.int8_t,mode="c"] geno_int8):
        cdef int8_t* genobytes = &(geno_int8[0])
        BytesToGenoarrUnsafe(genobytes, SpgwGetSampleCt(self._state_ptr), self._genovec)
        cdef PglErr reterr = SpgwAppendBiallelicGenovec(self._genovec, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_biallelic() error " + str(reterr))
        return


    cpdef append_alleles(self, np.ndarray[np.int32_t,mode="c"] allele_int32, bint all_phased = False, object allele_ct = None):
        cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
        cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
        cdef uint32_t allele_ct_limit = self._allele_ct_limit
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* patch_01_set = self._patch_01_set
        cdef AlleleCode* patch_01_vals = self._patch_01_vals
        cdef uintptr_t* patch_10_set = self._patch_10_set
        cdef AlleleCode* patch_10_vals = self._patch_10_vals
        cdef uintptr_t* phaseinfo = NULL
        if all_phased:
            if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
                raise RuntimeError("append_alleles called with all_phased True, but PgenWriter was constructed with hardcall_phase_present False")
            phaseinfo = self._phaseinfo
        cdef uint32_t patch_01_ct
        cdef uint32_t patch_10_ct
        cdef int32_t observed_allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, NULL, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, NULL, phaseinfo)
        if observed_allele_ct == -1:
            raise RuntimeError("append_alleles called with invalid allele codes")
        cdef uint32_t write_allele_ct = <uint32_t>(observed_allele_ct)
        if write_allele_ct > allele_ct_limit:
            raise RuntimeError("append_alleles called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
        if allele_ct is not None:
            if allele_ct < write_allele_ct:
                raise RuntimeError("append_alleles called with allele codes >= allele_ct argument")
            if allele_ct > allele_ct_limit:
                raise RuntimeError("append_alleles called with allele_ct > allele_ct_limit")
            write_allele_ct = allele_ct
        cdef PglErr reterr
        if not all_phased:
            if (patch_01_ct == 0) and (patch_10_ct == 0):
                reterr = SpgwAppendBiallelicGenovec(genovec, self._state_ptr)
            else:
                reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
        else:
            if (patch_01_ct == 0) and (patch_10_ct == 0):
                reterr = SpgwAppendBiallelicGenovecHphase(genovec, NULL, phaseinfo, self._state_ptr)
            else:
                reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, NULL, phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_alleles() error " + str(reterr))
        return


    cpdef append_partially_phased(self, np.ndarray[np.int32_t,mode="c"] allele_int32, np.ndarray[np.uint8_t,cast=True] phasepresent, object allele_ct = None):
        if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
            raise RuntimeError("append_partially_phased cannot be called when PgenWriter was constructed with hardcall_phase_present False")
        cdef int32_t* allele_codes = <int32_t*>(&(allele_int32[0]))
        cdef unsigned char* phasepresent_bytes = <unsigned char*>(&(phasepresent[0]))
        cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
        cdef uint32_t allele_ct_limit = self._allele_ct_limit
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* phasepresent_buf = self._phasepresent
        cdef uintptr_t* phaseinfo = self._phaseinfo
        cdef uintptr_t* patch_01_set = self._patch_01_set
        cdef AlleleCode* patch_01_vals = self._patch_01_vals
        cdef uintptr_t* patch_10_set = self._patch_10_set
        cdef AlleleCode* patch_10_vals = self._patch_10_vals
        cdef uint32_t patch_01_ct
        cdef uint32_t patch_10_ct
        cdef int32_t observed_allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, phasepresent_bytes, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, phasepresent_buf, phaseinfo)
        if observed_allele_ct == -1:
            raise RuntimeError("append_partially_phased called with invalid allele codes")
        cdef uint32_t write_allele_ct = <uint32_t>(observed_allele_ct)
        if write_allele_ct > allele_ct_limit:
            raise RuntimeError("append_partially_phased called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
        if allele_ct is not None:
            if allele_ct < write_allele_ct:
                raise RuntimeError("append_partially_phased called with allele codes >= allele_ct argument")
            if allele_ct > allele_ct_limit:
                raise RuntimeError("append_partially_phased called with allele_ct > allele_ct_limit")
            write_allele_ct = allele_ct
        cdef PglErr reterr
        if (patch_01_ct == 0) and (patch_10_ct == 0):
            reterr = SpgwAppendBiallelicGenovecHphase(genovec, phasepresent_buf, phaseinfo, self._state_ptr)
        else:
            reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent_buf, phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_partially_phased() error " + str(reterr))
        return


    cdef append_dosages_internal32(self, np.ndarray[np.float32_t,mode="c"] floatarr):
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_main = self._dosage_main
        cdef uint32_t dosage_ct
        FloatsToDosage16(<float*>(&(floatarr[0])), SpgwGetSampleCt(self._state_ptr), 6554, genovec, dosage_present, dosage_main, &dosage_ct)
        cdef PglErr reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_dosages() error " + str(reterr))
        return

    cdef append_dosages_internal64(self, np.ndarray[np.float64_t,mode="c"] doublearr):
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_main = self._dosage_main
        cdef uint32_t dosage_ct
        DoublesToDosage16(<double*>(&(doublearr[0])), SpgwGetSampleCt(self._state_ptr), 6554, genovec, dosage_present, dosage_main, &dosage_ct)
        cdef PglErr reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, self._state_ptr)
        if reterr != kPglRetSuccess:
            raise RuntimeError("append_dosages() error " + str(reterr))
        return

    cpdef append_dosages(self, np.ndarray floatarr):
        if (self._phase_dosage_gflags & kfPgenGlobalDosagePresent) == 0:
            raise RuntimeError("append_dosages cannot be called when PgenWriter was constructed with dosage_present False")
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
        cdef PglErr reterr
        for uii in range(batch_size):
            genobytes = &(geno_int8_batch[uii, 0])
            BytesToGenoarrUnsafe(genobytes, SpgwGetSampleCt(self._state_ptr), self._genovec)
            reterr = SpgwAppendBiallelicGenovec(self._genovec, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_biallelic_batch() error " + str(reterr))
        return


    cpdef append_alleles_batch(self, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_batch, bint all_phased = False, object allele_cts = None):
        cdef uint32_t batch_size = <uint32_t>allele_int32_batch.shape[0]
        cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
        cdef uint32_t allele_ct_limit = self._allele_ct_limit
        cdef uintptr_t* genovec = self._genovec
        cdef int32_t* allele_codes
        cdef uintptr_t* patch_01_set = self._patch_01_set
        cdef AlleleCode* patch_01_vals = self._patch_01_vals
        cdef uintptr_t* patch_10_set = self._patch_10_set
        cdef AlleleCode* patch_10_vals = self._patch_10_vals
        cdef uint32_t uii
        cdef int32_t observed_allele_ct
        cdef uint32_t write_allele_ct
        cdef uint32_t casted_allele_ct
        cdef uint32_t patch_01_ct
        cdef uint32_t patch_10_ct
        cdef PglErr reterr
        if not all_phased:
            for uii in range(batch_size):
                allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
                observed_allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, NULL, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, NULL, NULL)
                if observed_allele_ct == -1:
                    raise RuntimeError("append_alleles_batch called with invalid allele codes")
                write_allele_ct = <uint32_t>(observed_allele_ct)
                if write_allele_ct > allele_ct_limit:
                    raise RuntimeError("append_alleles_batch called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
                if allele_cts is not None:
                    casted_allele_ct = allele_cts[uii]
                    if casted_allele_ct < write_allele_ct:
                        raise RuntimeError("append_alleles_batch called with allele codes >= allele_cts[] value")
                    if casted_allele_ct > allele_ct_limit:
                        raise RuntimeError("append_alleles_batch called with allele_cts[] value > allele_ct_limit")
                    write_allele_ct = casted_allele_ct
                if (patch_01_ct == 0) and (patch_10_ct == 0):
                    reterr = SpgwAppendBiallelicGenovec(genovec, self._state_ptr)
                else:
                    reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("append_alleles_batch() error " + str(reterr))
        else:
            if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
                raise RuntimeError("append_alleles_batch called with all_phased True, but PgenWriter was constructed with hardcall_phase_present False")
            for uii in range(batch_size):
                allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
                observed_allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, NULL, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, NULL, self._phaseinfo)
                if observed_allele_ct == -1:
                    raise RuntimeError("append_alleles_batch called with invalid allele codes")
                write_allele_ct = <uint32_t>(observed_allele_ct)
                if write_allele_ct > allele_ct_limit:
                    raise RuntimeError("append_alleles_batch called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
                if allele_cts is not None:
                    casted_allele_ct = allele_cts[uii]
                    if casted_allele_ct < write_allele_ct:
                        raise RuntimeError("append_alleles_batch called with allele codes >= allele_cts[] value")
                    if casted_allele_ct > allele_ct_limit:
                        raise RuntimeError("append_alleles_batch called with allele_cts[] value > allele_ct_limit")
                    write_allele_ct = casted_allele_ct
                if (patch_01_ct == 0) and (patch_10_ct == 0):
                    reterr = SpgwAppendBiallelicGenovecHphase(genovec, NULL, self._phaseinfo, self._state_ptr)
                else:
                    reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, NULL, self._phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("append_alleles_batch() error " + str(reterr))
        return


    cpdef append_partially_phased_batch(self, np.ndarray[np.int32_t,mode="c",ndim=2] allele_int32_batch, np.ndarray[np.uint8_t,mode="c",cast=True,ndim=2] phasepresent_batch, object allele_cts = None):
        if (self._phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) == 0:
            raise RuntimeError("append_partially_phased_batch cannot be called when PgenWriter was constructed with hardcall_phase_present False")
        cdef uint32_t batch_size = <uint32_t>allele_int32_batch.shape[0]
        cdef uint32_t sample_ct = SpgwGetSampleCt(self._state_ptr)
        cdef uint32_t allele_ct_limit = self._allele_ct_limit
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* patch_01_set = self._patch_01_set
        cdef AlleleCode* patch_01_vals = self._patch_01_vals
        cdef uintptr_t* patch_10_set = self._patch_10_set
        cdef AlleleCode* patch_10_vals = self._patch_10_vals
        cdef uintptr_t* phasepresent_buf = self._phasepresent
        cdef uintptr_t* phaseinfo = self._phaseinfo
        cdef int32_t* allele_codes
        cdef unsigned char* phasepresent_bytes
        cdef uint32_t uii
        cdef int32_t observed_allele_ct
        cdef uint32_t write_allele_ct
        cdef uint32_t casted_allele_ct
        cdef uint32_t patch_01_ct
        cdef uint32_t patch_10_ct
        cdef PglErr reterr
        for uii in range(batch_size):
            allele_codes = <int32_t*>(&(allele_int32_batch[uii, 0]))
            phasepresent_bytes = <unsigned char*>(&(phasepresent_batch[uii, 0]))
            observed_allele_ct = ConvertMultiAlleleCodesUnsafe(allele_codes, phasepresent_bytes, sample_ct, genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, &patch_01_ct, &patch_10_ct, phasepresent_buf, phaseinfo)
            if observed_allele_ct == -1:
                raise RuntimeError("append_partially_phased_batch called with invalid allele codes")
            write_allele_ct = <uint32_t>(observed_allele_ct)
            if write_allele_ct > allele_ct_limit:
                raise RuntimeError("append_partially_phased_batch called with allele codes >= allele_ct_limit; you may need to construct the PgenWriter with a higher allele_ct_limit setting")
            if allele_cts is not None:
                casted_allele_ct = allele_cts[uii]
                if casted_allele_ct < write_allele_ct:
                    raise RuntimeError("append_alleles_batch called with allele codes >= allele_cts[] value")
                if casted_allele_ct > allele_ct_limit:
                    raise RuntimeError("append_alleles_batch called with allele_cts[] value > allele_ct_limit")
                write_allele_ct = casted_allele_ct
            if (patch_01_ct == 0) and (patch_10_ct == 0):
                reterr = SpgwAppendBiallelicGenovecHphase(genovec, phasepresent_buf, phaseinfo, self._state_ptr)
            else:
                reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent_buf, phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_partially_phased_batch() error " + str(reterr))
        return


    cdef append_dosages_batch_internal32(self, np.ndarray[np.float32_t,mode="c",ndim=2] floatarr_batch):
        cdef uint32_t batch_size = <uint32_t>floatarr_batch.shape[0]
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_main = self._dosage_main
        cdef uint32_t dosage_ct
        cdef uint32_t uii
        cdef PglErr reterr
        for uii in range(batch_size):
            FloatsToDosage16(<float*>(&(floatarr_batch[uii, 0])), SpgwGetSampleCt(self._state_ptr), 6554, genovec, dosage_present, dosage_main, &dosage_ct)
            reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_dosages_batch() error " + str(reterr))
        return

    cdef append_dosages_batch_internal64(self, np.ndarray[np.float64_t,mode="c",ndim=2] doublearr_batch):
        cdef uint32_t batch_size = <uint32_t>doublearr_batch.shape[0]
        cdef uintptr_t* genovec = self._genovec
        cdef uintptr_t* dosage_present = self._dosage_present
        cdef uint16_t* dosage_main = self._dosage_main
        cdef uint32_t dosage_ct
        cdef uint32_t uii
        cdef PglErr reterr
        for uii in range(batch_size):
            DoublesToDosage16(<double*>(&(doublearr_batch[uii, 0])), SpgwGetSampleCt(self._state_ptr), 6554, genovec, dosage_present, dosage_main, &dosage_ct)
            reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("append_dosages_batch() error " + str(reterr))
        return

    cpdef append_dosages_batch(self, np.ndarray floatarr_batch):
        if (self._phase_dosage_gflags & kfPgenGlobalDosagePresent) == 0:
            raise RuntimeError("append_dosages_batch cannot be called when PgenWriter was constructed with dosage_present False")
        if floatarr_batch.dtype == np.float32:
            self.append_dosages_batch_internal32(floatarr_batch)
        elif floatarr_batch.dtype == np.float64:
            self.append_dosages_batch_internal64(floatarr_batch)
        else:
            raise RuntimeError("Invalid append_dosages_batch() dosage array element type (float32 or float64 expected).")
        return


    cpdef close(self):
        cdef PglErr reterr = kPglRetSuccess
        if self._state_ptr:
            if SpgwGetVidx(self._state_ptr) != SpgwGetVariantCt(self._state_ptr):
                raise RuntimeError("PgenWriter.close() called when number of written variants (" + str(SpgwGetVidx(self._state_ptr)) + ") unequal to initially declared value (" + str(SpgwGetVariantCt(self._state_ptr)) + ").")
            reterr = SpgwFinish(self._state_ptr)
            if reterr != kPglRetSuccess:
                raise RuntimeError("PgenWriter.close(): SpgwFinish() error " + str(reterr))
            if self._nonref_flags:
                aligned_free(self._nonref_flags)
            PyMem_Free(self._state_ptr)
            self._state_ptr = NULL
        return


    cpdef __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return


    def __dealloc__(self):
        cdef PglErr reterr = kPglRetSuccess
        if self._state_ptr:
            if SpgwGetVidx(self._state_ptr) == SpgwGetVariantCt(self._state_ptr):
                reterr = SpgwFinish(self._state_ptr)
                if reterr != kPglRetSuccess:
                    raise RuntimeError("PgenWriter.__dealloc__(): SpgwFinish() error " + str(reterr))
            else:
                CleanupSpgw(self._state_ptr, &reterr)
            if self._nonref_flags:
                aligned_free(self._nonref_flags)
            PyMem_Free(self._state_ptr)
        return
