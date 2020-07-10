#include "pgenlib_ffi_support.h"
#include "include/pgenlib_read.h"
#include "pvar.h"  // includes Rcpp

class RPgenReader {
public:
  // imitates Python/pgenlib.pyx
  RPgenReader();

#if __cplusplus >= 201103L
  RPgenReader(const RPgenReader&) = delete;
  RPgenReader& operator=(const RPgenReader&) = delete;
#endif

  void Load(String filename, Nullable<List> pvar, Nullable<int> raw_sample_ct,
            Nullable<IntegerVector> sample_subset_1based);

  uint32_t GetRawSampleCt() const;

  uint32_t GetSubsetSize() const;

  uint32_t GetVariantCt() const;

  uint32_t GetAlleleCt(uint32_t variant_idx) const;

  uint32_t GetMaxAlleleCt() const;

  bool HardcallPhasePresent() const;

  void ReadIntHardcalls(IntegerVector buf, int variant_idx, int allele_idx);

  void ReadHardcalls(NumericVector buf, int variant_idx, int allele_idx);

  void Read(NumericVector buf, int variant_idx, int allele_idx);

  void ReadAlleles(IntegerMatrix acbuf,
                   Nullable<LogicalVector> phasepresent_buf, int variant_idx);

  void ReadAllelesNumeric(NumericMatrix acbuf,
                          Nullable<LogicalVector> phasepresent_buf,
                          int variant_idx);

  void ReadIntList(IntegerMatrix buf, IntegerVector variant_subset);

  void ReadList(NumericMatrix buf, IntegerVector variant_subset, bool meanimpute);

  void FillVariantScores(NumericVector result, NumericVector weights, Nullable<IntegerVector> variant_subset);

  void Close();

  ~RPgenReader();

private:
  plink2::PgenFileInfo* _info_ptr;
  plink2::RefcountedWptr* _allele_idx_offsetsp;
  plink2::RefcountedWptr* _nonref_flagsp;
  plink2::PgenReader* _state_ptr;
  uintptr_t* _subset_include_vec;
  uintptr_t* _subset_include_interleaved_vec;
  uint32_t* _subset_cumulative_popcounts;
  plink2::PgrSampleSubsetIndex _subset_index;
  uint32_t _subset_size;

  plink2::PgenVariant _pgv;

  plink2::VecW* _transpose_batch_buf;
  // kPglNypTransposeBatch (= 256) variants at a time, and then transpose
  uintptr_t* _multivar_vmaj_geno_buf;
  uintptr_t* _multivar_vmaj_phasepresent_buf;
  uintptr_t* _multivar_vmaj_phaseinfo_buf;
  uintptr_t* _multivar_smaj_geno_batch_buf;
  uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
  uintptr_t* _multivar_smaj_phasepresent_batch_buf;

  void SetSampleSubsetInternal(IntegerVector sample_subset_1based);

  void ReadAllelesPhasedInternal(int variant_idx);
};

RPgenReader::RPgenReader() : _info_ptr(nullptr),
                             _allele_idx_offsetsp(nullptr),
                             _nonref_flagsp(nullptr),
                             _state_ptr(nullptr) {
}

void RPgenReader::Load(String filename, Nullable<List> pvar,
                       Nullable<int> raw_sample_ct,
                       Nullable<IntegerVector> sample_subset_1based) {
  if (_info_ptr) {
    Close();
  }
  _info_ptr = static_cast<plink2::PgenFileInfo*>(malloc(sizeof(plink2::PgenFileInfo)));
  if (!_info_ptr) {
    stop("Out of memory");
  }
  plink2::PreinitPgfi(_info_ptr);
  uint32_t cur_sample_ct = UINT32_MAX;
  if (raw_sample_ct.isNotNull()) {
    cur_sample_ct = as<int>(raw_sample_ct.get());
  }
  uint32_t cur_variant_ct = UINT32_MAX;
  const char* fname = filename.get_cstring();
  plink2::PgenHeaderCtrl header_ctrl;
  uintptr_t pgfi_alloc_cacheline_ct;
  char errstr_buf[plink2::kPglErrstrBufBlen];
  if (PgfiInitPhase1(fname, cur_variant_ct, cur_sample_ct, 0, &header_ctrl, _info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
    stop(&(errstr_buf[7]));
  }
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  if (pvar.isNotNull()) {
    List pvarl = as<List>(pvar);
    if (strcmp_r_c(pvarl[0], "pvar")) {
      stop("pvar is not a pvar object");
    }
    XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvarl[1]);
    if (rp->GetVariantCt() != raw_variant_ct) {
      stop("pvar and pgen have different variant counts");
    }
    _allele_idx_offsetsp = rp->GetAlleleIdxOffsetsp();
    if (_allele_idx_offsetsp) {
      _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
    }
    _info_ptr->max_allele_ct = rp->GetMaxAlleleCt();
  } else {
    if (header_ctrl & 0x30) {
      // no need to zero-initialize this
      _allele_idx_offsetsp = plink2::CreateRefcountedWptr(raw_variant_ct + 1);
      _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
      // _info_ptr->max_allele_ct updated by PgfiInitPhase2() in this case
    }
    _info_ptr->max_allele_ct = 2;
  }
  if ((header_ctrl & 0xc0) == 0xc0) {
    // todo: load this in pvar, to enable consistency check.  we use a
    // (manually implemented) shared_ptr in preparation for this.
    const uintptr_t raw_variant_ctl = plink2::DivUp(raw_variant_ct, plink2::kBitsPerWord);
    // no need to zero-initialize this
    _nonref_flagsp = plink2::CreateRefcountedWptr(raw_variant_ctl + 1);
    _info_ptr->nonref_flags = _nonref_flagsp->p;
  }
  const uint32_t file_sample_ct = _info_ptr->raw_sample_ct;
  unsigned char* pgfi_alloc = nullptr;
  if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
    stop("Out of memory");
  }
  uint32_t max_vrec_width;
  uintptr_t pgr_alloc_cacheline_ct;
  if (PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, _info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
    if (pgfi_alloc && (!_info_ptr->vrtypes)) {
      plink2::aligned_free(pgfi_alloc);
    }
    stop(&(errstr_buf[7]));
  }
  if ((!_allele_idx_offsetsp) && (_info_ptr->gflags & 4)) {
    // Note that it's safe to be ignorant of multiallelic variants when
    // phase and dosage info aren't present; GetAlleleCt() then always returns
    // 2 when that isn't actually true, and all ALTs are treated as if they
    // were ALT1, but otherwise everything works properly.
    stop("Multiallelic variants and phase/dosage info simultaneously present; pvar required in this case");
  }
  _state_ptr = static_cast<plink2::PgenReader*>(malloc(sizeof(plink2::PgenReader)));
  if (!_state_ptr) {
    stop("Out of memory");
  }
  plink2::PreinitPgr(_state_ptr);
  plink2::PgrSetFreadBuf(nullptr, _state_ptr);
  const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
  const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
  const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
  const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
  const uintptr_t ac_byte_ct = plink2::RoundUpPow2(file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  uintptr_t multiallelic_hc_byte_ct = 0;
  if (_info_ptr->max_allele_ct != 2) {
    multiallelic_hc_byte_ct = 2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
  }
  const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
  unsigned char* pgr_alloc;
  if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
    stop("Out of memory");
  }
  plink2::PglErr reterr = PgrInit(fname, max_vrec_width, _info_ptr, _state_ptr, pgr_alloc);
  if (reterr != plink2::kPglRetSuccess) {
    if (!plink2::PgrGetFreadBuf(_state_ptr)) {
      plink2::aligned_free(pgr_alloc);
    }
    sprintf(errstr_buf, "PgrInit() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
  _subset_include_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _subset_include_interleaved_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);

#ifdef USE_AVX2
  _subset_include_interleaved_vec[-3] = 0;
  _subset_include_interleaved_vec[-2] = 0;
#endif
  _subset_include_interleaved_vec[-1] = 0;

  _subset_cumulative_popcounts = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
  _pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
  if (multiallelic_hc_byte_ct) {
    _pgv.patch_01_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.patch_01_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[ac_byte_ct]);
    _pgv.patch_10_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.patch_10_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[ac2_byte_ct]);
  } else {
    _pgv.patch_01_set = nullptr;
    _pgv.patch_01_vals = nullptr;
    _pgv.patch_10_set = nullptr;
    _pgv.patch_10_vals = nullptr;
  }
  _pgv.phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _pgv.dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);
  if (sample_subset_1based.isNotNull()) {
    SetSampleSubsetInternal(sample_subset_1based.get());
  } else {
    _subset_size = file_sample_ct;
  }
  _transpose_batch_buf = reinterpret_cast<plink2::VecW*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
  _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * genovec_byte_ct]);
  _multivar_vmaj_phasepresent_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
  _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
  _multivar_smaj_geno_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 4]);
  _multivar_smaj_phaseinfo_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
  _multivar_smaj_phasepresent_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  // pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
}

uint32_t RPgenReader::GetRawSampleCt() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return _info_ptr->raw_sample_ct;
}

uint32_t RPgenReader::GetSubsetSize() const {
  return _subset_size;
}

uint32_t RPgenReader::GetVariantCt() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return _info_ptr->raw_variant_ct;
}

uint32_t RPgenReader::GetAlleleCt(uint32_t variant_idx) const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (variant_idx >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (!_allele_idx_offsetsp) {
    return 2;
  }
  const uintptr_t* allele_idx_offsets = _allele_idx_offsetsp->p;
  return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
}

uint32_t RPgenReader::GetMaxAlleleCt() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return _info_ptr->max_allele_ct;
}

bool RPgenReader::HardcallPhasePresent() const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  return ((_info_ptr->gflags & plink2::kfPgenGlobalHardcallPhasePresent) != 0);
}

static const int32_t kGenoRInt32Quads[1024] ALIGNV16 = QUAD_TABLE256(0, 1, 2, NA_INTEGER);

void RPgenReader::ReadIntHardcalls(IntegerVector buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    char errstr_buf[256];
    sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", buf.size(), _subset_size);
    stop(errstr_buf);
  }
  plink2::PglErr reterr = PgrGet1(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoRInt32Quads, _subset_size, &buf[0]);
}

static const double kGenoRDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, NA_REAL);

void RPgenReader::ReadHardcalls(NumericVector buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    char errstr_buf[256];
    sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", buf.size(), _subset_size);
    stop(errstr_buf);
  }
  plink2::PglErr reterr = PgrGet1(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::GenoarrLookup16x8bx2(_pgv.genovec, kGenoRDoublePairs, _subset_size, &buf[0]);
}

void RPgenReader::Read(NumericVector buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    char errstr_buf[256];
    sprintf(errstr_buf, "buf has wrong length (%" PRIdPTR "; %u expected)", buf.size(), _subset_size);
    stop(errstr_buf);
  }
  uint32_t dosage_ct;
  plink2::PglErr reterr = PgrGet1D(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGet1D() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);
}

static const uint64_t kGenoToRIntcodeDPairs[32] ALIGNV16 = PAIR_TABLE16(0, 0x100000000LLU, 0x100000001LLU, 0x8000000080000000LLU);
static const int32_t kGenoToLogicalPhaseQuads[1024] ALIGNV16 = QUAD_TABLE256(1, 0, 1, NA_LOGICAL);

void RPgenReader::ReadAlleles(IntegerMatrix acbuf, Nullable<LogicalVector> phasepresent_buf, int variant_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if ((acbuf.nrow() != 2) || (acbuf.ncol() != static_cast<int>(_subset_size))) {
    char errstr_buf[256];
    sprintf(errstr_buf, "acbuf has wrong size (%dx%d; 2x%u expected)", acbuf.nrow(), acbuf.ncol(), _subset_size);
    stop(errstr_buf);
  }
  ReadAllelesPhasedInternal(variant_idx);
  plink2::GenoarrToAlleleCodes(kGenoToRIntcodeDPairs, _pgv.genovec, _subset_size, &acbuf[0]);
  const uintptr_t* allele_idx_offsets = _info_ptr->allele_idx_offsets;
  uint32_t cur_allele_ct = 2;
  if (allele_idx_offsets) {
    cur_allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
    if (cur_allele_ct != 2) {
      stop("multiallelic support under development");
    }
  }
  const uintptr_t* phasepresent = _pgv.phasepresent;
  const uintptr_t* phaseinfo = _pgv.phaseinfo;
  const uint32_t phasepresent_ct = _pgv.phasepresent_ct;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = phasepresent[0];
  if (!phasepresent_buf.isNotNull()) {
    if (cur_allele_ct == 2) {
      uint64_t* allele_codes_alias64 = R_CAST(uint64_t*, &acbuf[0]);
      for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
        const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
        if (plink2::IsSet(phaseinfo, sample_uidx)) {
          // 1|0
          allele_codes_alias64[sample_uidx] = 1;
        }
      }
    } else {
      int32_t* allele_codes = &acbuf[0];
      for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
        const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
        if (plink2::IsSet(phaseinfo, sample_uidx)) {
          const int32_t tmpval = allele_codes[2 * sample_uidx];
          allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
          allele_codes[2 * sample_uidx + 1] = tmpval;
        }
      }
    }
    return;
  }
  // Unfortunately, we can't use GenoarrPhasedToAlleleCodes directly, since
  // it's written for Python 1-byte bools instead of R 4-byte logical values.
  // (probable todo: allow the no-phasepresent_buf part to be called
  // separately)
  //
  // 0, 2 -> automatically phased.  3 -> NA_LOGICAL.
  // 1 -> assume unphased; then change to phased as necessary when iterating
  //      over phasepresent.
  int32_t* phasepresent_wbuf = &(as<LogicalVector>(phasepresent_buf)[0]);
  plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoToLogicalPhaseQuads, _subset_size, phasepresent_wbuf);
  if (cur_allele_ct == 2) {
    uint64_t* allele_codes_alias64 = R_CAST(uint64_t*, &acbuf[0]);
    for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
      const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
      phasepresent_wbuf[sample_uidx] = 1;
      if (plink2::IsSet(phaseinfo, sample_uidx)) {
        allele_codes_alias64[sample_uidx] = 1;
      }
    }
  } else {
    int32_t* allele_codes = &acbuf[0];
    for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
      const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
      phasepresent_wbuf[sample_uidx] = 1;
      if (plink2::IsSet(phaseinfo, sample_uidx)) {
        const int32_t tmpval = allele_codes[2 * sample_uidx];
        allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
        allele_codes[2 * sample_uidx + 1] = tmpval;
      }
    }
  }
}

static const double kGenoToRNumcodePairs[8] ALIGNV16 = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, NA_REAL, NA_REAL};

void RPgenReader::ReadAllelesNumeric(NumericMatrix acbuf, Nullable<LogicalVector> phasepresent_buf, int variant_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if ((acbuf.nrow() != 2) || (acbuf.ncol() != static_cast<int>(_subset_size))) {
    char errstr_buf[256];
    sprintf(errstr_buf, "acbuf has wrong size (%dx%d; 2x%u expected)", acbuf.nrow(), acbuf.ncol(), _subset_size);
    stop(errstr_buf);
  }
  ReadAllelesPhasedInternal(variant_idx);
  double* allele_codes = &acbuf[0];
  plink2::GenoarrLookup4x16b(_pgv.genovec, kGenoToRNumcodePairs, _subset_size, allele_codes);
  const uintptr_t* allele_idx_offsets = _info_ptr->allele_idx_offsets;
  uint32_t cur_allele_ct = 2;
  if (allele_idx_offsets) {
    cur_allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
    if (cur_allele_ct != 2) {
      stop("multiallelic support under development");
    }
  }
  const uintptr_t* phasepresent = _pgv.phasepresent;
  const uintptr_t* phaseinfo = _pgv.phaseinfo;
  const uint32_t phasepresent_ct = _pgv.phasepresent_ct;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = phasepresent[0];
  if (!phasepresent_buf.isNotNull()) {
    if (cur_allele_ct == 2) {
      for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
        const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
        if (plink2::IsSet(phaseinfo, sample_uidx)) {
          // 1|0
          allele_codes[2 * sample_uidx] = 1.0;
          allele_codes[2 * sample_uidx + 1] = 0.0;
        }
      }
    } else {
      for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
        const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
        if (plink2::IsSet(phaseinfo, sample_uidx)) {
          const double tmpval = allele_codes[2 * sample_uidx];
          allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
          allele_codes[2 * sample_uidx + 1] = tmpval;
        }
      }
    }
    return;
  }
  int32_t* phasepresent_wbuf = &(as<LogicalVector>(phasepresent_buf)[0]);
  plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoToLogicalPhaseQuads, _subset_size, phasepresent_wbuf);
  if (cur_allele_ct == 2) {
    for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
      const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
      phasepresent_wbuf[sample_uidx] = 1;
      if (plink2::IsSet(phaseinfo, sample_uidx)) {
        allele_codes[2 * sample_uidx] = 1.0;
        allele_codes[2 * sample_uidx + 1] = 0.0;
      }
    }
  } else {
    for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
      const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
      phasepresent_wbuf[sample_uidx] = 1;
      if (plink2::IsSet(phaseinfo, sample_uidx)) {
        const double tmpval = allele_codes[2 * sample_uidx];
        allele_codes[2 * sample_uidx] = allele_codes[2 * sample_uidx + 1];
        allele_codes[2 * sample_uidx + 1] = tmpval;
      }
    }
  }
}

void RPgenReader::ReadIntList(IntegerMatrix buf, IntegerVector variant_subset) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  // assume that buf has the correct dimensions
  const uintptr_t vsubset_size = variant_subset.size();
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  int32_t* buf_iter = &buf[0];
  for (uintptr_t col_idx = 0; col_idx != vsubset_size; ++col_idx) {
    uint32_t variant_idx = variant_subset[col_idx] - 1;
    if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
      char errstr_buf[256];
      sprintf(errstr_buf, "variant_subset element out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
      stop(errstr_buf);
    }
    plink2::PglErr reterr = PgrGet(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, _pgv.genovec);
    if (reterr != plink2::kPglRetSuccess) {
      char errstr_buf[256];
      sprintf(errstr_buf, "PgrGet() error %d", static_cast<int>(reterr));
      stop(errstr_buf);
    }
    plink2::GenoarrLookup256x4bx4(_pgv.genovec, kGenoRInt32Quads, _subset_size, buf_iter);
    buf_iter = &(buf_iter[_subset_size]);
  }
}

void RPgenReader::ReadList(NumericMatrix buf, IntegerVector variant_subset, bool meanimpute) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  // assume that buf has the correct dimensions
  const uintptr_t vsubset_size = variant_subset.size();
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  double* buf_iter = &buf[0];
  for (uintptr_t col_idx = 0; col_idx != vsubset_size; ++col_idx) {
    uint32_t variant_idx = variant_subset[col_idx] - 1;
    if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
      char errstr_buf[256];
      sprintf(errstr_buf, "variant_subset element out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
      stop(errstr_buf);
    }
    uint32_t dosage_ct;
    plink2::PglErr reterr = PgrGetD(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
    if (reterr != plink2::kPglRetSuccess) {
      char errstr_buf[256];
      sprintf(errstr_buf, "PgrGetD() error %d", static_cast<int>(reterr));
      stop(errstr_buf);
    }
    if (!meanimpute) {
      plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter);
    } else {
      plink2::ZeroTrailingNyps(_subset_size, _pgv.genovec);
      if (plink2::Dosage16ToDoublesMeanimpute(_pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter)) {
        char errstr_buf[256];
        sprintf(errstr_buf, "variant %d has only missing dosages", variant_idx + 1);
        stop(errstr_buf);
      }
    }
    buf_iter = &(buf_iter[_subset_size]);
  }
}

void RPgenReader::FillVariantScores(NumericVector result, NumericVector weights, Nullable<IntegerVector> variant_subset) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (weights.size() != _subset_size) {
    char errstr_buf[256];
    sprintf(errstr_buf, "weights.size()=%td doesn't match pgen sample-subset size=%d", weights.size(), _subset_size);
    stop(errstr_buf);
  }
  const int raw_variant_ct = _info_ptr->raw_variant_ct;
  const int* variant_idx_ints = nullptr;
  uintptr_t variant_ct = raw_variant_ct;
  if (variant_subset.isNotNull()) {
    IntegerVector vs = as<IntegerVector>(variant_subset);
    variant_idx_ints = &(vs[0]);
    variant_ct = vs.size();
  }
  for (uintptr_t ulii = 0; ulii != variant_ct; ++ulii) {
    int variant_idx = ulii;
    if (variant_idx_ints) {
      variant_idx = variant_idx_ints[ulii] - 1;
      if ((variant_idx < 0) || (variant_idx >= raw_variant_ct)) {
        char errstr_buf[256];
        sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
        stop(errstr_buf);
      }
    }
    uint32_t dosage_ct;
    plink2::PglErr reterr = plink2::PgrGetD(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
    if (reterr != plink2::kPglRetSuccess) {
      char errstr_buf[256];
      sprintf(errstr_buf, "PgrGetD() error %d", static_cast<int>(reterr));
      stop(errstr_buf);
    }
    plink2::ZeroTrailingNyps(_subset_size, _pgv.genovec);
    const double* wts = &(weights[0]);
    result[ulii] = plink2::LinearCombinationMeanimpute(wts, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct);
  }
}

void RPgenReader::Close() {
  // don't bother propagating file close errors for now
  if (_info_ptr) {
    CondReleaseRefcountedWptr(&_allele_idx_offsetsp);
    CondReleaseRefcountedWptr(&_nonref_flagsp);
    if (_info_ptr->vrtypes) {
      plink2::aligned_free(_info_ptr->vrtypes);
    }
    plink2::PglErr reterr = plink2::kPglRetSuccess;
    plink2::CleanupPgfi(_info_ptr, &reterr);
    free(_info_ptr);
    _info_ptr = nullptr;
  }
  if (_state_ptr) {
    plink2::PglErr reterr = plink2::kPglRetSuccess;
    plink2::CleanupPgr(_state_ptr, &reterr);
    if (PgrGetFreadBuf(_state_ptr)) {
      plink2::aligned_free(PgrGetFreadBuf(_state_ptr));
    }
    free(_state_ptr);
    _state_ptr = nullptr;
  }
  _subset_size = 0;
}

void RPgenReader::SetSampleSubsetInternal(IntegerVector sample_subset_1based) {
  const uint32_t raw_sample_ct = _info_ptr->raw_sample_ct;
  const uint32_t raw_sample_ctv = plink2::DivUp(raw_sample_ct, plink2::kBitsPerVec);
  const uint32_t raw_sample_ctaw = raw_sample_ctv * plink2::kWordsPerVec;
  uintptr_t* sample_include = _subset_include_vec;
  plink2::ZeroWArr(raw_sample_ctaw, sample_include);
  const uint32_t subset_size = sample_subset_1based.size();
  if (subset_size == 0) {
    stop("Empty sample_subset is not currently permitted");
  }
  uint32_t sample_uidx = sample_subset_1based[0] - 1;
  uint32_t idx = 0;
  uint32_t next_uidx;
  while (1) {
    if (sample_uidx >= raw_sample_ct) {
      char errstr_buf[256];
      sprintf(errstr_buf, "sample number out of range (%d; must be 1..%u)", static_cast<int>(sample_uidx + 1), raw_sample_ct);
      stop(errstr_buf);
    }
    plink2::SetBit(sample_uidx, sample_include);
    if (++idx == subset_size) {
      break;
    }
    next_uidx = sample_subset_1based[idx] - 1;

    // prohibit this since it implies that the caller expects genotypes to be
    // returned in a different order
    if (next_uidx <= sample_uidx) {
      stop("sample_subset is not in strictly increasing order");
    }
    sample_uidx = next_uidx;
  }
  plink2::FillInterleavedMaskVec(sample_include, raw_sample_ctv, _subset_include_interleaved_vec);
  const uint32_t raw_sample_ctl = plink2::DivUp(raw_sample_ct, plink2::kBitsPerWord);
  plink2::FillCumulativePopcounts(sample_include, raw_sample_ctl, _subset_cumulative_popcounts);
  plink2::PgrSetSampleSubsetIndex(_subset_cumulative_popcounts, _state_ptr, &_subset_index);
  _subset_size = subset_size;
}

void RPgenReader::ReadAllelesPhasedInternal(int variant_idx) {
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  plink2::PglErr reterr = plink2::PgrGetMP(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, &_pgv);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "PgrGetMP() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

RPgenReader::~RPgenReader() {
  Close();
}

// [[Rcpp::export]]
SEXP NewPgen(String filename, Nullable<List> pvar = R_NilValue,
             Nullable<int> raw_sample_ct = R_NilValue,
             Nullable<IntegerVector> sample_subset = R_NilValue) {
  XPtr<class RPgenReader> pgen(new RPgenReader(), true);
  pgen->Load(filename, pvar, raw_sample_ct, sample_subset);
  return List::create(_["class"] = "pgen", _["pgen"] = pgen);
}

// [[Rcpp::export]]
int GetRawSampleCt(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return rp->GetRawSampleCt();
}

// [[Rcpp::export]]
int GetVariantCt(List pvar_or_pgen) {
  const char* c_str = as<String>(pvar_or_pgen[0]).get_cstring();
  if (!strcmp(c_str, "pvar")) {
    XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar_or_pgen[1]);
    return rp->GetVariantCt();
  } else if (!strcmp(c_str, "pgen")) {
    XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pvar_or_pgen[1]);
    return rp->GetVariantCt();
  }
  stop("pvar_or_pgen is not a pvar or pgen object");
}

// [[Rcpp::export]]
int GetAlleleCt(List pvar_or_pgen, int variant_num) {
  const char* c_str = as<String>(pvar_or_pgen[0]).get_cstring();
  const uint32_t variant_idx = variant_num - 1;
  if (!strcmp(c_str, "pvar")) {
    XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar_or_pgen[1]);
    return rp->GetAlleleCt(variant_idx);
  } else if (!strcmp(c_str, "pgen")) {
    XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pvar_or_pgen[1]);
    return rp->GetAlleleCt(variant_idx);
  }
  stop("pvar_or_pgen is not a pvar or pgen object");
}

// [[Rcpp::export]]
int GetMaxAlleleCt(List pvar_or_pgen) {
  const char* c_str = as<String>(pvar_or_pgen[0]).get_cstring();
  if (!strcmp(c_str, "pvar")) {
    XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar_or_pgen[1]);
    return rp->GetMaxAlleleCt();
  } else if (!strcmp(c_str, "pgen")) {
    XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pvar_or_pgen[1]);
    return rp->GetMaxAlleleCt();
  }
  stop("pvar_or_pgen is not a pvar or pgen object");
}

// [[Rcpp::export]]
bool HardcallPhasePresent(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return rp->HardcallPhasePresent();
}

// [[Rcpp::export]]
NumericVector Buf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return NumericVector(rp->GetSubsetSize());
}

// [[Rcpp::export]]
NumericVector AlleleCodeBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return NumericMatrix(2, rp->GetSubsetSize());
}

// [[Rcpp::export]]
IntegerVector IntBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return IntegerVector(rp->GetSubsetSize());
}

// [[Rcpp::export]]
IntegerVector IntAlleleCodeBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return IntegerMatrix(2, rp->GetSubsetSize());
}

// [[Rcpp::export]]
LogicalVector BoolBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return LogicalVector(rp->GetSubsetSize());
}

// _num indicates R 1-based indexing.

// [[Rcpp::export]]
void ReadHardcalls(List pgen, SEXP buf, int variant_num, int allele_num = 2) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  if (Rf_isMatrix(buf)) {
    // otherwise the original buffer is not modified by Read[Int]Hardcalls
    stop("buf must be a non-matrix vector");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  const int allele_idx = allele_num - 1;
  if (TYPEOF(buf) == REALSXP) {
    rp->ReadHardcalls(buf, variant_idx, allele_idx);
  } else if (TYPEOF(buf) == INTSXP) {
    rp->ReadIntHardcalls(buf, variant_idx, allele_idx);
  } else {
    stop("Unsupported buf type");
  }
}

// [[Rcpp::export]]
void Read(List pgen, NumericVector buf, int variant_num, int allele_num = 2) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  if (Rf_isMatrix(buf)) {
    stop("buf must be a non-matrix vector");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  rp->Read(buf, variant_num - 1, allele_num - 1);
}

// [[Rcpp::export]]
void ReadAlleles(List pgen, SEXP acbuf, int variant_num, Nullable<LogicalVector> phasepresent_buf = R_NilValue) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  // in this case, integer may be a more appropriate default than numeric?
  if (TYPEOF(acbuf) == INTSXP) {
    rp->ReadAlleles(acbuf, phasepresent_buf, variant_idx);
  } else if (TYPEOF(acbuf) == REALSXP) {
    rp->ReadAllelesNumeric(acbuf, phasepresent_buf, variant_idx);
  } else {
    stop("Unsupported acbuf type");
  }
}

// [[Rcpp::export]]
IntegerMatrix ReadIntList(List pgen, IntegerVector variant_subset) {
  // return value: rows = samples, columns = variants (from R's perspective)
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  IntegerMatrix result(rp->GetSubsetSize(), variant_subset.size());
  rp->ReadIntList(result, variant_subset);
  return result;
}

// [[Rcpp::export]]
NumericMatrix ReadList(List pgen, IntegerVector variant_subset, bool meanimpute = false) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  NumericMatrix result(rp->GetSubsetSize(), variant_subset.size());
  rp->ReadList(result, variant_subset, meanimpute);
  return result;
}

// [[Rcpp::export]]
NumericVector VariantScores(List pgen, NumericVector weights,
                            Nullable<IntegerVector> variant_subset = R_NilValue) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  int variant_ct;
  if (variant_subset.isNotNull()) {
    variant_ct = as<IntegerVector>(variant_subset).size();
  } else {
    variant_ct = rp->GetVariantCt();
  }
  NumericVector result(variant_ct);
  rp->FillVariantScores(result, weights, variant_subset);
  return result;
}

// [[Rcpp::export]]
void ClosePgen(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  rp->Close();
}
