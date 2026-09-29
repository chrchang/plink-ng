#include "include/pgenlib_ffi_support.h"
#include "include/pgenlib_read.h"
#include "include/pgenlib_write.h"
#include "pvar.h"  // includes Rcpp

class RPgenReader {
public:
  // similar to Python/pgenlib.pyx ; has a bit more functionality as of Feb
  // 2025
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

  uint32_t GetVrtype(uint32_t variant_idx) const;

  bool HardcallPhasePresent() const;

  void ReadIntHardcalls(IntegerVector buf, int variant_idx, int allele_idx);

  void ReadHardcalls(NumericVector buf, int variant_idx, int allele_idx);

  void ReadIntMaybeSparseHardcalls(IntegerVector buf, int variant_idx, int allele_idx, int max_difflist_len, IntegerVector* sample_nums_ptr, IntegerVector* allele_counts_ptr);

  void ReadMaybeSparseHardcalls(NumericVector buf, int variant_idx, int allele_idx, int max_difflist_len, IntegerVector* sample_nums_ptr, NumericVector* allele_dosages_ptr);

  void Read(NumericVector buf, int variant_idx, int allele_idx);

  void ReadMaybeSparse(NumericVector buf, int variant_idx, int allele_idx, int max_difflist_len, IntegerVector* sample_nums_ptr, NumericVector* allele_dosages_ptr);

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

  uintptr_t* _raregeno_buf;
  uint32_t* _difflist_sample_ids_buf;

  plink2::VecW* _transpose_batch_buf;
  // kPglNypTransposeBatch (= 256) variants at a time, and then transpose
  uintptr_t* _multivar_vmaj_geno_buf;
  uintptr_t* _multivar_vmaj_phasepresent_buf;
  uintptr_t* _multivar_vmaj_phaseinfo_buf;
  uintptr_t* _multivar_smaj_geno_batch_buf;
  uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
  uintptr_t* _multivar_smaj_phasepresent_batch_buf;

  void SetSampleSubsetInternal(IntegerVector sample_subset_1based);

  void ReadMaybeSparseHardcallsInternal(int variant_idx, int max_simple_difflist_len, uint32_t* difflist_common_geno_ptr, uint32_t* difflist_len_ptr);

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
  if (plink2::PgfiInitPhase1(fname, nullptr, cur_variant_ct, cur_sample_ct, &header_ctrl, _info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
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
  if (plink2::PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, _info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
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
  const uint32_t is_not_plink1_bed = (_info_ptr->vrtypes != nullptr);
  uintptr_t raregeno_byte_ct = 0;
  uintptr_t difflist_sample_ids_byte_ct = 0;
  if (is_not_plink1_bed) {
    const uint32_t max_stored_single_difflist_len = file_sample_ct / plink2::kPglMaxDifflistLenDivisor;
    raregeno_byte_ct = plink2::DivUp(2 * max_stored_single_difflist_len, plink2::kNypsPerVec) * plink2::kBytesPerVec;
    // ReadSparse() uses max_difflist_len = 3 * max_stored_single_difflist_len.
    difflist_sample_ids_byte_ct = plink2::RoundUpPow2(3 * max_stored_single_difflist_len * sizeof(int32_t), plink2::kBytesPerVec);
  }
  const uintptr_t ac_byte_ct = plink2::RoundUpPow2(file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
  uintptr_t multiallelic_hc_byte_ct = 0;
  if (_info_ptr->max_allele_ct != 2) {
    multiallelic_hc_byte_ct = 2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
  }
  const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
  unsigned char* pgr_alloc;
  if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + raregeno_byte_ct + difflist_sample_ids_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
    stop("Out of memory");
  }
  plink2::PglErr reterr = plink2::PgrInit(fname, max_vrec_width, _info_ptr, _state_ptr, pgr_alloc);
  if (reterr != plink2::kPglRetSuccess) {
    if (!plink2::PgrGetFreadBuf(_state_ptr)) {
      plink2::aligned_free(pgr_alloc);
    }
    snprintf(errstr_buf, plink2::kPglErrstrBufBlen, "PgrInit() error %d", static_cast<int>(reterr));
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
  if (is_not_plink1_bed) {
    _raregeno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[raregeno_byte_ct]);
    _difflist_sample_ids_buf = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[difflist_sample_ids_byte_ct]);
  } else {
    _raregeno_buf = nullptr;
    _difflist_sample_ids_buf = nullptr;
  }
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
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
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

uint32_t RPgenReader::GetVrtype(uint32_t variant_idx) const {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (variant_idx >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  return plink2::PgrGetVrtype(_state_ptr, variant_idx);
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
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    using namespace plink2;
    char errstr_buf[256];
    char* write_iter = strcpya_k(errstr_buf, "buf has wrong length (");
    write_iter = wtoa(buf.size(), write_iter);
    write_iter = strcpya_k(write_iter, "; ");
    write_iter = u32toa(_subset_size, write_iter);
    strcpy_k(write_iter, " expected)");
    stop(errstr_buf);
  }
  plink2::PglErr reterr = plink2::PgrGet1(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "PgrGet1() error %d", static_cast<int>(reterr));
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
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    using namespace plink2;
    char errstr_buf[256];
    char* write_iter = strcpya_k(errstr_buf, "buf has wrong length (");
    write_iter = wtoa(buf.size(), write_iter);
    write_iter = strcpya_k(write_iter, "; ");
    write_iter = u32toa(_subset_size, write_iter);
    strcpy_k(write_iter, " expected)");
    stop(errstr_buf);
  }
  plink2::PglErr reterr = plink2::PgrGet1(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "PgrGet1() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::GenoarrLookup16x8bx2(_pgv.genovec, kGenoRDoublePairs, _subset_size, &buf[0]);
}

static const int32_t kGenoRInt32QuadsFlipped[1024] ALIGNV16 = QUAD_TABLE256(2, 1, 0, NA_INTEGER);

void RPgenReader::ReadIntMaybeSparseHardcalls(IntegerVector buf, int variant_idx, int allele_idx, int max_difflist_len, IntegerVector* sample_nums_ptr, IntegerVector* allele_counts_ptr) {
  uint32_t difflist_common_geno;
  uint32_t difflist_len;
  ReadMaybeSparseHardcallsInternal(variant_idx, max_difflist_len, &difflist_common_geno, &difflist_len);
  const int32_t* quad_table = (allele_idx == 0)? kGenoRInt32QuadsFlipped : kGenoRInt32Quads;
  if (((allele_idx == 0) && (difflist_common_geno != 2)) ||
      ((allele_idx == 1) && (difflist_common_geno != 0)) ||
      (static_cast<uint32_t>(allele_idx) > 1)) {
    if (buf.size() != _subset_size) {
      // Note that buf is *not* required to be the expected size when we return
      // the sparse representation.
      using namespace plink2;
      char errstr_buf[256];
      char* write_iter = strcpya_k(errstr_buf, "buf has wrong length (");
      write_iter = wtoa(buf.size(), write_iter);
      write_iter = strcpya_k(write_iter, "; ");
      write_iter = u32toa(_subset_size, write_iter);
      strcpy_k(write_iter, " expected)");
      stop(errstr_buf);
    }
    if (difflist_common_geno != UINT32_MAX) {
      // Sparse, but not w.r.t. the correct allele.  Just return dense
      // representation.
      plink2::PgrDifflistToGenovecUnsafe(_raregeno_buf, _difflist_sample_ids_buf, difflist_common_geno, _subset_size, difflist_len, _pgv.genovec);
    }
    plink2::GenoarrLookup256x4bx4(_pgv.genovec, quad_table, _subset_size, &buf[0]);
    return;
  }
  *sample_nums_ptr = IntegerVector(difflist_len);
  int* sample_nums_data = &((*sample_nums_ptr)[0]);
  for (uint32_t uii = 0; uii != difflist_len; ++uii) {
    sample_nums_data[uii] = _difflist_sample_ids_buf[uii] + 1;
  }
  *allele_counts_ptr = IntegerVector(difflist_len);
  int* allele_counts_data = &((*allele_counts_ptr)[0]);
  plink2::GenoarrLookup256x4bx4(_raregeno_buf, quad_table, difflist_len, allele_counts_data);
}

static const double kGenoRDoublePairsFlipped[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, NA_REAL);

void RPgenReader::ReadMaybeSparseHardcalls(NumericVector buf, int variant_idx, int allele_idx, int max_difflist_len, IntegerVector* sample_nums_ptr, NumericVector* allele_dosages_ptr) {
  uint32_t difflist_common_geno;
  uint32_t difflist_len;
  ReadMaybeSparseHardcallsInternal(variant_idx, max_difflist_len, &difflist_common_geno, &difflist_len);
  const double* pair_table = (allele_idx == 0)? kGenoRDoublePairsFlipped : kGenoRDoublePairs;
  if (((allele_idx == 0) && (difflist_common_geno != 2)) ||
      ((allele_idx == 1) && (difflist_common_geno != 0)) ||
      (static_cast<uint32_t>(allele_idx) > 1)) {
    if (buf.size() != _subset_size) {
      using namespace plink2;
      char errstr_buf[256];
      char* write_iter = strcpya_k(errstr_buf, "buf has wrong length (");
      write_iter = wtoa(buf.size(), write_iter);
      write_iter = strcpya_k(write_iter, "; ");
      write_iter = u32toa(_subset_size, write_iter);
      strcpy_k(write_iter, " expected)");
      stop(errstr_buf);
    }
    if (difflist_common_geno != UINT32_MAX) {
      // Sparse, but not w.r.t. the correct allele.  Just return dense
      // representation.
      plink2::PgrDifflistToGenovecUnsafe(_raregeno_buf, _difflist_sample_ids_buf, difflist_common_geno, _subset_size, difflist_len, _pgv.genovec);
    }
    plink2::GenoarrLookup16x8bx2(_pgv.genovec, pair_table, _subset_size, &buf[0]);
    return;
  }
  *sample_nums_ptr = IntegerVector(difflist_len);
  int* sample_nums_data = &((*sample_nums_ptr)[0]);
  for (uint32_t uii = 0; uii != difflist_len; ++uii) {
    sample_nums_data[uii] = _difflist_sample_ids_buf[uii] + 1;
  }
  *allele_dosages_ptr = NumericVector(difflist_len);
  double* allele_dosages_data = &((*allele_dosages_ptr)[0]);
  plink2::GenoarrLookup16x8bx2(_raregeno_buf, pair_table, difflist_len, allele_dosages_data);
}

void RPgenReader::Read(NumericVector buf, int variant_idx, int allele_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  if (buf.size() != _subset_size) {
    using namespace plink2;
    char errstr_buf[256];
    char* write_iter = strcpya_k(errstr_buf, "buf has wrong length (");
    write_iter = wtoa(buf.size(), write_iter);
    write_iter = strcpya_k(write_iter, "; ");
    write_iter = u32toa(_subset_size, write_iter);
    strcpy_k(write_iter, " expected)");
    stop(errstr_buf);
  }
  uint32_t dosage_ct;
  plink2::PglErr reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "PgrGet1D() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);
}

void RPgenReader::ReadMaybeSparse(NumericVector buf, int variant_idx, int allele_idx, int max_difflist_len, IntegerVector* sample_nums_ptr, NumericVector* allele_dosages_ptr) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  uint32_t dosage_ct;
  uint16_t difflist_common_dosage;
  plink2::PglErr reterr = plink2::PgrGetDMaybeSparse(_subset_include_vec, _subset_index, _subset_size, variant_idx, max_difflist_len, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct, &difflist_common_dosage, _difflist_sample_ids_buf);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    char* log = plink2::PglReturnLog();
    if (log == nullptr) {
      snprintf(errstr_buf, 256, "PgrGetDMaybeSparse() error %d", static_cast<int>(reterr));
    } else {
      snprintf(errstr_buf, 256, "PgrGetDMaybeSparse() error %d: %s", static_cast<int>(reterr), log);
    }
    stop(errstr_buf);
  }
  if (difflist_common_dosage == 1) {
    if (buf.size() != _subset_size) {
      using namespace plink2;
      char errstr_buf[256];
      char* write_iter = strcpya_k(errstr_buf, "buf has wrong length (");
      write_iter = wtoa(buf.size(), write_iter);
      write_iter = strcpya_k(write_iter, "; ");
      write_iter = u32toa(_subset_size, write_iter);
      strcpy_k(write_iter, " expected)");
      stop(errstr_buf);
    }
    const double* pair_table = (allele_idx == 0)? kGenoRDoublePairsFlipped : kGenoRDoublePairs;
    plink2::Dosage16ToDoubles(pair_table, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);
    return;
  }
  *sample_nums_ptr = IntegerVector(dosage_ct);
  int* sample_nums_data = &((*sample_nums_ptr)[0]);
  for (uint32_t uii = 0; uii != dosage_ct; ++uii) {
    sample_nums_data[uii] = _difflist_sample_ids_buf[uii] + 1;
  }
  uint16_t* dosage_main = _pgv.dosage_main;
  *allele_dosages_ptr = NumericVector(dosage_ct);
  double* allele_dosages_data = &((*allele_dosages_ptr)[0]);
  if (allele_idx == 0) {
    plink2::BiallelicDosage16Invert(dosage_ct, dosage_main);
  }
  plink2::DenseDosage16ToDoubles(dosage_main, dosage_ct, NA_REAL, allele_dosages_data);
}

static const uint64_t kGenoToRIntcodeDPairs[32] ALIGNV16 = PAIR_TABLE16(0, 0x100000000LLU, 0x100000001LLU, 0x8000000080000000LLU);
static const int32_t kGenoToLogicalPhaseQuads[1024] ALIGNV16 = QUAD_TABLE256(1, 0, 1, NA_LOGICAL);

void RPgenReader::ReadAlleles(IntegerMatrix acbuf, Nullable<LogicalVector> phasepresent_buf, int variant_idx) {
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  if ((acbuf.nrow() != 2) || (acbuf.ncol() != static_cast<int>(_subset_size))) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "acbuf has wrong size (%dx%d; 2x%u expected)", acbuf.nrow(), acbuf.ncol(), _subset_size);
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
    int32_t* allele_codes = &acbuf[0];
    if (cur_allele_ct == 2) {
      const uint64_t one = 1;
      for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
        const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
        if (plink2::IsSet(phaseinfo, sample_uidx)) {
          // 1|0
          memcpy(&(allele_codes[2 * sample_uidx]), &one, sizeof(int64_t));
        }
      }
    } else {
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
  int32_t* allele_codes = &acbuf[0];
  if (cur_allele_ct == 2) {
    const uint64_t one = 1;
    for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
      const uintptr_t sample_uidx = plink2::BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
      phasepresent_wbuf[sample_uidx] = 1;
      if (plink2::IsSet(phaseinfo, sample_uidx)) {
        memcpy(&(allele_codes[2 * sample_uidx]), &one, sizeof(int64_t));
      }
    }
  } else {
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
    snprintf(errstr_buf, 256, "acbuf has wrong size (%dx%d; 2x%u expected)", acbuf.nrow(), acbuf.ncol(), _subset_size);
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
      snprintf(errstr_buf, 256, "variant_subset element out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
      stop(errstr_buf);
    }
    plink2::PglErr reterr = plink2::PgrGet(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, _pgv.genovec);
    if (reterr != plink2::kPglRetSuccess) {
      char errstr_buf[256];
      snprintf(errstr_buf, 256, "PgrGet() error %d", static_cast<int>(reterr));
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
      snprintf(errstr_buf, 256, "variant_subset element out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
      stop(errstr_buf);
    }
    uint32_t dosage_ct;
    plink2::PglErr reterr = plink2::PgrGetD(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
    if (reterr != plink2::kPglRetSuccess) {
      char errstr_buf[256];
      snprintf(errstr_buf, 256, "PgrGetD() error %d", static_cast<int>(reterr));
      stop(errstr_buf);
    }
    if (!meanimpute) {
      plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter);
    } else {
      plink2::ZeroTrailingNyps(_subset_size, _pgv.genovec);
      if (plink2::Dosage16ToDoublesMeanimpute(_pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, buf_iter)) {
        char errstr_buf[256];
        snprintf(errstr_buf, 256, "variant %d has only missing dosages", variant_idx + 1);
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
    using namespace plink2;
    char errstr_buf[256];
    char* write_iter = strcpya_k(errstr_buf, "weights.size()=");
    write_iter = wtoa(weights.size(), write_iter);
    write_iter = strcpya_k(write_iter, " doesn't match pgen sample-subset size=");
    write_iter = wtoa(_subset_size, write_iter);
    *write_iter = '\0';
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
        snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, raw_variant_ct);
        stop(errstr_buf);
      }
    }
    uint32_t dosage_ct;
    plink2::PglErr reterr = plink2::PgrGetD(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
    if (reterr != plink2::kPglRetSuccess) {
      char errstr_buf[256];
      snprintf(errstr_buf, 256, "PgrGetD() error %d", static_cast<int>(reterr));
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
      snprintf(errstr_buf, 256, "sample number out of range (%d; must be 1..%u)", static_cast<int>(sample_uidx + 1), raw_sample_ct);
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

void RPgenReader::ReadMaybeSparseHardcallsInternal(int variant_idx, int max_difflist_len, uint32_t* difflist_common_geno_ptr, uint32_t* difflist_len_ptr) {
  // Fills {_raregeno_buf, _difflist_sample_ids_buf, *difflist_common_geno_ptr,
  // *difflist_len_ptr} iff hardcalls are either (i) stored as a simple
  // difflist no longer than max_difflist_len, or (ii) stored as a list of
  // differences from an earlier variant, which is itself stored as a simple
  // difflist, and the sum of the two list lengths is <= max_difflist_len.
  //
  // Otherwise, fills _pgv.genovec and sets *difflist_common_geno_ptr to
  // UINT32_MAX.
  if (!_info_ptr) {
    stop("pgen is closed");
  }
  const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
  if (static_cast<uint32_t>(variant_idx) >= raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  plink2::PglErr reterr = plink2::PgrGetDifflistOrGenovec(_subset_include_vec, _subset_index, _subset_size, max_difflist_len, variant_idx, _state_ptr, _pgv.genovec, difflist_common_geno_ptr, _raregeno_buf, _difflist_sample_ids_buf, difflist_len_ptr);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "PgrGetDifflistOrGenovec() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

void RPgenReader::ReadAllelesPhasedInternal(int variant_idx) {
  if (static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _info_ptr->raw_variant_ct);
    stop(errstr_buf);
  }
  plink2::PglErr reterr = plink2::PgrGetMP(_subset_include_vec, _subset_index, _subset_size, variant_idx, _state_ptr, &_pgv);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "PgrGetMP() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

RPgenReader::~RPgenReader() {
  Close();
}

//' Opens a .pgen or PLINK 1 .bed file.
//'
//' @param filename .pgen/.bed file path.
//' @param pvar Object (see NewPvar()) corresponding to the .pgen's companion
//' .pvar; technically optional, but necessary for some functionality.  In
//' particular, at multiallelic variants, all ALT alleles may be collapsed
//' together when .pvar information is not available.
//' @param raw_sample_ct Number of samples in file; required if it's a PLINK 1
//' .bed file, otherwise optional.
//' @param sample_subset List of 1-based positions of samples to load;
//' optional, all samples are loaded if this is not specified.
//' @return A pgen object, which can be queried for genotype/dosage data.
//' @export
// [[Rcpp::export]]
SEXP NewPgen(String filename, Nullable<List> pvar = R_NilValue,
             Nullable<int> raw_sample_ct = R_NilValue,
             Nullable<IntegerVector> sample_subset = R_NilValue) {
  XPtr<class RPgenReader> pgen(new RPgenReader(), true);
  pgen->Load(filename, pvar, raw_sample_ct, sample_subset);
  return List::create(_["class"] = "pgen", _["pgen"] = pgen);
}

//' Returns the number of samples in the file.
//'
//' @param pgen Object returned by NewPgen().
//' @return Number of samples.
//' @export
// [[Rcpp::export]]
int GetRawSampleCt(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return rp->GetRawSampleCt();
}

//' Returns the number of variants in the file.
//'
//' @param pvar_or_pgen Object returned by NewPvar() or NewPgen().
//' @return Number of variants.
//' @export
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

// Throughout this codebase, _num indicates R 1-based indexing.

//' Returns the effective number of alleles for a variant.  Note that if no
//' pvar was provided to the NewPgen() call, this function may return 2 even at
//' multiallelic variants, since the .pgen may not store allele-count
//' information.
//'
//' @param pvar_or_pgen Object returned by NewPvar() or NewPgen().
//' @param variant_num Variant index (1-based).
//' @return max(2, <number of alleles the variant_numth variant is known to
//' have>).  Note that if no
//' @export
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

//' Returns the maximum GetAlleleCt() value across all variants in the file.
//'
//' @param pvar_or_pgen Object returned by NewPvar() or NewPgen().
//' @return Maximum GetAlleleCt() value across all variants.
//' @export
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

//' Returns whether explicitly phased hardcalls are present.
//'
//' @param pgen Object returned by NewPgen().
//' @return TRUE if the file contains at least one phased heterozygous
//' hardcall, FALSE otherwise.
//' @export
// [[Rcpp::export]]
bool HardcallPhasePresent(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return rp->HardcallPhasePresent();
}

//' Returns a numeric buffer that Read() or ReadHardcalls() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Numeric vector with appropriate length for Read() and
//' ReadHardcalls().
//' @export
// [[Rcpp::export]]
NumericVector Buf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return NumericVector(rp->GetSubsetSize());
}

//' Returns an empty two-row numeric matrix that ReadAlleles() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Numeric matrix with two rows, and appropriate number of columns for
//' ReadAlleles().
//' @export
// [[Rcpp::export]]
NumericVector AlleleCodeBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return NumericMatrix(2, rp->GetSubsetSize());
}

//' Returns an integer buffer that ReadHardcalls() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Integer vector with appropriate length for ReadHardcalls().
//' @export
// [[Rcpp::export]]
IntegerVector IntBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return IntegerVector(rp->GetSubsetSize());
}

//' Returns an empty two-row integer matrix that ReadAlleles() can load to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Integer matrix with two rows, and appropriate number of columns for
//' ReadAlleles().
//' @export
// [[Rcpp::export]]
IntegerVector IntAlleleCodeBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return IntegerMatrix(2, rp->GetSubsetSize());
}

//' Returns a bool buffer that ReadAlleles() can load phasing information to.
//'
//' @param pgen Object returned by NewPgen().
//' @return Logical vector with appropriate length for ReadAlleles().
//' @export
// [[Rcpp::export]]
LogicalVector BoolBuf(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  return LogicalVector(rp->GetSubsetSize());
}

//' Loads the variant_numth variant, and then fills buf with \{0, 1, 2, NA\}
//' values indicating the number of copies of the first ALT (or user-specified)
//' allele each sample has.
//'
//' This function treats the data as diploid; you can divide by 2, and then
//' treat 0.5 as NA, if it's actually haploid.
//'
//' @param pgen Object returned by NewPgen().
//' @param buf Buffer returned by Buf() or IntBuf().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' to 2.
//' @return No return value, called for buf-filling side-effect.
//' @export
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

//' Returns whether hardcalls for the variant_numth variant and given allele
//' are represented in a sparse manner that is supported by
//' ReadSparseHardcalls().
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' to 2.
//' @return True iff the (variant, allele) pair has a sparse representation
//' that can be returned by ReadSparseHardcalls().
//' @export
// [[Rcpp::export]]
bool HasSparseHardcalls(List pgen, int variant_num, int allele_num = 2) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  const uint32_t vrtype = rp->GetVrtype(variant_idx);
  // Don't support multiallelic variants outside the trivial REF case for now.
  if (allele_num == 1) {
    return ((vrtype & 7) == 6);
  } else if (allele_num == 2) {
    return ((vrtype & 15) == 4);
  } else {
    return false;
  }
}

//' If HasSparseHardcalls() is true, returns a sparse representation for the
//' (variant, allele) pair.  If HasSparseHardcalls() is false, the function
//' fails.
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' to 2.
//' @param return_ints Whether to make the "counts" component of the return
//' value an IntegerVector instead of a NumericVector; defaults to false.
//' @return An object where "sample_nums" is an increasing sequence of positive
//' integers listing which samples have the allele, and "counts" is a vector
//' listing the allele counts for those samples.
//' @export
// [[Rcpp::export]]
List ReadSparseHardcalls(List pgen, int variant_num, int allele_num = 2, bool return_ints = false) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  const uint32_t vrtype = rp->GetVrtype(variant_idx);
  const bool is_supported_sparse =
    ((allele_num == 1) && ((vrtype & 7) == 6)) ||
    ((allele_num == 2) && ((vrtype & 15) == 4));
  if (!is_supported_sparse) {
    stop("(variant, allele) does not have supported sparse representation");
  }
  IntegerVector sample_nums(0);
  const int allele_idx = allele_num - 1;
  const uint32_t file_sample_ct = rp->GetRawSampleCt();
  // Should be impossible for this to fail.
  const uint32_t max_difflist_len = 2 * (file_sample_ct / plink2::kPglMaxDifflistLenDivisor);
  if (return_ints) {
    IntegerVector unused_buf(0);
    IntegerVector integer_counts(0);
    rp->ReadIntMaybeSparseHardcalls(unused_buf, variant_idx, allele_idx, max_difflist_len, &sample_nums, &integer_counts);
    return List::create(
                        _["sample_nums"] = sample_nums,
                        _["counts"] = integer_counts
                        );
  } else {
    NumericVector unused_buf(0);
    NumericVector numeric_counts(0);
    rp->ReadMaybeSparseHardcalls(unused_buf, variant_idx, allele_idx, max_difflist_len, &sample_nums, &numeric_counts);
    return List::create(
                        _["sample_nums"] = sample_nums,
                        _["counts"] = numeric_counts
                        );
  }
}

//' Loads the variant_numth variant, and then fills buf with numeric dosages
//' in [0, 2] indicating the dosages of the first ALT (or user-specified)
//' allele for each sample, with missing values represented by NA.
//'
//' This function treats the data as diploid; divide by 2 to obtain haploid
//' dosages.
//'
//' @param pgen Object returned by NewPgen().
//' @param buf Buffer returned by Buf().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' to 2.
//' @return No return value, called for buf-filling side-effect.
//' @export
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

//' Returns whether dosages for the variant_numth variant and given allele
//' are represented in a sparse manner that is supported by
//' ReadSparse(), under the current sample subset.
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' to 2.
//' @return True iff the (variant, allele) pair has a sparse representation
//' that can be returned by ReadSparse().
//' @export
// [[Rcpp::export]]
bool HasSparse(List pgen, int variant_num, int allele_num = 2) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  const uint32_t vrtype = rp->GetVrtype(variant_idx);
  // Don't support multiallelic variants outside the trivial REF case for now.
  if (allele_num == 1) {
    if ((vrtype & 7) != 6) {
      return false;
    }
  } else if (allele_num == 2) {
    if ((vrtype & 15) != 4) {
      return false;
    }
  } else {
    return false;
  }
  const uint32_t vrtype_dosage = vrtype & 0x60;
  if (vrtype_dosage != 0x20) {
    return (vrtype_dosage == 0);
  }
  // Dosage-list.  If hardcall-phase information is present and we're
  // extracting a sample-subset, the current implementation doesn't support
  // sparse return.  Otherwise we're ok.
  return (!(vrtype & 0x10)) || (rp->GetRawSampleCt() == rp->GetSubsetSize());
}

//' If HasSparse() is true, returns a sparse representation for the
//' (variant, allele) pair.  If HasSparse() is false, the function fails.
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_num Variant index (1-based).
//' @param allele_num Allele index; 1 corresponds to REF, 2 to the first ALT
//' allele, 3 to the second ALT allele if it exists, etc.  Optional, defaults
//' to 2.
//' @return An object where "sample_nums" is an increasing sequence of positive
//' integers listing which samples have the allele, and "dosages" is a vector
//' listing the dosages (on a 0-2 scale) for those samples.
//' @export
// [[Rcpp::export]]
List ReadSparse(List pgen, int variant_num, int allele_num = 2) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  const int variant_idx = variant_num - 1;
  const uint32_t vrtype = rp->GetVrtype(variant_idx);
  const int allele_idx = allele_num - 1;
  const uint32_t file_sample_ct = rp->GetRawSampleCt();
  // Keep this in sync with HasSparse().
  if (allele_idx == 0) {
    if ((vrtype & 7) != 6) {
    ReadSparse_fail:
      stop("(variant, allele) does not have supported sparse representation");
    }
  } else if (allele_idx == 1) {
    if ((vrtype & 15) != 4) {
      goto ReadSparse_fail;
    }
  } else {
    goto ReadSparse_fail;
  }
  const uint32_t vrtype_dosage = vrtype & 0x60;
  if (vrtype_dosage == 0x20) {
    if ((vrtype & 0x10) && (file_sample_ct != rp->GetSubsetSize())) {
      goto ReadSparse_fail;
    }
  } else if (vrtype_dosage != 0) {
    goto ReadSparse_fail;
  }

  IntegerVector sample_nums(0);

  // Should be impossible for this to fail.
  const uint32_t max_difflist_len = 3 * (file_sample_ct / plink2::kPglMaxDifflistLenDivisor);

  NumericVector unused_buf(0);
  NumericVector dosages(0);
  rp->ReadMaybeSparse(unused_buf, variant_idx, allele_idx, max_difflist_len, &sample_nums, &dosages);
  return List::create(
                      _["sample_nums"] = sample_nums,
                      _["dosages"] = dosages
                      );
}

//' Loads the variant_numth variant, and then fills acbuf with integer allele
//' codes, where each column of the buffer corresponds to a sample.  An allele
//' code of 0 corresponds to the REF allele, 1 to the first ALT, 2 to the
//' second ALT, etc.  Missing hardcalls are represented by a pair of NA codes.
//'
//' This function treats the data as diploid.  If it's really haploid, you may
//' want to compare the two rows, and then treat samples where the allele codes
//' differ as missing values.
//'
//' @param pgen Object returned by NewPgen().
//' @param acbuf Buffer returned by AlleleCodeBuf() or IntAlleleCodeBuf().
//' @param variant_num Variant index (1-based).
//' @param phasepresent_buf Buffer returned by BoolBuf().  Optional; if
//' provided, elements are set to true when the sample has known phase.  Most
//' of these values will be TRUE even when the raw data is unphased, because
//' homozygous genotypes always have known phase.  (Missing genotypes are
//' considered to have unknown phase.)
//' @return No return value, called for acbuf-filling side-effect.
//' @export
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

//' Load hardcalls for multiple variants as an integer matrix.
//'
//' This function treats the data as diploid; you can divide by 2, and then
//' treat 0.5 as NA, if it's actually haploid.
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_subset Integer vector containing 1-based indexes of variants
//' to load.
//' @return Integer matrix, where rows correspond to samples, columns
//' correspond to variant_subset, and values are in \{0, 1, 2, NA\} indicating
//' the number of hardcall ALT allele copies.  For multiallelic variants, all
//' ALT alleles are combined.
//' @export
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

//' Load dosages for multiple variants as a numeric matrix.
//'
//' This function treats the data as diploid; divide by 2 to obtain haploid
//' dosages.
//'
//' @param pgen Object returned by NewPgen().
//' @param variant_subset Integer vector containing 1-based indexes of variants
//' to load.
//' @param meanimpute Optional; if true, missing values are mean-imputed
//' instead of being represented by NA.
//' @return Numeric matrix, where rows correspond to samples, and columns
//' correspond to variant_subset.  Values are in [0, 2] indicating ALT
//' allele dosages, or NA for missing dosages.  For multiallelic variants, all
//' ALT alelles are combined.
//' @export
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

//' Compute variant scores.
//'
//' This function treats the data as diploid; divide by 2 to obtain scores
//' based on a haploid dosage matrix.
//'
//' @param pgen Object returned by NewPgen().
//' @param weights Sample weights.
//' @param variant_subset Integer vector containing 1-based indexes of variants
//' to include in the dosage matrix.  Optional; by default, all variants are
//' included.
//' @return Numeric vector, containing product of sample-weight vector and the
//' specified subset of the dosage matrix.
//' @export
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

//' Closes a pgen object, releasing resources.
//'
//' @param pgen Object returned by NewPgen().
//' @return No return value, called for side-effect.
//' @export
// [[Rcpp::export]]
void ClosePgen(List pgen) {
  if (strcmp_r_c(pgen[0], "pgen")) {
    stop("pgen is not a pgen object");
  }
  XPtr<class RPgenReader> rp = as<XPtr<class RPgenReader> >(pgen[1]);
  rp->Close();
}

class RPgenWriter {
public:
  // similar to PgenWriter in Python/pgenlib.pyx
  RPgenWriter();

#if __cplusplus >= 201103L
  RPgenWriter(const RPgenWriter&) = delete;
  RPgenWriter& operator=(const RPgenWriter&) = delete;
#endif

  void Open(String filename, int sample_ct, int variant_ct,
            SEXP nonref_flags, int allele_ct_limit,
            bool hardcall_phase_present, bool dosage_present);

  uint32_t GetWrittenVariantCt() const;

  void AppendBiallelic(SEXP geno);

  void AppendAlleles(SEXP acbuf, Nullable<LogicalVector> phasepresent,
                     bool all_phased, Nullable<int> allele_ct);

  void AppendDosages(NumericVector dosages);

  void Close();

  ~RPgenWriter();

private:
  plink2::STPgenWriter* _state_ptr;
  unsigned char* _spgw_alloc;
  uintptr_t* _nonref_flags;
  plink2::PgenGlobalFlags _phase_dosage_gflags;
  uint32_t _allele_ct_limit;
  uintptr_t* _genovec;
  uintptr_t* _patch_01_set;
  plink2::AlleleCode* _patch_01_vals;
  uintptr_t* _patch_10_set;
  plink2::AlleleCode* _patch_10_vals;
  uintptr_t* _phasepresent;
  uintptr_t* _phaseinfo;
  uintptr_t* _dosage_present;
  uint16_t* _dosage_main;
  int32_t* _allele_codes;
  unsigned char* _phasepresent_bytes;

  void CheckAppendOk() const;

  // Returns SpgwFinish() error code if finish is true, otherwise
  // kPglRetSuccess.  Never throws.
  plink2::PglErr Release(bool finish);
};

RPgenWriter::RPgenWriter() : _state_ptr(nullptr), _spgw_alloc(nullptr),
                             _nonref_flags(nullptr) {
}

void RPgenWriter::Open(String filename, int sample_ct, int variant_ct,
                       SEXP nonref_flags, int allele_ct_limit,
                       bool hardcall_phase_present, bool dosage_present) {
  if ((sample_ct <= 0) || (static_cast<uint32_t>(sample_ct) > plink2::kPglMaxSampleCt)) {
    stop("sample_ct must be positive, and less than ~2^31");
  }
  if ((variant_ct <= 0) || (static_cast<uint32_t>(variant_ct) > plink2::kPglMaxVariantCt)) {
    stop("variant_ct must be positive, and less than ~2^31");
  }
  if ((allele_ct_limit < 2) || (static_cast<uint32_t>(allele_ct_limit) > plink2::kPglMaxAlleleCt)) {
    stop("allele_ct_limit must be in [2, 255]");
  }
  const uint32_t sample_ct_u = sample_ct;
  const uint32_t variant_ct_u = variant_ct;
  // Same choice as the Python PgenWriter: the multiallelic case uses
  // write-and-copy mode, since the variant record width isn't known upfront.
  const plink2::PgenWriteMode write_mode = (allele_ct_limit == 2)? plink2::kPgenWriteBackwardSeek : plink2::kPgenWriteAndCopy;
  uint32_t nonref_flags_storage = 0;
  if (!Rf_isNull(nonref_flags)) {
    if (TYPEOF(nonref_flags) != LGLSXP) {
      stop("nonref_flags must be NULL, TRUE, FALSE, or a logical vector of length variant_ct");
    }
    LogicalVector nonref_lgl(nonref_flags);
    if ((nonref_lgl.size() == 1) && (variant_ct != 1)) {
      if (nonref_lgl[0] == NA_LOGICAL) {
        stop("nonref_flags cannot be NA");
      }
      nonref_flags_storage = nonref_lgl[0]? 2 : 1;
    } else if (nonref_lgl.size() == variant_ct) {
      for (uint32_t variant_idx = 0; variant_idx != variant_ct_u; ++variant_idx) {
        if (nonref_lgl[variant_idx] == NA_LOGICAL) {
          stop("nonref_flags cannot contain NA values");
        }
      }
      nonref_flags_storage = 3;
      const uintptr_t nonref_flags_cacheline_ct = plink2::DivUp(variant_ct_u, plink2::kBitsPerCacheline);
      if (plink2::cachealigned_malloc(nonref_flags_cacheline_ct * plink2::kCacheline, &_nonref_flags)) {
        stop("Out of memory");
      }
      plink2::ZeroWArr(nonref_flags_cacheline_ct * plink2::kWordsPerCacheline, _nonref_flags);
      for (uint32_t variant_idx = 0; variant_idx != variant_ct_u; ++variant_idx) {
        if (nonref_lgl[variant_idx]) {
          plink2::SetBit(variant_idx, _nonref_flags);
        }
      }
    } else {
      stop("nonref_flags must be NULL, TRUE, FALSE, or a logical vector of length variant_ct");
    }
  }
  plink2::PgenGlobalFlags phase_dosage_gflags = plink2::kfPgenGlobal0;
  if (hardcall_phase_present) {
    phase_dosage_gflags |= plink2::kfPgenGlobalHardcallPhasePresent;
  }
  if (dosage_present) {
    phase_dosage_gflags |= plink2::kfPgenGlobalDosagePresent;
  }
  _state_ptr = static_cast<plink2::STPgenWriter*>(malloc(sizeof(plink2::STPgenWriter)));
  if (!_state_ptr) {
    Release(false);
    stop("Out of memory");
  }
  plink2::PreinitSpgw(_state_ptr);
  const char* fname = filename.get_cstring();
  uintptr_t alloc_cacheline_ct;
  uint32_t max_vrec_len;
  plink2::PglErr reterr = plink2::SpgwInitPhase1(fname, nullptr, _nonref_flags, variant_ct_u, sample_ct_u, allele_ct_limit, write_mode, phase_dosage_gflags, nonref_flags_storage, _state_ptr, &alloc_cacheline_ct, &max_vrec_len);
  if (reterr != plink2::kPglRetSuccess) {
    Release(false);
    char errstr_buf[256];
    if (reterr == plink2::kPglRetOpenFail) {
      snprintf(errstr_buf, 256, "Failed to open %s for writing", fname);
    } else {
      snprintf(errstr_buf, 256, "SpgwInitPhase1() error %d", static_cast<int>(reterr));
    }
    stop(errstr_buf);
  }
  const uintptr_t bitvec_cacheline_ct = plink2::DivUp(sample_ct_u, plink2::kBitsPerCacheline);
  const uintptr_t genovec_cacheline_ct = plink2::DivUp(sample_ct_u, plink2::kNypsPerCacheline);
  const uintptr_t patch_01_vals_cacheline_ct = plink2::DivUp(sample_ct_u * sizeof(plink2::AlleleCode), plink2::kCacheline);
  const uintptr_t patch_10_vals_cacheline_ct = plink2::DivUp(sample_ct_u * 2 * sizeof(plink2::AlleleCode), plink2::kCacheline);
  const uintptr_t dosage_main_cacheline_ct = plink2::DivUp(sample_ct_u, 2 * plink2::kInt32PerCacheline);
  const uintptr_t allele_codes_cacheline_ct = plink2::DivUp(sample_ct_u * 2 * sizeof(int32_t), plink2::kCacheline);
  const uintptr_t phasepresent_bytes_cacheline_ct = plink2::DivUp(sample_ct_u, plink2::kCacheline);
  const uintptr_t buf_cacheline_ct = genovec_cacheline_ct + 5 * bitvec_cacheline_ct + patch_01_vals_cacheline_ct + patch_10_vals_cacheline_ct + dosage_main_cacheline_ct + allele_codes_cacheline_ct + phasepresent_bytes_cacheline_ct;
  if (plink2::cachealigned_malloc((alloc_cacheline_ct + buf_cacheline_ct) * plink2::kCacheline, &_spgw_alloc)) {
    Release(false);
    stop("Out of memory");
  }
  plink2::SpgwInitPhase2(max_vrec_len, _state_ptr, _spgw_alloc);
  unsigned char* spgw_alloc_iter = &(_spgw_alloc[alloc_cacheline_ct * plink2::kCacheline]);
  // Trailing bits of the genotype and phasepresent buffers must be clear.
  memset(spgw_alloc_iter, 0, buf_cacheline_ct * plink2::kCacheline);
  _genovec = reinterpret_cast<uintptr_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[genovec_cacheline_ct * plink2::kCacheline]);
  _patch_01_set = reinterpret_cast<uintptr_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
  _patch_01_vals = reinterpret_cast<plink2::AlleleCode*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[patch_01_vals_cacheline_ct * plink2::kCacheline]);
  _patch_10_set = reinterpret_cast<uintptr_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
  _patch_10_vals = reinterpret_cast<plink2::AlleleCode*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[patch_10_vals_cacheline_ct * plink2::kCacheline]);
  _phasepresent = reinterpret_cast<uintptr_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
  _phaseinfo = reinterpret_cast<uintptr_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
  _dosage_present = reinterpret_cast<uintptr_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[bitvec_cacheline_ct * plink2::kCacheline]);
  _dosage_main = reinterpret_cast<uint16_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[dosage_main_cacheline_ct * plink2::kCacheline]);
  _allele_codes = reinterpret_cast<int32_t*>(spgw_alloc_iter);
  spgw_alloc_iter = &(spgw_alloc_iter[allele_codes_cacheline_ct * plink2::kCacheline]);
  _phasepresent_bytes = spgw_alloc_iter;
  _phase_dosage_gflags = phase_dosage_gflags;
  _allele_ct_limit = allele_ct_limit;
}

void RPgenWriter::CheckAppendOk() const {
  if (!_state_ptr) {
    stop("pgen_writer is closed");
  }
  if (plink2::SpgwGetVidx(_state_ptr) == plink2::SpgwGetVariantCt(_state_ptr)) {
    stop("all variant_ct variants have already been written");
  }
}

uint32_t RPgenWriter::GetWrittenVariantCt() const {
  if (!_state_ptr) {
    stop("pgen_writer is closed");
  }
  return plink2::SpgwGetVidx(_state_ptr);
}

template <class T> static void CheckWriterBufLen(const T& vec, uint32_t expected_len, const char* vec_name) {
  if (static_cast<uintptr_t>(vec.size()) != expected_len) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "%s has wrong length (%lld; %u expected)", vec_name, static_cast<long long>(vec.size()), expected_len);
    stop(errstr_buf);
  }
}

void RPgenWriter::AppendBiallelic(SEXP geno) {
  CheckAppendOk();
  const uint32_t sample_ct = plink2::SpgwGetSampleCt(_state_ptr);
  uintptr_t* genovec = _genovec;
  const uint32_t word_ct = plink2::NypCtToWordCt(sample_ct);
  plink2::ZeroWArr(word_ct, genovec);
  if (TYPEOF(geno) == INTSXP) {
    IntegerVector geno_int(geno);
    CheckWriterBufLen(geno_int, sample_ct, "geno");
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const int cur_val = geno_int[sample_idx];
      uintptr_t cur_geno;
      if ((cur_val >= 0) && (cur_val <= 2)) {
        cur_geno = cur_val;
      } else if ((cur_val == NA_INTEGER) || (cur_val == -9)) {
        cur_geno = 3;
      } else {
        stop("geno values must be in {0, 1, 2, NA, -9}");
      }
      genovec[sample_idx / plink2::kBitsPerWordD2] |= cur_geno << (2 * (sample_idx % plink2::kBitsPerWordD2));
    }
  } else if (TYPEOF(geno) == REALSXP) {
    NumericVector geno_dbl(geno);
    CheckWriterBufLen(geno_dbl, sample_ct, "geno");
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const double cur_val = geno_dbl[sample_idx];
      uintptr_t cur_geno;
      if ((cur_val == 0.0) || (cur_val == 1.0) || (cur_val == 2.0)) {
        cur_geno = static_cast<uintptr_t>(cur_val);
      } else if (ISNAN(cur_val) || (cur_val == -9.0)) {
        cur_geno = 3;
      } else {
        stop("geno values must be in {0, 1, 2, NA, -9}");
      }
      genovec[sample_idx / plink2::kBitsPerWordD2] |= cur_geno << (2 * (sample_idx % plink2::kBitsPerWordD2));
    }
  } else {
    stop("Unsupported geno type (integer or numeric vector expected)");
  }
  const plink2::PglErr reterr = plink2::SpgwAppendBiallelicGenovec(genovec, _state_ptr);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "SpgwAppendBiallelicGenovec() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

void RPgenWriter::AppendAlleles(SEXP acbuf, Nullable<LogicalVector> phasepresent,
                                bool all_phased, Nullable<int> allele_ct) {
  CheckAppendOk();
  const uint32_t sample_ct = plink2::SpgwGetSampleCt(_state_ptr);
  if (phasepresent.isNotNull() && all_phased) {
    stop("phasepresent and all_phased cannot both be specified");
  }
  if ((phasepresent.isNotNull() || all_phased) && (!(_phase_dosage_gflags & plink2::kfPgenGlobalHardcallPhasePresent))) {
    stop("phased calls cannot be appended when pgen_writer was created with hardcall_phase_present = FALSE");
  }
  // Translate to the -9-for-missing int32 representation expected by
  // ConvertMultiAlleleCodesUnsafe().
  int32_t* allele_codes = _allele_codes;
  const uint32_t allele_code_ct = 2 * sample_ct;
  if (TYPEOF(acbuf) == INTSXP) {
    IntegerVector acbuf_int(acbuf);
    CheckWriterBufLen(acbuf_int, allele_code_ct, "acbuf");
    for (uint32_t uii = 0; uii != allele_code_ct; ++uii) {
      const int cur_val = acbuf_int[uii];
      if ((cur_val == NA_INTEGER) || (cur_val == -9)) {
        allele_codes[uii] = -9;
      } else if ((cur_val < 0) || (static_cast<uint32_t>(cur_val) >= _allele_ct_limit)) {
        stop("acbuf allele codes must be in [0, allele_ct_limit - 1], NA, or -9");
      } else {
        allele_codes[uii] = cur_val;
      }
    }
  } else if (TYPEOF(acbuf) == REALSXP) {
    NumericVector acbuf_dbl(acbuf);
    CheckWriterBufLen(acbuf_dbl, allele_code_ct, "acbuf");
    for (uint32_t uii = 0; uii != allele_code_ct; ++uii) {
      const double cur_val = acbuf_dbl[uii];
      if (ISNAN(cur_val) || (cur_val == -9.0)) {
        allele_codes[uii] = -9;
      } else if ((cur_val < 0.0) || (cur_val >= static_cast<double>(_allele_ct_limit)) || (cur_val != static_cast<double>(static_cast<int32_t>(cur_val)))) {
        stop("acbuf allele codes must be in [0, allele_ct_limit - 1], NA, or -9");
      } else {
        allele_codes[uii] = static_cast<int32_t>(cur_val);
      }
    }
  } else {
    stop("Unsupported acbuf type (integer or numeric matrix expected)");
  }
  unsigned char* phasepresent_bytes = nullptr;
  uintptr_t* phasepresent_buf = nullptr;
  uintptr_t* phaseinfo = nullptr;
  if (phasepresent.isNotNull()) {
    LogicalVector phasepresent_lgl(phasepresent.get());
    CheckWriterBufLen(phasepresent_lgl, sample_ct, "phasepresent");
    phasepresent_bytes = _phasepresent_bytes;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const int cur_val = phasepresent_lgl[sample_idx];
      phasepresent_bytes[sample_idx] = (cur_val != NA_LOGICAL) && cur_val;
    }
    phasepresent_buf = _phasepresent;
    phaseinfo = _phaseinfo;
  } else if (all_phased) {
    phaseinfo = _phaseinfo;
  }
  uint32_t patch_01_ct;
  uint32_t patch_10_ct;
  const int32_t observed_allele_ct = plink2::ConvertMultiAlleleCodesUnsafe(allele_codes, phasepresent_bytes, sample_ct, _genovec, _patch_01_set, _patch_01_vals, _patch_10_set, _patch_10_vals, &patch_01_ct, &patch_10_ct, phasepresent_buf, phaseinfo);
  if (observed_allele_ct == -1) {
    stop("invalid allele codes (missing codes must occur in pairs)");
  }
  uint32_t write_allele_ct = observed_allele_ct;
  if (allele_ct.isNotNull()) {
    const int allele_ct_int = as<int>(allele_ct.get());
    if ((allele_ct_int < 2) || (static_cast<uint32_t>(allele_ct_int) > _allele_ct_limit)) {
      stop("allele_ct must be in [2, allele_ct_limit]");
    }
    if (static_cast<uint32_t>(allele_ct_int) < write_allele_ct) {
      stop("acbuf contains allele codes >= allele_ct");
    }
    write_allele_ct = allele_ct_int;
  }
  plink2::PglErr reterr;
  const bool is_biallelic = (patch_01_ct == 0) && (patch_10_ct == 0);
  if (!phaseinfo) {
    if (is_biallelic) {
      reterr = plink2::SpgwAppendBiallelicGenovec(_genovec, _state_ptr);
    } else {
      reterr = plink2::SpgwAppendMultiallelicSparse(_genovec, _patch_01_set, _patch_01_vals, _patch_10_set, _patch_10_vals, write_allele_ct, patch_01_ct, patch_10_ct, _state_ptr);
    }
  } else {
    if (is_biallelic) {
      reterr = plink2::SpgwAppendBiallelicGenovecHphase(_genovec, phasepresent_buf, phaseinfo, _state_ptr);
    } else {
      reterr = plink2::SpgwAppendMultiallelicGenovecHphase(_genovec, _patch_01_set, _patch_01_vals, _patch_10_set, _patch_10_vals, phasepresent_buf, phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, _state_ptr);
    }
  }
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "AppendAlleles() write error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

void RPgenWriter::AppendDosages(NumericVector dosages) {
  CheckAppendOk();
  if (!(_phase_dosage_gflags & plink2::kfPgenGlobalDosagePresent)) {
    stop("dosages cannot be appended when pgen_writer was created with dosage_present = FALSE");
  }
  const uint32_t sample_ct = plink2::SpgwGetSampleCt(_state_ptr);
  CheckWriterBufLen(dosages, sample_ct, "dosages");
  const double* dosages_ptr = &dosages[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double cur_val = dosages_ptr[sample_idx];
    if ((!ISNAN(cur_val)) && (!((cur_val >= 0.0) && (cur_val <= 2.0)))) {
      stop("dosages must be in [0, 2], or NA");
    }
  }
  uint32_t dosage_ct;
  // hard_call_halfdist = 6554 corresponds to plink2's default hard-call
  // threshold of 0.1; same as the Python PgenWriter.
  plink2::DoublesToDosage16(dosages_ptr, sample_ct, 6554, _genovec, _dosage_present, _dosage_main, &dosage_ct);
  const plink2::PglErr reterr = plink2::SpgwAppendBiallelicGenovecDosage16(_genovec, _dosage_present, _dosage_main, dosage_ct, _state_ptr);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "SpgwAppendBiallelicGenovecDosage16() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

plink2::PglErr RPgenWriter::Release(bool finish) {
  plink2::PglErr reterr = plink2::kPglRetSuccess;
  if (_state_ptr) {
    if (finish) {
      reterr = plink2::SpgwFinish(_state_ptr);
    }
    // SpgwFinish() closes the file(s) on success; CleanupSpgw() handles the
    // incomplete and error cases.
    plink2::CleanupSpgw(_state_ptr, &reterr);
    free(_state_ptr);
    _state_ptr = nullptr;
  }
  if (_spgw_alloc) {
    plink2::aligned_free(_spgw_alloc);
    _spgw_alloc = nullptr;
  }
  if (_nonref_flags) {
    plink2::aligned_free(_nonref_flags);
    _nonref_flags = nullptr;
  }
  return reterr;
}

void RPgenWriter::Close() {
  if (!_state_ptr) {
    return;
  }
  const uint32_t vidx = plink2::SpgwGetVidx(_state_ptr);
  const uint32_t variant_ct = plink2::SpgwGetVariantCt(_state_ptr);
  if (vidx != variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "ClosePgenWriter() called when number of written variants (%u) is unequal to declared variant_ct (%u)", vidx, variant_ct);
    stop(errstr_buf);
  }
  const plink2::PglErr reterr = Release(true);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "SpgwFinish() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
}

RPgenWriter::~RPgenWriter() {
  // Runs when the pgen_writer object is garbage-collected.  As with the Python
  // PgenWriter, a complete file is finalized, while an incomplete one is
  // closed without its header being filled in.
  const bool finish = _state_ptr && (plink2::SpgwGetVidx(_state_ptr) == plink2::SpgwGetVariantCt(_state_ptr));
  Release(finish);
}

static XPtr<class RPgenWriter> GetPgenWriterXPtr(List pgen_writer) {
  if (strcmp_r_c(pgen_writer[0], "pgen_writer")) {
    stop("pgen_writer is not a pgen_writer object");
  }
  return as<XPtr<class RPgenWriter> >(pgen_writer[1]);
}

//' Creates a new .pgen file for writing.
//'
//' Only the .pgen is written; you are responsible for writing the companion
//' .pvar (or .bim) and .psam (or .fam) files yourself, with samples in the
//' same order as in the vectors passed to the append functions, and variants
//' in append order.  After all variant_ct variants have been appended (with
//' AppendBiallelic(), AppendAlleles(), or AppendDosages()), ClosePgenWriter()
//' must be called to finalize the file.
//'
//' @param filename .pgen file path.
//' @param sample_ct Number of samples.
//' @param variant_ct Number of variants that will be written.
//' @param nonref_flags TRUE when the data is from an ordinary PLINK 1 .bed
//' file (where the A2 allele is major rather than consistently reference),
//' FALSE when the REF allele is consistently reference, or a logical vector of
//' length variant_ct when this is mixed.  NULL is also permitted (this
//' delegates tracking of nonref information to the .pvar file), but it's
//' discouraged since .pgen+.bim+.fam is a useful data representation with
//' direct plink2 support.  Defaults to TRUE, as in the Python PgenWriter.
//' @param allele_ct_limit Maximum number of alleles (REF + ALTs) that any
//' variant can have; defaults to 2.  Values larger than 2 are only useful with
//' AppendAlleles().
//' @param hardcall_phase_present Whether phased hardcalls may be written (with
//' AppendAlleles()).  Defaults to FALSE.
//' @param dosage_present Whether dosages may be written (with
//' AppendDosages()).  Defaults to FALSE.
//' @return A pgen_writer object.
//' @export
// [[Rcpp::export]]
SEXP NewPgenWriter(String filename, int sample_ct, int variant_ct,
                   SEXP nonref_flags = LogicalVector::create(TRUE), int allele_ct_limit = 2,
                   bool hardcall_phase_present = false,
                   bool dosage_present = false) {
  XPtr<class RPgenWriter> pgen_writer(new RPgenWriter(), true);
  pgen_writer->Open(filename, sample_ct, variant_ct, nonref_flags, allele_ct_limit, hardcall_phase_present, dosage_present);
  return List::create(_["class"] = "pgen_writer", _["pgen_writer"] = pgen_writer);
}

//' Appends an unphased biallelic hardcall variant to a .pgen being written.
//'
//' @param pgen_writer Object returned by NewPgenWriter().
//' @param geno Integer or numeric vector of length sample_ct, with values in
//' \{0, 1, 2, NA\} indicating the number of ALT allele copies; -9 is also
//' accepted as a missing-value code.  (A buffer filled by ReadHardcalls() is
//' in this format.)
//' @return No return value, called for side-effect.
//' @export
// [[Rcpp::export]]
void AppendBiallelic(List pgen_writer, SEXP geno) {
  GetPgenWriterXPtr(pgen_writer)->AppendBiallelic(geno);
}

//' Appends a hardcall variant, specified as allele codes, to a .pgen being
//' written.  Unlike AppendBiallelic(), this supports phased and multiallelic
//' hardcalls.
//'
//' @param pgen_writer Object returned by NewPgenWriter().
//' @param acbuf Integer or numeric matrix with 2 rows and sample_ct columns,
//' in the format filled by ReadAlleles(): 0 corresponds to the REF allele, 1
//' to the first ALT, 2 to the second ALT, etc., and a missing hardcall is
//' represented by a pair of NA (or -9) codes.  For a phased heterozygous call,
//' the first row holds the allele on the first haplotype; unphased
//' heterozygous calls can be stored in either order.
//' @param phasepresent Logical vector of length sample_ct indicating which
//' calls are phased; entries for non-heterozygous calls are ignored.
//' Optional.  (A buffer filled by ReadAlleles() is in this format.)  Requires
//' hardcall_phase_present = TRUE in the NewPgenWriter() call.
//' @param all_phased Whether to treat all calls as phased; defaults to FALSE.
//' Cannot be combined with phasepresent, and requires
//' hardcall_phase_present = TRUE in the NewPgenWriter() call.  If neither
//' phasepresent nor all_phased is specified, all calls are unphased.
//' @param allele_ct Number of alleles (REF + ALTs) for this variant.
//' Optional; by default, this is inferred to be max(2, 1 + <max observed
//' allele code>).  It must be specified when the last ALT allele in the .pvar
//' is unobserved.
//' @return No return value, called for side-effect.
//' @export
// [[Rcpp::export]]
void AppendAlleles(List pgen_writer, SEXP acbuf, Nullable<LogicalVector> phasepresent = R_NilValue, bool all_phased = false, Nullable<int> allele_ct = R_NilValue) {
  GetPgenWriterXPtr(pgen_writer)->AppendAlleles(acbuf, phasepresent, all_phased, allele_ct);
}

//' Appends a biallelic dosage variant to a .pgen being written.
//'
//' Dosages are stored with 1/16384 precision.  Hardcalls are saved as well:
//' each is the nearest integer when the dosage is within 0.1 of it, and
//' missing otherwise.
//'
//' @param pgen_writer Object returned by NewPgenWriter(), with
//' dosage_present = TRUE.
//' @param dosages Numeric vector of length sample_ct, with ALT allele dosages
//' in [0, 2], or NA for missing.  (A buffer filled by Read() is in this
//' format.)
//' @return No return value, called for side-effect.
//' @export
// [[Rcpp::export]]
void AppendDosages(List pgen_writer, NumericVector dosages) {
  GetPgenWriterXPtr(pgen_writer)->AppendDosages(dosages);
}

//' Returns the number of variants appended so far to a .pgen being written.
//'
//' @param pgen_writer Object returned by NewPgenWriter().
//' @return Number of variants appended so far.
//' @export
// [[Rcpp::export]]
int GetWrittenVariantCt(List pgen_writer) {
  return GetPgenWriterXPtr(pgen_writer)->GetWrittenVariantCt();
}

//' Finalizes a .pgen being written, and releases resources.  This must be
//' called after all variant_ct variants have been appended; the .pgen is not
//' valid before then.
//'
//' @param pgen_writer Object returned by NewPgenWriter().
//' @return No return value, called for side-effect.
//' @export
// [[Rcpp::export]]
void ClosePgenWriter(List pgen_writer) {
  GetPgenWriterXPtr(pgen_writer)->Close();
}
