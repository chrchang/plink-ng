#include "pgenlib_ffi_support.h"

#include <Rcpp.h>
using namespace Rcpp;

class RPgenReader {
public:
  // imitates Python/pgenlib.pyx
  RPgenReader(String filename, Nullable<int> raw_sample_ct,
              Nullable<int> variant_ct,
              Nullable<IntegerVector> sample_subset_1based);

  int GetRawSampleCt() const;

  int GetSubsetSize() const;

  int GetVariantCt() const;

  bool HardcallPhasePresent() const;

  void ReadIntHardcalls(IntegerVector buf, int variant_idx, int allele_idx);

  void ReadHardcalls(NumericVector buf, int variant_idx, int allele_idx);

  void Read(NumericVector buf, int variant_idx, int allele_idx);

  void Close();

  ~RPgenReader();

private:
  plink2::PgenFileInfo* _info_ptr;
  plink2::PgenReader* _state_ptr;
  uintptr_t* _subset_include_vec;
  uintptr_t* _subset_include_interleaved_vec;
  uint32_t* _subset_cumulative_popcounts;
  uint32_t _subset_size;
  // preallocate buffers we'll use repeatedly
  uintptr_t* _genovec;
  uintptr_t* _phasepresent;
  uintptr_t* _phaseinfo;
  uintptr_t* _dosage_present;
  uint16_t* _dosage_main;
  plink2::VecW* _transpose_batch_buf;
  // kPglQuaterTransposeBatch (= 256) variants at a time, and then transpose
  uintptr_t* _multivar_vmaj_geno_buf;
  uintptr_t* _multivar_vmaj_phasepresent_buf;
  uintptr_t* _multivar_vmaj_phaseinfo_buf;
  uintptr_t* _multivar_smaj_geno_batch_buf;
  uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
  uintptr_t* _multivar_smaj_phasepresent_batch_buf;

  void SetSampleSubsetInternal(IntegerVector sample_subset_1based);
};

RPgenReader::RPgenReader(String filename, Nullable<int> raw_sample_ct,
                         Nullable<int> variant_ct,
                         Nullable<IntegerVector> sample_subset_1based) {
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
  if (variant_ct.isNotNull()) {
    cur_variant_ct = as<int>(variant_ct.get());
  }
  const char* fname = filename.get_cstring();
  plink2::PgenHeaderCtrl header_ctrl;
  uintptr_t pgfi_alloc_cacheline_ct;
  char errstr_buf[plink2::kPglErrstrBufBlen];
  if (PgfiInitPhase1(fname, cur_variant_ct, cur_sample_ct, 0, &header_ctrl, _info_ptr, &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
    stop(&(errstr_buf[7]));
  }
  if (header_ctrl & 0x30) {
    stop("Explicit ALT allele counts not yet supported");
  }
  if ((header_ctrl & 0xc0) == 0xc0) {
    stop("Explicit nonref_flags not yet supported");
  }
  const uint32_t file_sample_ct = _info_ptr->raw_sample_ct;
  unsigned char* pgfi_alloc = nullptr;
  if (plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
    stop("Out of memory");
  }
  uint32_t max_vrec_width;
  uintptr_t pgr_alloc_cacheline_ct;
  if (PgfiInitPhase2(header_ctrl, 1, 1, 0, 0, _info_ptr->raw_variant_ct, &max_vrec_width, _info_ptr, pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
    if (pgfi_alloc && (!_info_ptr->vrtypes)) {
      plink2::aligned_free(pgfi_alloc);
    }
    stop(&(errstr_buf[7]));
  }
  _state_ptr = static_cast<plink2::PgenReader*>(malloc(sizeof(plink2::PgenReader)));
  if (!_state_ptr) {
    stop("Out of memory");
  }
  plink2::PreinitPgr(_state_ptr);
  _state_ptr->fread_buf = nullptr;
  const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
  const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
  const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
  const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kQuatersPerVec) * plink2::kBytesPerVec;
  const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
  unsigned char* pgr_alloc;
  if (plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglQuaterTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglQuaterTransposeBatch) * genovec_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglQuaterTransposeBatch * plink2::kPglQuaterTransposeBatch / 8), &pgr_alloc)) {
    stop("Out of memory");
  }
  plink2::PglErr reterr = PgrInit(fname, max_vrec_width, _info_ptr, _state_ptr, pgr_alloc);
  if (reterr != plink2::kPglRetSuccess) {
    if (!_state_ptr->fread_buf) {
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

  // assumes kWordsPerVec <= 2
  _subset_include_interleaved_vec[-1] = 0;

  _subset_cumulative_popcounts = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
  _genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
  _phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
  _dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);
  if (sample_subset_1based.isNotNull()) {
    SetSampleSubsetInternal(sample_subset_1based.get());
  } else {
    _subset_size = file_sample_ct;
  }
  _transpose_batch_buf = reinterpret_cast<plink2::VecW*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
  _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglQuaterTransposeBatch * genovec_byte_ct]);
  _multivar_vmaj_phasepresent_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglQuaterTransposeBatch * sample_subset_byte_ct]);
  _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglQuaterTransposeBatch * sample_subset_byte_ct]);
  _multivar_smaj_geno_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglQuaterTransposeBatch * plink2::kPglQuaterTransposeBatch / 4]);
  _multivar_smaj_phaseinfo_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglQuaterTransposeBatch * plink2::kPglQuaterTransposeBatch / 8]);
  _multivar_smaj_phasepresent_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
  // pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglQuaterTransposeBatch * plink2::kPglQuaterTransposeBatch / 8]);
}

int RPgenReader::GetRawSampleCt() const {
  return _info_ptr->raw_sample_ct;
}

int RPgenReader::GetSubsetSize() const {
  return _subset_size;
}

int RPgenReader::GetVariantCt() const {
  return _info_ptr->raw_variant_ct;
}

bool RPgenReader::HardcallPhasePresent() const {
  return ((_info_ptr->gflags & plink2::kfPgenGlobalHardcallPhasePresent) != 0);
}

void RPgenReader::ReadIntHardcalls(IntegerVector buf, int variant_idx, int allele_idx) {
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
  plink2::PglErr reterr = PgrGet1(_subset_include_vec, _subset_cumulative_popcounts, _subset_size, variant_idx, allele_idx, _state_ptr, _genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "ReadIntHardcalls() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::GenoarrToInt32sMinus9(_genovec, _subset_size, &buf[0]);
}

void RPgenReader::ReadHardcalls(NumericVector buf, int variant_idx, int allele_idx) {
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
  plink2::PglErr reterr = PgrGet1(_subset_include_vec, _subset_cumulative_popcounts, _subset_size, variant_idx, allele_idx, _state_ptr, _genovec);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "ReadHardcalls() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::GenoarrToDoublesMinus9(_genovec, _subset_size, &buf[0]);
}

void RPgenReader::Read(NumericVector buf, int variant_idx, int allele_idx) {
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
  plink2::PglErr reterr = PgrGet1D(_subset_include_vec, _subset_cumulative_popcounts, _subset_size, variant_idx, allele_idx, _state_ptr, _genovec, _dosage_present, _dosage_main, &dosage_ct);
  if (reterr != plink2::kPglRetSuccess) {
    char errstr_buf[256];
    sprintf(errstr_buf, "Read() error %d", static_cast<int>(reterr));
    stop(errstr_buf);
  }
  plink2::Dosage16ToDoublesMinus9(_genovec, _dosage_present, _dosage_main, _subset_size, dosage_ct, &buf[0]);
}

void RPgenReader::Close() {
  // don't bother propagating file close errors for now
  if (_info_ptr) {
    plink2::CleanupPgfi(_info_ptr);
    if (_info_ptr->vrtypes) {
      plink2::aligned_free(_info_ptr->vrtypes);
      if (_state_ptr) {
        plink2::CleanupPgr(_state_ptr);
        if (_state_ptr->fread_buf) {
          plink2::aligned_free(_state_ptr->fread_buf);
        }
        free(_state_ptr);
        _state_ptr = nullptr;
      }
    }
    free(_info_ptr);
    _info_ptr = nullptr;
  }
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
  _subset_size = subset_size;
}

RPgenReader::~RPgenReader() {
  if (_info_ptr) {
    plink2::CleanupPgfi(_info_ptr);
    if (_info_ptr->vrtypes) {
      plink2::aligned_free(_info_ptr->vrtypes);
      if (_state_ptr) {
        plink2::CleanupPgr(_state_ptr);
        if (_state_ptr->fread_buf) {
          plink2::aligned_free(_state_ptr->fread_buf);
        }
        free(_state_ptr);
      }
    }
    free(_info_ptr);
  }
}

// [[Rcpp::export]]
SEXP NewReader(String filename, Nullable<int> raw_sample_ct = R_NilValue,
               Nullable<int> variant_ct = R_NilValue,
               Nullable<IntegerVector> sample_subset = R_NilValue) {
  XPtr<class RPgenReader> rp(new RPgenReader(filename, raw_sample_ct,
                                             variant_ct, sample_subset), true);
  return rp;
}

// [[Rcpp::export]]
int GetRawSampleCt(SEXP xp) {
  XPtr<class RPgenReader> rp(xp);
  return rp->GetRawSampleCt();
}

// [[Rcpp::export]]
int GetVariantCt(SEXP xp) {
  XPtr<class RPgenReader> rp(xp);
  return rp->GetVariantCt();
}

// [[Rcpp::export]]
bool HardcallPhasePresent(SEXP xp) {
  XPtr<class RPgenReader> rp(xp);
  return rp->HardcallPhasePresent();
}

// [[Rcpp::export]]
IntegerVector IntBuf(SEXP xp) {
  XPtr<class RPgenReader> rp(xp);
  return IntegerVector(rp->GetSubsetSize());
}

// [[Rcpp::export]]
NumericVector Buf(SEXP xp) {
  XPtr<class RPgenReader> rp(xp);
  return NumericVector(rp->GetSubsetSize());
}

// _num indicates R 1-based indexing.

// [[Rcpp::export]]
void ReadHardcalls(SEXP xp, SEXP buf, int variant_num, int allele_num = 2) {
  const int variant_idx = variant_num - 1;
  const int allele_idx = allele_num - 1;
  XPtr<class RPgenReader> rp(xp);
  if (TYPEOF(buf) == REALSXP) {
    rp->ReadHardcalls(buf, variant_idx, allele_idx);
  } else if (TYPEOF(buf) == INTSXP) {
    rp->ReadIntHardcalls(buf, variant_idx, allele_idx);
  } else {
    stop("Unsupported buf type");
  }
}

// [[Rcpp::export]]
void Read(SEXP xp, NumericVector buf, int variant_num, int allele_num = 2) {
  XPtr<class RPgenReader> rp(xp);
  rp->Read(buf, variant_num - 1, allele_num - 1);
}

// [[Rcpp::export]]
void Close(SEXP xp) {
  XPtr<class RPgenReader> rp(xp);
  rp->Close();
}
