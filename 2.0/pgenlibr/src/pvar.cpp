#include "pvar.h"  // includes Rcpp

RPvarReader::RPvarReader(String filename) {
  char errbuf[plink2::kPglErrstrBufBlen];
  PreinitMinimalPvar(&_mp);
  plink2::PglErr reterr = LoadMinimalPvar(filename.get_cstring(), &_mp, errbuf);
  if (reterr != plink2::kPglRetSuccess) {
    if (reterr == plink2::kPglRetNomem) {
      stop("Out of memory");
    } else if (reterr == plink2::kPglRetReadFail) {
      stop("File read failure");
    } else {
      stop(&errbuf[7]);
    }
  }
}

uint32_t RPvarReader::GetVariantCt() const {
  return _mp.variant_ct;
}

const char* RPvarReader::GetVariantId(uint32_t variant_idx) const {
  if (variant_idx >= _mp.variant_ct) {
    char errbuf[256];
    if (_mp.variant_ct) {
      sprintf(errbuf, "variant_num out of range (%d; must be 1..%d)", variant_idx + 1, _mp.variant_ct);
    } else {
      strcpy(errbuf, "pvar closed");
    }
    stop(errbuf);
  }
  return _mp.variant_ids[variant_idx];
}

uint32_t RPvarReader::GetAlleleCt(uint32_t variant_idx) const {
  if (!_mp.allele_idx_offsetsp) {
    return 2;
  } else {
    const uintptr_t* allele_idx_offsets = _mp.allele_idx_offsetsp->p;
    return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
  }
}

const char* RPvarReader::GetAlleleCode(uint32_t variant_idx, uint32_t allele_idx) const {
  if (variant_idx >= _mp.variant_ct) {
    char errbuf[256];
    if (_mp.variant_ct) {
      sprintf(errbuf, "variant_num out of range (%d; must be 1..%d)", variant_idx + 1, _mp.variant_ct);
    } else {
      strcpy(errbuf, "pvar closed");
    }
    stop(errbuf);
  }
  uintptr_t allele_idx_offset_base = 2 * variant_idx;
  uint32_t allele_ct = 2;
  if (_mp.allele_idx_offsetsp) {
    const uintptr_t* allele_idx_offsets = _mp.allele_idx_offsetsp->p;
    allele_idx_offset_base = allele_idx_offsets[variant_idx];
    allele_ct = allele_idx_offsets[variant_idx + 1] - allele_idx_offset_base;
  }
  if (allele_idx >= allele_ct) {
    char errbuf[256];
    sprintf(errbuf, "allele_num out of range (%d; must be 1..%d)", allele_idx + 1, allele_ct);
    stop(errbuf);
  }
  return _mp.allele_storage[allele_idx_offset_base + allele_idx];
}

plink2::RefcountedWptr* RPvarReader::GetAlleleIdxOffsetsp() {
  if (_mp.allele_idx_offsetsp) {
    _mp.allele_idx_offsetsp->ref_ct += 1;
  }
  return _mp.allele_idx_offsetsp;
}

uint32_t RPvarReader::GetMaxAlleleCt() const {
  return _mp.max_allele_ct;
}

void RPvarReader::Close() {
  plink2::CleanupMinimalPvar(&_mp);
}

RPvarReader::~RPvarReader() {
  plink2::CleanupMinimalPvar(&_mp);
}

// [[Rcpp::export]]
SEXP NewPvar(String filename) {
  XPtr<class RPvarReader> rp(new RPvarReader(filename), true);
  return rp;
}

// [[Rcpp::export]]
int GetPvarCt(SEXP pvar) {
  XPtr<class RPvarReader> rp(pvar);
  return rp->GetVariantCt();
}

// [[Rcpp::export]]
String GetVariantId(SEXP pvar, int variant_num) {
  XPtr<class RPvarReader> rp(pvar);
  String ss(rp->GetVariantId(variant_num - 1));
  return ss;
}

// [[Rcpp::export]]
int GetAlleleCt(SEXP pvar, int variant_num) {
  XPtr<class RPvarReader> rp(pvar);
  return rp->GetAlleleCt(variant_num - 1);
}

// [[Rcpp::export]]
String GetAlleleCode(SEXP pvar, int variant_num, int allele_num) {
  XPtr<class RPvarReader> rp(pvar);
  String ss(rp->GetAlleleCode(variant_num - 1, allele_num - 1));
  return ss;
}

// [[Rcpp::export]]
int GetMaxAlleleCt(SEXP pvar) {
  XPtr<class RPvarReader> rp(pvar);
  return rp->GetMaxAlleleCt();
}

// [[Rcpp::export]]
void ClosePvar(SEXP pvar) {
  XPtr<class RPvarReader> rp(pvar);
  rp->Close();
}
