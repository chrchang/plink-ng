#include "pvar.h"  // includes Rcpp

RPvar::RPvar() {
  PreinitMinimalPvar(&_mp);
}

void RPvar::Load(String filename) {
  char errbuf[plink2::kPglErrstrBufBlen];
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

uint32_t RPvar::GetVariantCt() const {
  return _mp.variant_ct;
}

const char* RPvar::GetVariantId(uint32_t variant_idx) const {
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

uint32_t RPvar::GetAlleleCt(uint32_t variant_idx) const {
  if (variant_idx >= _mp.variant_ct) {
    char errstr_buf[256];
    sprintf(errstr_buf, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _mp.variant_ct);
    stop(errstr_buf);
  }
  if (!_mp.allele_idx_offsetsp) {
    return 2;
  }
  const uintptr_t* allele_idx_offsets = _mp.allele_idx_offsetsp->p;
  return allele_idx_offsets[variant_idx + 1] - allele_idx_offsets[variant_idx];
}

const char* RPvar::GetAlleleCode(uint32_t variant_idx, uint32_t allele_idx) const {
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

plink2::RefcountedWptr* RPvar::GetAlleleIdxOffsetsp() {
  if (_mp.allele_idx_offsetsp) {
    _mp.allele_idx_offsetsp->ref_ct += 1;
  }
  return _mp.allele_idx_offsetsp;
}

uint32_t RPvar::GetMaxAlleleCt() const {
  return _mp.max_allele_ct;
}

void RPvar::Close() {
  plink2::CleanupMinimalPvar(&_mp);
}

RPvar::~RPvar() {
  plink2::CleanupMinimalPvar(&_mp);
}

// [[Rcpp::export]]
SEXP NewPvar(String filename) {
  XPtr<class RPvar> pvar(new RPvar(), true);
  pvar->Load(filename);
  return List::create(_["class"] = "pvar", _["pvar"] = pvar);
}

// [[Rcpp::export]]
String GetVariantId(List pvar, int variant_num) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  String ss(rp->GetVariantId(variant_num - 1));
  return ss;
}

// [[Rcpp::export]]
String GetAlleleCode(List pvar, int variant_num, int allele_num) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  String ss(rp->GetAlleleCode(variant_num - 1, allele_num - 1));
  return ss;
}

// [[Rcpp::export]]
void ClosePvar(List pvar) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  rp->Close();
}
