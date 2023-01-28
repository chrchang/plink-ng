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
      snprintf(errbuf, 256, "variant_num out of range (%d; must be 1..%d)", variant_idx + 1, _mp.variant_ct);
    } else {
      strcpy(errbuf, "pvar closed");
    }
    stop(errbuf);
  }
  return _mp.variant_ids[variant_idx];
}

std::pair<std::multimap<const char*, int, classcomp>::iterator, std::multimap<const char*, int, classcomp>::iterator> RPvar::GetVariantsById(const char* id) {
  if (_nameToIdxs.empty()) {
    const uint32_t len = _mp.variant_ct;
    for (uint32_t variant_idx = 0; variant_idx != len; ++variant_idx) {
      _nameToIdxs.insert(std::pair<const char*, int>(_mp.variant_ids[variant_idx], variant_idx));
    }
  }
  return _nameToIdxs.equal_range(id);
}

uint32_t RPvar::GetAlleleCt(uint32_t variant_idx) const {
  if (variant_idx >= _mp.variant_ct) {
    char errstr_buf[256];
    snprintf(errstr_buf, 256, "variant_num out of range (%d; must be 1..%u)", variant_idx + 1, _mp.variant_ct);
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
      snprintf(errbuf, 256, "variant_num out of range (%d; must be 1..%d)", variant_idx + 1, _mp.variant_ct);
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
    snprintf(errbuf, 256, "allele_num out of range (%d; must be 1..%d)", allele_idx + 1, allele_ct);
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
  _nameToIdxs.clear();
  plink2::CleanupMinimalPvar(&_mp);
}

RPvar::~RPvar() {
  _nameToIdxs.clear();
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
IntegerVector GetVariantsById(List pvar, String id) {
  if (strcmp_r_c(pvar[0], "pvar")) {
    stop("pvar is not a pvar object");
  }
  XPtr<class RPvar> rp = as<XPtr<class RPvar> >(pvar[1]);
  std::pair<std::multimap<const char*, int, classcomp>::iterator, std::multimap<const char*, int, classcomp>::iterator> equal_range = rp->GetVariantsById(id.get_cstring());
  std::multimap<const char*, int, classcomp>::iterator i1 = equal_range.first;
  std::multimap<const char*, int, classcomp>::iterator i2 = equal_range.second;
  const uint32_t len = std::distance(i1, i2);
  IntegerVector iv = IntegerVector(len);
  for (uint32_t uii = 0; uii != len; ++uii) {
    iv[uii] = i1->second + 1;
    ++i1;
  }
  return iv;
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
