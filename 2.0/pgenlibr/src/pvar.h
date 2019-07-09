#ifndef __PVAR_H__
#define __PVAR_H__

#include "pvar_ffi_support.h"

#include <Rcpp.h>
using namespace Rcpp;

class RPvar {
public:
  // only tracks variant IDs and allele codes for now
  RPvar();

#if __cplusplus >= 201103L
  RPvar(const RPvar&) = delete;
  RPvar& operator=(const RPvar&) = delete;
#endif

  void Load(String filename);

  uint32_t GetVariantCt() const;

  const char* GetVariantId(uint32_t variant_idx) const;

  uint32_t GetAlleleCt(uint32_t variant_idx) const;

  const char* GetAlleleCode(uint32_t variant_idx, uint32_t allele_idx) const;

  plink2::RefcountedWptr* GetAlleleIdxOffsetsp();

  uint32_t GetMaxAlleleCt() const;

  void Close();

  ~RPvar();

private:
  plink2::MinimalPvar _mp;
};

HEADER_INLINE int strcmp_r_c(String r_string, const char* cstr) {
  return strcmp(r_string.get_cstring(), cstr);
}

#endif  // __PVAR_H__
