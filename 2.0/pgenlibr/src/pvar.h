#ifndef __PVAR_H__
#define __PVAR_H__

#include "pvar_ffi_support.h"

#include <Rcpp.h>
using namespace Rcpp;

class RPvarReader {
public:
  // only tracks variant IDs and allele codes for now
  RPvarReader(String filename);

  uint32_t GetVariantCt() const;

  const char* GetVariantId(uint32_t variant_idx) const;

  uint32_t GetAlleleCt(uint32_t variant_idx) const;

  const char* GetAlleleCode(uint32_t variant_idx, uint32_t allele_idx) const;

  plink2::RefcountedWptr* GetAlleleIdxOffsetsp();

  uint32_t GetMaxAlleleCt() const;

  void Close();

  ~RPvarReader();

private:
  plink2::MinimalPvar _mp;
};

#endif  // __PVAR_H__
