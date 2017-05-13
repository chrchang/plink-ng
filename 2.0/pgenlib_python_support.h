#ifndef __PGENLIB_PYTHON_SUPPORT_H__
#define __PGENLIB_PYTHON_SUPPORT_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "pgenlib_internal.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Could define a slightly-more-efficient version of this function which uses a
// missing code of 3 instead of -9.  But let's play well with existing scripts
// first.
void genoarr_to_bytes_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int8_t* genobytes);

void genoarr_to_int32s_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* geno_int32);

void genoarr_to_int64s_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int64_t* geno_int64);

// For Python interface, allele_codes is always int32_t.  Python programmers
// should not need to worry about whether pgenlib was compiled with 1-, 2-, or
// 4-byte alt_allele_ct_t.
void genoarr_to_allele_codes(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* allele_codes);

// phasebytes can be nullptr; if it isn't, entry is 1 iff genotype is an
// explicitly phased het, OR genotype is homozygous
// phasepresent cannot be nullptr
void genoarr_phased_to_allele_codes(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t sample_ct, uint32_t phasepresent_ct, unsigned char* phasebytes, int32_t* allele_codes);

// assumes transposed genoarr, phaseinfo
void genoarr_phased_to_hap_codes(const uintptr_t* genoarr, const uintptr_t* phaseinfo, uint32_t variant_batch_size, int32_t* hap0_codes_iter, int32_t* hap1_codes_iter);

void dosage16_to_floats_minus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_vals, uint32_t sample_ct, uint32_t dosage_ct, float* geno_float);

void dosage16_to_doubles_minus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_vals, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double);

// Does not zero out trailing bits of bitarr.
void bytes_to_bits_unsafe(const uint8_t* boolbytes, uint32_t sample_ct, uintptr_t* bitarr);

// Bottom 2 bits are extracted from every byte.  Conveniently, -9 and 3 are
// treated identically.
// Does not zero out trailing bits of genoarr.
void bytes_to_genoarr_unsafe(const int8_t* genobytes, uint32_t sample_ct, uintptr_t* genoarr);

// - If phasepresent_bytes is nullptr, phasepresent is not updated.  In this
//   case, phaseinfo is updated iff it's not nullptr.  It's okay for both
//   phasepresent and phaseinfo to be nullptr here.
// - Otherwise, phasepresent and phaseinfo are always updated; neither can be
//   nullptr.
void allele_codes_to_genoarr_unsafe(const int32_t* allele_codes, const unsigned char* phasepresent_bytes, uint32_t sample_ct, uintptr_t* genoarr, uintptr_t* phasepresent, uintptr_t* phaseinfo);

void floats_to_dosage16(const float* floatarr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr);

void doubles_to_dosage16(const double* doublearr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PGENLIB_PYTHON_SUPPORT_H__
