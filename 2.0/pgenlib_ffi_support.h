#ifndef __PGENLIB_FFI_SUPPORT_H__
#define __PGENLIB_FFI_SUPPORT_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2023 Shaun Purcell,
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

#include "include/pgenlib_misc.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Could define a slightly-more-efficient version of this function which uses a
// missing code of 3 instead of -9.  But let's play well with existing scripts
// first.
void GenoarrToBytesMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int8_t* genobytes);

void GenoarrToInt32sMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* geno_int32);

void GenoarrToInt64sMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int64_t* geno_int64);

// May want to use STD_ARRAY_INIT_{START,END}... though it may not be worth the
// additional compilation headaches here.
extern const double kGenoDoublePairs[32];

HEADER_INLINE void GenoarrToDoublesMinus9(const uintptr_t* genoarr, uint32_t sample_ct, double* geno_double) {
  GenoarrLookup16x8bx2(genoarr, kGenoDoublePairs, sample_ct, geno_double);
}

HEADER_INLINE void GenoarrToAlleleCodes(const uint64_t* geno_to_intcode_pair_table, const uintptr_t* genoarr, uint32_t sample_ct, int32_t* allele_codes) {
  GenoarrLookup16x8bx2(genoarr, geno_to_intcode_pair_table, sample_ct, allele_codes);
}

extern const uint64_t kGenoToIntcodeDPairs[32];

// For FFI, allele_codes is always int32_t.  Python/R programmers should not
// need to worry about whether pgenlib was compiled with 1-, 2-, or 4-byte
// AlleleCode.
//
// phasebytes can be nullptr; if it isn't, entry is 1 iff genotype is an
// explicitly phased het, OR genotype is homozygous
// phasepresent cannot be nullptr
void GenoarrPhasedToAlleleCodes(const uint64_t* geno_to_intcode_dpair_table, const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t sample_ct, uint32_t phasepresent_ct, unsigned char* phasebytes, int32_t* allele_codes);

HEADER_INLINE void GenoarrPhasedToAlleleCodesMinus9(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t sample_ct, uint32_t phasepresent_ct, unsigned char* phasebytes, int32_t* allele_codes) {
  GenoarrPhasedToAlleleCodes(kGenoToIntcodeDPairs, genoarr, phasepresent, phaseinfo, sample_ct, phasepresent_ct, phasebytes, allele_codes);
}

// assumes transposed genoarr, phaseinfo
void GenoarrPhasedToHapCodes(const uintptr_t* genoarr, const uintptr_t* phaseinfo, uint32_t variant_batch_size, int32_t* hap0_codes_iter, int32_t* hap1_codes_iter);

void Dosage16ToFloatsMinus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, float* geno_float);

void Dosage16ToDoubles(const double* geno_double_pair_table, const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double);

// If all samples are missing, this errors out.
BoolErr Dosage16ToDoublesMeanimpute(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double);

// Currently requires trailing bits of genoarr to be zeroed out.
double LinearCombinationMeanimpute(const double* weights, const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct);

HEADER_INLINE void Dosage16ToDoublesMinus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double) {
  Dosage16ToDoubles(kGenoDoublePairs, genoarr, dosage_present, dosage_main, sample_ct, dosage_ct, geno_double);
}

// Does not zero out trailing bits of bitarr.
void BytesToBitsUnsafe(const uint8_t* boolbytes, uint32_t sample_ct, uintptr_t* bitarr);

// Bottom 2 bits are extracted from every byte.  Conveniently, -9 and 3 are
// treated identically.
// Does not zero out trailing bits of genoarr.
void BytesToGenoarrUnsafe(const int8_t* genobytes, uint32_t sample_ct, uintptr_t* genoarr);

// - If phasepresent_bytes is nullptr, phasepresent is not updated.  In this
//   case, phaseinfo is updated iff it's not nullptr.  It's okay for both
//   phasepresent and phaseinfo to be nullptr here.
// - Otherwise, phasepresent and phaseinfo are always updated; neither can be
//   nullptr.
void AlleleCodesToGenoarrUnsafe(const int32_t* allele_codes, const unsigned char* phasepresent_bytes, uint32_t sample_ct, uintptr_t* genoarr, uintptr_t* phasepresent, uintptr_t* phaseinfo);

void FloatsToDosage16(const float* floatarr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

void DoublesToDosage16(const double* doublearr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PGENLIB_FFI_SUPPORT_H__
