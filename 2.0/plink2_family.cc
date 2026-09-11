// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_family.h"

#include <assert.h>

#include "include/pgenlib_misc.h"
#include "include/pgenlib_write.h"
#include "include/plink2_bits.h"
#include "include/plink2_float.h"
#include "include/plink2_htable.h"
#include "include/plink2_simd.h"
#include "include/plink2_string.h"
#include "include/plink2_thread.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"
#include "plink2_compress_stream.h"
#include "include/plink2_stats.h"
#include "plink2_perm.h"
#include "plink2_random.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitMendel(MendelInfo* mendel_info_ptr) {
  mendel_info_ptr->flags = kfMendel0;
  mendel_info_ptr->max_trio_error = 1.0;
  mendel_info_ptr->max_var_error = 1.0;
  mendel_info_ptr->exclude_one_ratio = 0.0;
}

void PreinitFamilyInfo(FamilyInfo* fip) {
  fip->family_list = nullptr;
  fip->trio_list = nullptr;
  fip->trio_lookup = nullptr;
  fip->trio_fids = nullptr;
  fip->iids = nullptr;
  fip->sids = nullptr;
  fip->parent_to_trio_idxs = nullptr;
  fip->parent_to_trio_offsets = nullptr;
  fip->family_ct = 0;
  fip->trio_ct = 0;
  fip->duo_exists = 0;
  fip->ordered_erase_ok = 0;
  fip->max_fid_blen = 0;
  fip->max_iid_blen = 0;
  fip->max_sid_blen = 0;
}

// Main difference from the corresponding PLINK 1.9 function is that subsetted
// instead of raw sample-indexes are returned (since PgenReader takes care of
// subsetting).
// If trio_sample_include is non-null, it is initialized to the set of samples
// included in at least one trio(/duo), then *sample_ct_ptr is updated, and
// subsetted indexes are w.r.t. trio_sample_include.
PglErr GetTriosAndFamilies(const uintptr_t* orig_sample_include, const PedigreeIdInfo* piip, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, uint32_t raw_sample_ct, TrioFlags flags, uint32_t* sample_ct_ptr, uintptr_t* trio_sample_include, FamilyInfo* fip) {
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    uint32_t sample_ct = *sample_ct_ptr;
    const uint32_t include_duos = (flags / kfTrioDuos) & 1;
    const uint32_t raw_sample_ct_wduo = raw_sample_ct + include_duos;
    const uint32_t raw_sample_ct_wduol = BitCtToWordCt(raw_sample_ct_wduo);
    const uint32_t orig_sample_ct_wduo = sample_ct + include_duos;
    const uint32_t orig_sample_ct_wduol = BitCtToWordCt(orig_sample_ct_wduo);
    const char* sample_ids = piip->sii.sample_ids;
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    const uint32_t sample_id_htable_size = GetHtableFastSize(sample_ct);
    uintptr_t* founder_info2;
    uint64_t* trio_list_tmp;
    uint32_t* sample_id_htable;
    uint32_t* sample_include_cumulative_popcounts;
    char* idbuf;
    if (unlikely(bigstack_end_alloc_w(orig_sample_ct_wduol, &founder_info2) ||
                 bigstack_end_alloc_u64(sample_ct, &trio_list_tmp) ||
                 bigstack_end_alloc_u32(sample_id_htable_size, &sample_id_htable) ||
                 bigstack_end_alloc_u32(raw_sample_ct_wduol, &sample_include_cumulative_popcounts) ||
                 bigstack_end_alloc_c(max_sample_id_blen, &idbuf))) {
      goto GetTriosAndFamilies_ret_NOMEM;
    }
    {
      const uint32_t dup_sample_uidx = PopulateStrboxSubsetHtable(sample_ids, orig_sample_include, sample_ct, max_sample_id_blen, 0, sample_id_htable_size, sample_id_htable);
      if (unlikely(dup_sample_uidx)) {
        char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate FID+IID \"");
        const char* dup_sample_id = &(sample_ids[dup_sample_uidx * max_sample_id_blen]);
        const char* fid_end = AdvToDelim(dup_sample_id, '\t');
        write_iter = memcpyax(write_iter, dup_sample_id, fid_end - dup_sample_id, ' ');
        const char* iid_start = &(fid_end[1]);
        write_iter = strcpya(write_iter, iid_start);
        strcpy_k(write_iter, " in family analysis. (--select-sid-representatives may be useful.)\n");
        goto GetTriosAndFamilies_ret_INCONSISTENT_INPUT_WW;
      }
    }
    const uintptr_t* sample_include = orig_sample_include;
    const char* paternal_ids = piip->parental_id_info.paternal_ids;
    const char* maternal_ids = piip->parental_id_info.maternal_ids;
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (trio_sample_include) {
      ZeroWArr(raw_sample_ctl, trio_sample_include);
      // This is highly redundant with the next loop, but that's not a big
      // deal.
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = orig_sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(orig_sample_include, &sample_uidx_base, &cur_bits);
        if (IsSet(founder_info, sample_uidx)) {
          continue;
        }
        const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
        const char* iid_start = AdvPastDelim(cur_sample_id, '\t');
        const uintptr_t fid_blen = iid_start - cur_sample_id;
        memcpy(idbuf, cur_sample_id, fid_blen);
        const char* dad_iid = &(paternal_ids[sample_uidx * max_paternal_id_blen]);
        const char* mom_iid = &(maternal_ids[sample_uidx * max_maternal_id_blen]);
        uintptr_t iid_slen = strlen(dad_iid);
        uint32_t dad_uidx = UINT32_MAX;
        if (fid_blen + iid_slen < max_sample_id_blen) {
          memcpy(&(idbuf[fid_blen]), dad_iid, iid_slen + 1);
          dad_uidx = StrboxHtableFind(idbuf, sample_ids, sample_id_htable, max_sample_id_blen, fid_blen + iid_slen, sample_id_htable_size);
        }
        if (dad_uidx == UINT32_MAX) {
          if (!include_duos) {
            continue;
          }
          dad_uidx = raw_sample_ct;
        }
        // Don't need to repeat sanity checks that will always be performed in
        // the next loop.
        iid_slen = strlen(mom_iid);
        uint32_t mom_uidx = UINT32_MAX;
        if (fid_blen + iid_slen < max_sample_id_blen) {
          memcpy(&(idbuf[fid_blen]), mom_iid, iid_slen + 1);
          mom_uidx = StrboxHtableFind(idbuf, sample_ids, sample_id_htable, max_sample_id_blen, fid_blen + iid_slen, sample_id_htable_size);
        }
        if (mom_uidx == UINT32_MAX) {
          if ((!include_duos) || (dad_uidx == raw_sample_ct)) {
            continue;
          }
          SetBit(dad_uidx, trio_sample_include);
        } else {
          SetBit(mom_uidx, trio_sample_include);
          if (dad_uidx != raw_sample_ct) {
            SetBit(dad_uidx, trio_sample_include);
          }
        }
        SetBit(sample_uidx, trio_sample_include);
      }
      sample_include = trio_sample_include;
      sample_ct = PopcountWords(sample_include, raw_sample_ctl);
      *sample_ct_ptr = sample_ct;
      if (sample_ct == 0) {
        // Assumes PreinitFamilyInfo(fip) was called.
        goto GetTriosAndFamilies_ret_1;
      }
    }
    FillCumulativePopcounts(sample_include, raw_sample_ct_wduol, sample_include_cumulative_popcounts);
    CopyBitarrSubset(founder_info, sample_include, sample_ct, founder_info2);
    if (include_duos) {
      if (sample_ct % kBitsPerWord) {
        SetBit(sample_ct, founder_info2);
      } else {
        founder_info2[sample_ct / kBitsPerWord] = 1;
      }
    }
    uint32_t* child_to_trio_idxs;
    uint32_t* sample_child_cts;
    // over-allocate here, we shrink family_list later when we know how many
    // families there are
    uint64_t* family_list;
    uint64_t* family_htable;
    uint32_t* family_idxs;
    const uint32_t family_htable_size = GetHtableMinSize(sample_ct);
    if (unlikely(bigstack_alloc_u32(sample_ct, &child_to_trio_idxs) ||
                 bigstack_calloc_u32(sample_ct + 1, &sample_child_cts) ||
                 bigstack_alloc_u64(sample_ct, &family_list) ||
                 bigstack_calloc_u64(family_htable_size, &family_htable) ||
                 bigstack_alloc_u32(family_htable_size, &family_idxs))) {
      goto GetTriosAndFamilies_ret_NOMEM;
    }
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      child_to_trio_idxs[sample_idx] = UINT32_MAX;
    }
    const char* orig_sids = piip->sii.sids;
    const uintptr_t orig_max_sid_blen = piip->sii.max_sid_blen;
    uint64_t* trio_write_iter = trio_list_tmp;
    uintptr_t max_fid_blen = 2;
    uintptr_t max_iid_blen = 2;
    uintptr_t max_sid_blen = 2;
    uint32_t family_ct = 0;
    uint32_t duo_exists = 0;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* iid_start = AdvPastDelim(cur_sample_id, '\t');
      const uintptr_t fid_blen = iid_start - cur_sample_id;
      if (fid_blen > max_fid_blen) {
        max_fid_blen = fid_blen;
      }
      uintptr_t iid_slen = strlen(iid_start);
      if (iid_slen >= max_iid_blen) {
        max_iid_blen = iid_slen + 1;
      }
      if (orig_sids) {
        uintptr_t sid_slen = strlen(&(orig_sids[sample_uidx * orig_max_sid_blen]));
        if (sid_slen >= max_sid_blen) {
          max_sid_blen = sid_slen + 1;
        }
      }
      if (IsSet(founder_info2, sample_idx)) {
        continue;
      }
      memcpy(idbuf, cur_sample_id, fid_blen);
      const char* dad_iid = &(paternal_ids[sample_uidx * max_paternal_id_blen]);
      const char* mom_iid = &(maternal_ids[sample_uidx * max_maternal_id_blen]);
      iid_slen = strlen(dad_iid);
      uint32_t dad_uidx = UINT32_MAX;
      if (fid_blen + iid_slen < max_sample_id_blen) {
        memcpy(&(idbuf[fid_blen]), dad_iid, iid_slen + 1);
        dad_uidx = StrboxHtableFind(idbuf, sample_ids, sample_id_htable, max_sample_id_blen, fid_blen + iid_slen, sample_id_htable_size);
      }
      if (dad_uidx == UINT32_MAX) {
        if (!include_duos) {
          SetBit(sample_idx, founder_info2);
          continue;
        }
        dad_uidx = raw_sample_ct;
      } else {
        if (unlikely(dad_uidx == sample_uidx)) {
          idbuf[fid_blen - 1] = ' ';
          snprintf(g_logbuf, kLogbufSize, "Error: \"%s\" is their own parent.\n", idbuf);
          goto GetTriosAndFamilies_ret_INCONSISTENT_INPUT_WW;
        }
        if (unlikely((!IsSet(sex_nm, dad_uidx)) || (!IsSet(sex_male, dad_uidx)))) {
          idbuf[fid_blen - 1] = ' ';
          snprintf(g_logbuf, kLogbufSize, "Error: Father \"%s\" has unspecified or inconsistent sex.\n", idbuf);
          goto GetTriosAndFamilies_ret_INCONSISTENT_INPUT_WW;
        }
      }
      iid_slen = strlen(mom_iid);
      uint32_t mom_uidx = UINT32_MAX;
      if (fid_blen + iid_slen < max_sample_id_blen) {
        memcpy(&(idbuf[fid_blen]), mom_iid, iid_slen + 1);
        mom_uidx = StrboxHtableFind(idbuf, sample_ids, sample_id_htable, max_sample_id_blen, fid_blen + iid_slen, sample_id_htable_size);
      }
      if (mom_uidx == UINT32_MAX) {
        if ((!include_duos) || (dad_uidx == raw_sample_ct)) {
          SetBit(sample_idx, founder_info2);
          continue;
        }
        mom_uidx = raw_sample_ct;
        duo_exists = 1;
      } else {
        if (unlikely(mom_uidx == sample_uidx)) {
          idbuf[fid_blen - 1] = ' ';
          snprintf(g_logbuf, kLogbufSize, "Error: \"%s\" is their own parent.\n", idbuf);
          goto GetTriosAndFamilies_ret_INCONSISTENT_INPUT_WW;
        }
        if (unlikely((!IsSet(sex_nm, mom_uidx)) || IsSet(sex_male, mom_uidx))) {
          idbuf[fid_blen - 1] = ' ';
          snprintf(g_logbuf, kLogbufSize, "Error: Mother \"%s\" has unspecified or inconsistent sex.\n", idbuf);
          goto GetTriosAndFamilies_ret_INCONSISTENT_INPUT_WW;
        }
      }
      const uintptr_t mom_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, mom_uidx);
      const uintptr_t dad_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, dad_uidx);
      if (mom_idx != sample_ct) {
        sample_child_cts[mom_idx] += 1;
      }
      if (dad_idx != sample_ct) {
        sample_child_cts[dad_idx] += 1;
      } else {
        duo_exists = 1;
      }
      const uint64_t family_code = (S_CAST(uint64_t, mom_idx) << 32) | dad_idx;
      // Simple linear-probing setup.  Could make this fancier if adversarial
      // input is ever a concern.
      uint32_t initial_hashval;
      {
        uint64_t ullii = dad_idx + 2LLU * mom_idx;
        assert(ullii < 2LLU * family_htable_size);
        if (ullii >= family_htable_size) {
          ullii -= family_htable_size;
        }
        initial_hashval = ullii;
      }
      uint32_t family_idx;
      for (uint32_t hashval = initial_hashval; ; ) {
        uint64_t ullii = family_htable[hashval];
        if (ullii == 0) {
          family_idx = family_ct++;
          family_list[family_idx] = family_code;
          family_htable[hashval] = family_code;
          family_idxs[hashval] = family_idx;
          break;
        } else if (ullii == family_code) {
          family_idx = family_idxs[hashval];
          break;
        }
        ++hashval;
        if (hashval == family_htable_size) {
          hashval = 0;
        }
      }
      // Save sample_uidxs first, then convert to sample_idxs after
      // populate_ids step done.
      *trio_write_iter++ = (S_CAST(uint64_t, family_idx) << 32) | sample_uidx;
    }
    const uintptr_t trio_ct = trio_write_iter - trio_list_tmp;
    BigstackReset(family_list);
    BigstackFinalizeU64(family_list, trio_ct);
    BigstackEndReset(trio_list_tmp);
    uint64_t* trio_write;
    if (unlikely(bigstack_alloc_u64(trio_ct, &trio_write))) {
      goto GetTriosAndFamilies_ret_NOMEM;
    }
    memcpy(trio_write, trio_list_tmp, trio_ct * sizeof(int64_t));
    STD_SORT(trio_ct, u64cmp, trio_write);
    BigstackEndReset(founder_info2);
    uint32_t* parent_to_trio_offsets = sample_child_cts;
    {
      uint32_t cumsum = 0;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint32_t cur_child_ct = parent_to_trio_offsets[sample_idx];
        parent_to_trio_offsets[sample_idx] = cumsum;
        cumsum += cur_child_ct;
      }
      parent_to_trio_offsets[sample_ct] = cumsum;
    }
    uint32_t* parent_to_trio_idxs;
    // theoretically possible for trio_lookup[] to have size > 2^32
    uint32_t* trio_lookup;
    if (unlikely(bigstack_alloc_u32(parent_to_trio_offsets[sample_ct], &parent_to_trio_idxs) ||
                 bigstack_alloc_u32(trio_ct * 3, &trio_lookup))) {
      goto GetTriosAndFamilies_ret_NOMEM;
    }
    const uint32_t populate_ids = (flags / kfTrioPopulateIds) & 1;
    char* trio_fids = nullptr;
    char* iids = nullptr;
    char* sids = nullptr;
    if (populate_ids) {
      if (unlikely(bigstack_alloc_c(trio_ct * max_fid_blen, &trio_fids) ||
                   bigstack_alloc_c((sample_ct + duo_exists) * max_iid_blen, &iids))) {
        goto GetTriosAndFamilies_ret_NOMEM;
      }
      if (flags & kfTrioPopulateSids) {
        if (unlikely(bigstack_alloc_c((sample_ct + duo_exists) * max_sid_blen, &sids))) {
          goto GetTriosAndFamilies_ret_NOMEM;
        }
      }
    }
    if (populate_ids) {
      sample_uidx_base = 0;
      cur_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
        const char* iid_start = AdvPastDelim(cur_sample_id, '\t');
        strcpy(&(iids[sample_idx * max_iid_blen]), iid_start);
        if (sids) {
          char* dst = &(sids[sample_idx * max_sid_blen]);
          if (!orig_sids) {
            strcpy_k(dst, "0");
          } else {
            strcpy(dst, &(orig_sids[sample_uidx * orig_max_sid_blen]));
          }
        }
      }
      if (duo_exists) {
        strcpy_k(&(iids[sample_ct * max_iid_blen]), "0");
        if (sids) {
          strcpy_k(&(sids[sample_ct * max_sid_blen]), "0");
        }
      }
    }
    uint32_t* trio_lookup_write_iter = trio_lookup;
    uint32_t ordered_erase_ok = 1;
    for (uintptr_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
      const uint64_t trio_code = trio_write[trio_idx];
      const uint32_t sample_uidx = S_CAST(uint32_t, trio_code);
      const uintptr_t family_idx = trio_code >> 32;
      if (populate_ids) {
        const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
        const char* fid_end = AdvToDelim(cur_sample_id, '\t');
        const uint32_t fid_slen = fid_end - cur_sample_id;
        memcpyx(&(trio_fids[trio_idx * max_fid_blen]), cur_sample_id, fid_slen, '\0');
      }
      const uintptr_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_uidx);
      child_to_trio_idxs[sample_idx] = trio_idx;
      trio_write[trio_idx] = (S_CAST(uint64_t, family_idx) << 32) | sample_idx;
      const uint64_t family_code = family_list[family_idx];
      const uint32_t dad_idx = S_CAST(uint32_t, family_code);
      const uint32_t mom_idx = family_code >> 32;
      if (ordered_erase_ok) {
        ordered_erase_ok = (child_to_trio_idxs[dad_idx] == UINT32_MAX) && (child_to_trio_idxs[mom_idx] == UINT32_MAX);
      }
      if (dad_idx != sample_ct) {
        parent_to_trio_idxs[parent_to_trio_offsets[dad_idx]++] = trio_idx;
      }
      if (mom_idx != sample_ct) {
        parent_to_trio_idxs[parent_to_trio_offsets[mom_idx]++] = trio_idx;
      }
      *trio_lookup_write_iter++ = sample_idx;
      *trio_lookup_write_iter++ = dad_idx;
      *trio_lookup_write_iter++ = mom_idx;
    }
    for (uint32_t sample_idx = sample_ct; sample_idx; --sample_idx) {
      parent_to_trio_offsets[sample_idx] = parent_to_trio_offsets[sample_idx - 1];
    }
    parent_to_trio_offsets[0] = 0;
    fip->family_list = family_list;
    fip->trio_list = trio_write;
    fip->trio_lookup = trio_lookup;
    fip->trio_fids = trio_fids;
    fip->iids = iids;
    fip->sids = sids;
    fip->parent_to_trio_idxs = parent_to_trio_idxs;
    fip->parent_to_trio_offsets = parent_to_trio_offsets;
    fip->child_to_trio_idxs = child_to_trio_idxs;
    fip->family_ct = family_ct;
    fip->trio_ct = trio_ct;
    fip->duo_exists = duo_exists;
    fip->ordered_erase_ok = ordered_erase_ok;
    fip->max_fid_blen = max_fid_blen;
    fip->max_iid_blen = max_iid_blen;
    fip->max_sid_blen = max_sid_blen;
  }
  while (0) {
  GetTriosAndFamilies_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GetTriosAndFamilies_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 GetTriosAndFamilies_ret_1:
  BigstackEndReset(bigstack_end_mark);
  return reterr;
}

void InitMultiallelicWideCodes(const uint16_t* hc_to_allele_codes, const uintptr_t* sex_male_collapsed, const uintptr_t* sex_female_collapsed, const uintptr_t* genoarr, const uintptr_t* patch_01_set, const AlleleCode* patch_01_vals, const uintptr_t* patch_10_set, const AlleleCode* patch_10_vals, uint32_t sample_ct, uint32_t male_ct, uint32_t female_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, uint32_t is_x, uint32_t is_y, AlleleCode* wide_codes) {
  // Similar to PglMultiallelicSparseToDenseMiss().
  GenoarrLookup256x2bx4(genoarr, hc_to_allele_codes, sample_ct, wide_codes);
  if (!is_y) {
    if (patch_01_ct) {
      uintptr_t sample_idx_base = 0;
      uintptr_t cur_bits = patch_01_set[0];
      AlleleCode* wide_codes1 = &(wide_codes[1]);
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        wide_codes1[2 * sample_idx] = patch_01_vals[uii];
      }
    }
  }
  if (patch_10_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = patch_10_set[0];
    if (!is_y) {
      const DoubleAlleleCode* patch_10_vals_alias = R_CAST(const DoubleAlleleCode*, patch_10_vals);
      DoubleAlleleCode* __attribute__((may_alias)) wide_codes_alias = R_CAST(DoubleAlleleCode*, wide_codes);
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
        wide_codes_alias[sample_idx] = patch_10_vals_alias[uii];
      }
    } else {
      const AlleleCode* patch_10_vals_stop = &(patch_10_vals[2 * patch_10_ct]);
      for (const AlleleCode* patch_10_vals_iter = patch_10_vals; patch_10_vals_iter != patch_10_vals_stop; ) {
        const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
        AlleleCode ac0 = *patch_10_vals_iter++;
        AlleleCode ac1 = *patch_10_vals_iter++;
        if (ac0 != ac1) {
          ac0 = kMissingAlleleCode;
          ac1 = kMissingAlleleCode;
        }
        wide_codes[2 * sample_idx] = ac0;
        wide_codes[2 * sample_idx + 1] = ac1;
      }
    }
  }
  // could vectorize these loops with logic similar to InverseMovemaskFF?
  if (is_x) {
    // male hets -> missing
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = sex_male_collapsed[0];
    for (uint32_t male_idx = 0; male_idx != male_ct; ++male_idx) {
      const uintptr_t sample_idx = BitIter1(sex_male_collapsed, &sample_idx_base, &cur_bits);
      if (wide_codes[2 * sample_idx] != wide_codes[2 * sample_idx + 1]) {
        wide_codes[2 * sample_idx] = kMissingAlleleCode;
        wide_codes[2 * sample_idx + 1] = kMissingAlleleCode;
      }
    }
  } else if (is_y) {
    // female -> missing
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = sex_female_collapsed[0];
    DoubleAlleleCode* __attribute__((may_alias)) wide_codes_alias = R_CAST(DoubleAlleleCode*, wide_codes);
    for (uint32_t female_idx = 0; female_idx != female_ct; ++female_idx) {
      const uintptr_t sample_idx = BitIter1(sex_female_collapsed, &sample_idx_base, &cur_bits);
      wide_codes_alias[sample_idx] = kMissingDoubleAlleleCode;
    }
  }
}

typedef struct MendelErrorScanCtxStruct {
  const FamilyInfo* fip;
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  // On chrX, male het -> missing.
  const uintptr_t* sex_male_collapsed;
  const uintptr_t* sex_male_collapsed_interleaved;
  // On chrY, female or het -> missing.
  const uintptr_t* sex_female_collapsed;
  const uintptr_t* sex_female_collapsed_interleaved;

  MendelFlags flags;
  uint32_t sample_ct;
  uint32_t male_ct;
  uint32_t female_ct;
  uint32_t max_difflist_len;
  double max_var_error;

  PgenReader** pgr_ptrs;
  // in duo_exists case, each genovec must have space for trailing sample
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_size;

  uintptr_t** thread_trio_include_bufs;
  uint32_t** thread_trio_error_ct_acc_bufs;
  uint32_t** thread_trio_missing_ct_acc_bufs;
  // Only relevant in var_first case, where we don't proceed to add to
  // trio_{error,missing}_ct_acc when the variant is filtered out.
  uint32_t** thread_trio_error_ct_cur_bufs;
  uint32_t** thread_trio_missing_ct_cur_bufs;
  // Only need this for multiallelic variants.
  AlleleCode** thread_wide_code_bufs;

  uint32_t* variant_error_cts[2];
  uint32_t* variant_obs_cts[2];
  // One bitarray (trio_ctl words) per variant, if needed.
  uintptr_t* variant_errors[2];
  // Variant-major, trio-minor.  Each entry takes 1, 6, or 7 AlleleCodes,
  // depending on CODE/ERROR column output settings.
  //   [0]: error code
  //   [1-2]: paternal genotype to render
  //            kMissingAlleleCode/0 = <blank> (chrY female)
  //            kMissingAlleleCode/(kMissingAlleleCode-1) = *
  //            kMissingAlleleCode/kMissingAlleleCode = */*
  //            <non-missing>/kMissingAlleleCode = <haploid>
  //   [3-4]: maternal genotype to render
  //   [5-6]: child genotype to render
  AlleleCode* edescrips[2];

  // length (3 * trio_ct): child, paternal, maternal
  uint32_t** thread_trio_error_cts;
  uint32_t** thread_trio_missing_cts;

  uint64_t err_info;
} MendelErrorScanCtx;

const uint16_t kHetMissingToAlleleCodes[1024] = QUAD_TABLE256(0, 0xffff, 0x101, 0xffff);

// index bits:
//   0-1: dad_geno
//   2-3: mom_geno
//   4-5: child_geno
// result:
//   bit 0 set if error implicates child
//   bit 8 set if error implicates dad
//   bit 16 set if error implicates mom
//   error code in bits 24+
const uint32_t kBiallelicMendelErrorTableAutosomalOrX[48] = {
  0, 0, 0x6000101, 0,
  0, 0, 0x6000101, 0,
  0x7010001, 0x7010001, 0x8000001, 0x7010001,
  0, 0, 0x6000101, 0,
  0x2010101, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0x1010101, 0,
  0, 0, 0, 0,
  0x5000001, 0x4010001, 0x4010001, 0x4010001,
  0x3000101, 0, 0, 0,
  0x3000101, 0, 0, 0,
  0x3000101, 0, 0, 0};

const uint32_t kBiallelicMendelErrorTableChrY[48] = {
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0xb000101, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0xc000101, 0, 0, 0};

const uint32_t kBiallelicMendelErrorTableChrM[48] = {
  0, 0, 0, 0,
  0, 0, 0, 0,
  0x9010001, 0x9010001, 0x9010001, 0x9010001,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0xa010001, 0xa010001, 0xa010001, 0xa010001,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0};

// index bits:
//   0: dad has (or at least may have) child_ac0
//   1: dad has child_ac1
//   2: mom has child_ac0
//   3: mom has child_ac1
//   4-5: (child_ac0 == 0) + (child_ac1 == 0)
// (impossible entries marked with UINT32_MAX)
const uint32_t kMultiallelicMendelErrorTableAutosomalOrX[48] = {
  0x5000001, 0x4010001, 0x4010001, 0x4010001,
  0x3000101, 0x1010101, 0, 0,
  0x3000101, 0, 0x1010101, 0,
  0x3000101, 0, 0, 0,
  0x5000001, 0x4010001, 0x4010001, 0x4010001,
  0x3000101, 0x2010101, 0, 0,
  0x3000101, 0, 0x1010101, 0,
  0x3000101, 0, 0, 0,
  0x8000001, UINT32_MAX, UINT32_MAX, 0x7010001,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  0x6000101, UINT32_MAX, UINT32_MAX, 0};

const uint32_t kMultiallelicMendelErrorTableChrY[48] = {
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  0xc000101, UINT32_MAX, UINT32_MAX, 0,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  0xb000101, UINT32_MAX, UINT32_MAX, 0};

// MendelEraseErrors() doesn't force d0=d1=1, and only looks at the first 16
// entries.
const uint32_t kMultiallelicMendelErrorTableChrM[48] = {
  0xa010001, 0xa010001, 0xa010001, 0xa010001,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, 0xa010001,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, 0,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, 0,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, 0,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, 0x9010001,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX,
  UINT32_MAX, UINT32_MAX, UINT32_MAX, 0};

const DoubleAlleleCode kGenoToDac[4] = {0, 1 << (8 * sizeof(AlleleCode)), 1 + (1 << 8 * sizeof(AlleleCode)), kMissingDoubleAlleleCode};

THREAD_FUNC_DECL MendelErrorScanThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  MendelErrorScanCtx* ctx = S_CAST(MendelErrorScanCtx*, arg->sharedp->context);

  const FamilyInfo* fip = ctx->fip;
  const uint32_t* trio_lookup = fip->trio_lookup;
  const uint32_t* parent_to_trio_idxs = fip->parent_to_trio_idxs;
  const uint32_t* parent_to_trio_offsets = fip->parent_to_trio_offsets;
  const uint32_t* child_to_trio_idxs = fip->child_to_trio_idxs;
  const uintptr_t trio_ct = fip->trio_ct;
  const uintptr_t trio_ctl = BitCtToWordCt(trio_ct);
  const uint32_t trio_vec_ct = Int32CtToVecCt(trio_ct);
  const uint32_t trio_ctav = trio_vec_ct * kInt32PerVec;
  const uint32_t duo_exists = fip->duo_exists;
  const uintptr_t* variant_include = ctx->variant_include;
  const ChrInfo* cip = ctx->cip;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t* sample_include = ctx->sample_include;
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  const uintptr_t* sex_male_collapsed = ctx->sex_male_collapsed;
  const uintptr_t* sex_female_collapsed = ctx->sex_female_collapsed;
  const uintptr_t* sex_male_collapsed_interleaved = ctx->sex_male_collapsed_interleaved;
  const uintptr_t* sex_female_collapsed_interleaved = ctx->sex_female_collapsed_interleaved;
  const uint32_t exclude_missing_from_denom = !(ctx->flags & kfMendelMissingInDenom);
  const uint32_t var_first = (ctx->flags / kfMendelFilterVarFirst) & 1;
  const uint32_t ecode_col = (ctx->flags / kfMendelRptColCode) & 1;
  const uint32_t edescrip_col = (ctx->flags / kfMendelRptColError) & 1;
  const uintptr_t edescrip_stride = trio_ct * (ecode_col + 6 * edescrip_col);
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctv2 = NypCtToVecCt(sample_ct);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t male_ct = ctx->male_ct;
  const uint32_t female_ct = ctx->female_ct;
  const uint32_t max_difflist_len = ctx->max_difflist_len;
  const double max_var_error = ctx->max_var_error;
  uintptr_t* genovec = ctx->genovecs[tidx];
  PgenVariant pgv;
  pgv.genovec = genovec;
  SetPgvThreadMhcNull(sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  uintptr_t* raregeno = ctx->raregenos[tidx];
  uint32_t* difflist_sample_ids = ctx->difflist_sample_id_bufs[tidx];

  uintptr_t* trio_include_buf = ctx->thread_trio_include_bufs[tidx];
  uint32_t* trio_error_ct_acc_buf = ctx->thread_trio_error_ct_acc_bufs[tidx];
  uint32_t* trio_missing_ct_acc_buf = ctx->thread_trio_missing_ct_acc_bufs[tidx];
  ZeroU32Arr(trio_ctav, trio_error_ct_acc_buf);
  ZeroU32Arr(trio_ctav, trio_missing_ct_acc_buf);
  uint32_t* trio_error_ct_cur_buf = trio_error_ct_acc_buf;
  uint32_t* trio_missing_ct_cur_buf = trio_missing_ct_acc_buf;
  if (var_first) {
    trio_error_ct_cur_buf = ctx->thread_trio_error_ct_cur_bufs[tidx];
    trio_missing_ct_cur_buf = ctx->thread_trio_missing_ct_cur_bufs[tidx];
  }
  AlleleCode* wide_codes = nullptr;
  if (allele_idx_offsets) {
    wide_codes = ctx->thread_wide_code_bufs[tidx];
  }
  uint32_t* trio_error_cts = ctx->thread_trio_error_cts[tidx];
  uint32_t* trio_missing_cts = ctx->thread_trio_missing_cts[tidx];
  ZeroU32Arr(3 * trio_ct, trio_error_cts);
  ZeroU32Arr(3 * trio_ct, trio_missing_cts);

  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);

  uint32_t cur_allele_ct = 2;
  uint64_t new_err_info = 0;
  uint32_t parity = 0;
  uint32_t error_ct_fill = 0;
  uint32_t missing_ct_fill = 0;
  do {
    const uint32_t cur_block_size = ctx->cur_block_size;
    const uint32_t cur_variant_bidx_offset = (tidx * cur_block_size) / calc_thread_ct;
    uint32_t* variant_error_cts = &(ctx->variant_error_cts[parity][cur_variant_bidx_offset]);
    uint32_t* variant_obs_cts = &(ctx->variant_obs_cts[parity][cur_variant_bidx_offset]);
    uintptr_t* cur_variant_errors = nullptr;
    AlleleCode* cur_edescrips = nullptr;
    if (ctx->variant_errors[parity]) {
      cur_variant_errors = &(ctx->variant_errors[parity][cur_variant_bidx_offset * trio_ctl]);
      if (ctx->edescrips[parity]) {
        cur_edescrips = &(ctx->edescrips[parity][cur_variant_bidx_offset * edescrip_stride]);
      }
    }
    const uint32_t cur_variant_bidx_ct = (((tidx + 1) * cur_block_size) / calc_thread_ct) - cur_variant_bidx_offset;
    uintptr_t variant_uidx_base;
    uintptr_t cur_variant_include_bits;
    BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &cur_variant_include_bits);
    const uint16_t* hc_to_allele_codes = nullptr;
    const uint32_t* biallelic_mendel_error_table = nullptr;
    const uint32_t* multiallelic_mendel_error_table = nullptr;
    uint32_t chr_end = 0;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t is_mt = 0;
    uint32_t is_xymt = 0;
    for (uint32_t cur_variant_bidx = 0; cur_variant_bidx != cur_variant_bidx_ct; ++cur_variant_bidx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_variant_include_bits);
      if (variant_uidx >= chr_end) {
        const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        is_x = (chr_idx == x_code);
        is_y = (chr_idx == y_code);
        is_mt = (chr_idx == mt_code);
        is_xymt = is_x || is_y || is_mt;
        hc_to_allele_codes = kHcToAlleleCodes;
        biallelic_mendel_error_table = kBiallelicMendelErrorTableAutosomalOrX;
        multiallelic_mendel_error_table = kMultiallelicMendelErrorTableAutosomalOrX;
        if (is_y) {
          hc_to_allele_codes = kHetMissingToAlleleCodes;
          biallelic_mendel_error_table = kBiallelicMendelErrorTableChrY;
          multiallelic_mendel_error_table = kMultiallelicMendelErrorTableChrY;
        } else if (is_mt) {
          biallelic_mendel_error_table = kBiallelicMendelErrorTableChrM;
          multiallelic_mendel_error_table = kMultiallelicMendelErrorTableChrM;
        }
      }
      if (allele_idx_offsets) {
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
      }
      if (cur_variant_errors) {
        ZeroWArr(trio_ctl, cur_variant_errors);
      }
      uint32_t variant_error_ct = 0;
      uint32_t variant_error_ct_denom = trio_ct;
      if (var_first) {
        ZeroU32Arr(trio_ctav, trio_error_ct_cur_buf);
        ZeroU32Arr(trio_ctav, trio_missing_ct_cur_buf);
      }
      AlleleCode* edescrips_write_iter = cur_edescrips;
      if (cur_allele_ct == 2) {
        uint32_t difflist_common_geno;
        uint32_t difflist_len;
        const PglErr reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_difflist_len, variant_uidx, pgrp, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto MendelErrorScanThread_err;
        }
        if (difflist_common_geno != UINT32_MAX) {
          PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, sample_ct, difflist_len, genovec);
        }
        if (duo_exists) {
          SetNyparrEntryTo3(sample_ct, genovec);
        }
        if (((difflist_common_geno == 0) || (difflist_common_geno == 2)) && (!is_xymt)) {
          if (difflist_len > 0) {
            ZeroWArr(trio_ctl, trio_include_buf);
            const uint32_t* difflist_sample_ids_iter = difflist_sample_ids;
            const uint32_t word_ct_m1 = (difflist_len - 1) / kBitsPerWordD2;
            uint32_t loop_len = kBitsPerWordD2;
            for (uint32_t widx = 0; ; ++widx) {
              if (widx >= word_ct_m1) {
                if (widx > word_ct_m1) {
                  break;
                }
                loop_len = ModNz(difflist_len, kBitsPerWordD2);
              }
              uintptr_t raregeno_word = raregeno[widx];
              for (uint32_t difflist_idx_lowbits = 0; difflist_idx_lowbits != loop_len; ++difflist_idx_lowbits) {
                const uint32_t sample_idx = *difflist_sample_ids_iter++;
                const uint32_t cur_geno = raregeno_word & 3;
                SetBit(child_to_trio_idxs[sample_idx], trio_include_buf);
                if (cur_geno != 1) {
                  const uint32_t offset_start = parent_to_trio_offsets[sample_idx];
                  const uint32_t offset_end = parent_to_trio_offsets[sample_idx + 1];
                  for (uint32_t uii = offset_start; uii != offset_end; ++uii) {
                    SetBit(parent_to_trio_idxs[uii], trio_include_buf);
                  }
                }
                raregeno_word = raregeno_word >> 2;
              }
            }
            uintptr_t trio_idx_offset = 0;
            for (uint32_t widx = 0; widx != trio_ctl; ++widx, trio_idx_offset += kBitsPerWord) {
              for (uintptr_t trio_include_bits = trio_include_buf[widx]; trio_include_bits; trio_include_bits &= trio_include_bits - 1) {
                const uintptr_t trio_idx = trio_idx_offset + ctzw(trio_include_bits);
                const uint32_t* cur_trio_lookup = &(trio_lookup[trio_idx * 3]);
                const uint32_t child_idx = cur_trio_lookup[0];
                const uint32_t child_geno = GetNyparrEntry(genovec, child_idx);
                if (child_geno == 3) {
                  variant_error_ct_denom -= exclude_missing_from_denom;
                  trio_missing_ct_cur_buf[trio_idx] += 0x10101;
                  continue;
                }
                const uint32_t dad_idx = cur_trio_lookup[1];
                const uint32_t mom_idx = cur_trio_lookup[2];
                const uint32_t dad_geno = GetNyparrEntry(genovec, dad_idx);
                const uint32_t mom_geno = GetNyparrEntry(genovec, mom_idx);
                if (dad_geno == 3) {
                  if (mom_geno == 3) {
                    variant_error_ct_denom -= exclude_missing_from_denom;
                    trio_missing_ct_cur_buf[trio_idx] += 0x10101;
                    continue;
                  }
                  trio_missing_ct_cur_buf[trio_idx] += 0x100;
                } else if (mom_geno == 3) {
                  trio_missing_ct_cur_buf[trio_idx] += 0x10000;
                }
                const uint32_t error_result = biallelic_mendel_error_table[dad_geno + mom_geno * 4 + child_geno * 16];
                if (!error_result) {
                  continue;
                }
                trio_error_ct_cur_buf[trio_idx] += error_result & 0xffffff;
                ++variant_error_ct;
                if (!cur_variant_errors) {
                  continue;
                }
                SetBit(trio_idx, cur_variant_errors);
                if (!cur_edescrips) {
                  continue;
                }
                AlleleCode mendel_ecode = error_result >> 24;
                if (ecode_col) {
                  *edescrips_write_iter++ = mendel_ecode;
                }
                if (!edescrip_col) {
                  continue;
                }
                DoubleAlleleCode dac = kGenoToDac[dad_geno];
                memcpy_k(edescrips_write_iter, &dac, 2 * sizeof(AlleleCode));
                dac = kGenoToDac[mom_geno];
                memcpy_k(&(edescrips_write_iter[2]), &dac, 2 * sizeof(AlleleCode));
                dac = kGenoToDac[child_geno];
                memcpy_k(&(edescrips_write_iter[4]), &dac, 2 * sizeof(AlleleCode));
                switch (mendel_ecode) {
                case 3:
                case 6:
                  edescrips_write_iter[2] = kMissingAlleleCode;
                  edescrips_write_iter[3] = kMissingAlleleCode;
                  break;
                case 4:
                case 7:
                  edescrips_write_iter[0] = kMissingAlleleCode;
                  edescrips_write_iter[1] = kMissingAlleleCode;
                  break;
                }
                edescrips_write_iter = &(edescrips_write_iter[6]);
              }
            }
          }
        } else {
          // chrX: male hets -> missing
          // chrY: female -> missing, het -> missing
          if (is_x) {
            SetMaleHetMissing(sex_male_collapsed_interleaved, sample_ctv2, genovec);
          } else if (is_y) {
            InterleavedSetMissing(sex_female_collapsed_interleaved, sample_ctv2, genovec);
            SetHetMissing(sample_ctl2, genovec);
          }
          // possible todo: try transposing and processing in batches of
          // kBitsPerWordD2 variants
          const uint32_t* trio_lookup_iter = trio_lookup;
          uint32_t ignore_father = is_mt;
          for (uintptr_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
            const uint32_t child_idx = trio_lookup_iter[0];
            trio_lookup_iter = &(trio_lookup_iter[3]);
            const uint32_t child_geno = GetNyparrEntry(genovec, child_idx);
            if (child_geno == 3) {
              variant_error_ct_denom -= exclude_missing_from_denom;
              trio_missing_ct_cur_buf[trio_idx] += 0x10101;
              continue;
            }
            const uint32_t dad_idx = trio_lookup_iter[-2];
            const uint32_t mom_idx = trio_lookup_iter[-1];
            uint32_t dad_geno = GetNyparrEntry(genovec, dad_idx);
            const uint32_t mom_geno = GetNyparrEntry(genovec, mom_idx);
            if (is_x) {
              ignore_father = IsSet(sex_male_collapsed, child_idx);
            }
            if (ignore_father || (dad_geno == 3)) {
              if (mom_geno == 3) {
                variant_error_ct_denom -= exclude_missing_from_denom;
                trio_missing_ct_cur_buf[trio_idx] += 0x10101;
                continue;
              }
              dad_geno = 3;
              trio_missing_ct_cur_buf[trio_idx] += 0x100;
            } else if (mom_geno == 3) {
              trio_missing_ct_cur_buf[trio_idx] += 0x10000;
            }
            const uint32_t error_result = biallelic_mendel_error_table[dad_geno + mom_geno * 4 + child_geno * 16];
            if (!error_result) {
              continue;
            }
            trio_error_ct_cur_buf[trio_idx] += error_result & 0xffffff;
            ++variant_error_ct;
            if (!cur_variant_errors) {
              continue;
            }
            SetBit(trio_idx, cur_variant_errors);
            if (!cur_edescrips) {
              continue;
            }
            AlleleCode mendel_ecode = error_result >> 24;
            if (ignore_father) {
              if (mendel_ecode == 7) {
                mendel_ecode = 9;
              } else if (mendel_ecode == 4) {
                mendel_ecode = 10;
              }
            }
            if (ecode_col) {
              *edescrips_write_iter++ = mendel_ecode;
            }
            if (!edescrip_col) {
              continue;
            }
            DoubleAlleleCode dac = kGenoToDac[dad_geno];
            memcpy_k(edescrips_write_iter, &dac, 2 * sizeof(AlleleCode));
            dac = kGenoToDac[mom_geno];
            memcpy_k(&(edescrips_write_iter[2]), &dac, 2 * sizeof(AlleleCode));
            dac = kGenoToDac[child_geno];
            memcpy_k(&(edescrips_write_iter[4]), &dac, 2 * sizeof(AlleleCode));
            if (is_x) {
              // render dad as haploid
              edescrips_write_iter[1] = kMissingAlleleCode - (dad_geno == 3);
            }
            switch (mendel_ecode) {
            case 3:
            case 6:
              edescrips_write_iter[2] = kMissingAlleleCode;
              edescrips_write_iter[3] = kMissingAlleleCode;
              break;
            case 4:
            case 7:
              edescrips_write_iter[0] = kMissingAlleleCode;
              edescrips_write_iter[1] = kMissingAlleleCode;
              break;
            case 9:
            case 10:
              // chrX male, or any chrM: render dad as '*', child as haploid
              // unless heterozygous.
              edescrips_write_iter[0] = kMissingAlleleCode;
              edescrips_write_iter[1] = kMissingAlleleCode - 1;
              // mom is always rendered as diploid on chrX; try to render as
              // haploid on chrM.
              if (is_mt && (mom_geno != 1)) {
                edescrips_write_iter[3] = kMissingAlleleCode;
              }
              if (child_geno != 1) {
                edescrips_write_iter[5] = kMissingAlleleCode;
              }
              break;
            case 11:
            case 12:
              // chrY non-female: render dad and child as haploid, don't render
              // mom at all (don't even include preceding "x")
              edescrips_write_iter[1] = kMissingAlleleCode;
              edescrips_write_iter[2] = kMissingAlleleCode;
              edescrips_write_iter[3] = 0;
              edescrips_write_iter[5] = kMissingAlleleCode;
              break;
            }
            edescrips_write_iter = &(edescrips_write_iter[6]);
          }
        }
      } else {
        const PglErr reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, pgrp, &pgv);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto MendelErrorScanThread_err;
        }
        InitMultiallelicWideCodes(hc_to_allele_codes, sex_male_collapsed, sex_female_collapsed, genovec, pgv.patch_01_set, pgv.patch_01_vals, pgv.patch_10_set, pgv.patch_10_vals, sample_ct, male_ct, female_ct, pgv.patch_01_ct, pgv.patch_10_ct, is_x, is_y, wide_codes);
        if (duo_exists) {
          wide_codes[2 * sample_ct] = kMissingAlleleCode;
          wide_codes[2 * sample_ct + 1] = kMissingAlleleCode;
        }
        const uint32_t* trio_lookup_iter = trio_lookup;
        uint32_t ignore_father = is_mt;  // chrM, or chrX male
        for (uintptr_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
          const uint32_t child_idx = trio_lookup_iter[0];
          trio_lookup_iter = &(trio_lookup_iter[3]);
          const AlleleCode child_ac0 = wide_codes[child_idx * 2];
          if (child_ac0 == kMissingAlleleCode) {
            variant_error_ct_denom -= exclude_missing_from_denom;
            trio_missing_ct_cur_buf[trio_idx] += 0x10101;
            continue;
          }
          const AlleleCode child_ac1 = wide_codes[child_idx * 2 + 1];
          const uint32_t dad_idx = trio_lookup_iter[-2];
          const uint32_t mom_idx = trio_lookup_iter[-1];
          const AlleleCode dad_ac0 = wide_codes[dad_idx * 2];
          const AlleleCode dad_ac1 = wide_codes[dad_idx * 2 + 1];
          uint32_t d0 = (dad_ac0 == child_ac0) || (dad_ac1 == child_ac0);
          uint32_t d1 = (dad_ac0 == child_ac1) || (dad_ac1 == child_ac1);
          if (is_x) {
            ignore_father = IsSet(sex_male_collapsed, child_idx);
          }
          const AlleleCode mom_ac0 = wide_codes[mom_idx * 2];
          if (ignore_father || (dad_ac0 == kMissingAlleleCode)) {
            if (mom_ac0 == kMissingAlleleCode) {
              variant_error_ct_denom -= exclude_missing_from_denom;
              trio_missing_ct_cur_buf[trio_idx] += 0x10101;
              continue;
            }
            d0 = 1;
            d1 = 1;
            trio_missing_ct_cur_buf[trio_idx] += 0x100;
          } else if (mom_ac0 == kMissingAlleleCode) {
            trio_missing_ct_cur_buf[trio_idx] += 0x10000;
          }
          const AlleleCode mom_ac1 = wide_codes[mom_idx * 2 + 1];
          const uint32_t m0 = (mom_ac0 == child_ac0) || (mom_ac1 == child_ac0) || (mom_ac0 == kMissingAlleleCode);
          const uint32_t m1 = (mom_ac0 == child_ac1) || (mom_ac1 == child_ac1) || (mom_ac0 == kMissingAlleleCode);
          const uint32_t child_ref_ct = (child_ac0 == 0) + (child_ac1 == 0);
          const uint32_t error_result = multiallelic_mendel_error_table[d0 + d1 * 2 + m0 * 4 + m1 * 8 + child_ref_ct * 16];
          if (!error_result) {
            continue;
          }
          trio_error_ct_cur_buf[trio_idx] += error_result & 0xffffff;
          ++variant_error_ct;
          if (!cur_variant_errors) {
            continue;
          }
          SetBit(trio_idx, cur_variant_errors);
          if (!cur_edescrips) {
            continue;
          }
          AlleleCode mendel_ecode = error_result >> 24;
          if (ignore_father) {
            if (mendel_ecode == 7) {
              mendel_ecode = 9;
            } else if (mendel_ecode == 4) {
              mendel_ecode = 10;
            }
          }
          if (ecode_col) {
            *edescrips_write_iter++ = mendel_ecode;
          }
          if (!edescrip_col) {
            continue;
          }
          edescrips_write_iter[0] = dad_ac0;
          edescrips_write_iter[1] = dad_ac1;
          edescrips_write_iter[2] = mom_ac0;
          edescrips_write_iter[3] = mom_ac1;
          edescrips_write_iter[4] = child_ac0;
          edescrips_write_iter[5] = child_ac1;
          if (is_x) {
            // render dad as haploid
            edescrips_write_iter[1] = kMissingAlleleCode - (dad_ac0 == kMissingAlleleCode);
          }
          switch (mendel_ecode) {
          case 3:
          case 6:
            edescrips_write_iter[2] = kMissingAlleleCode;
            edescrips_write_iter[3] = kMissingAlleleCode;
            break;
          case 4:
          case 7:
            edescrips_write_iter[0] = kMissingAlleleCode;
            edescrips_write_iter[1] = kMissingAlleleCode;
            break;
          case 9:
          case 10:
            // chrX male, or any chrM: render dad as '*', child as haploid
            // unless heterozygous.
            edescrips_write_iter[0] = kMissingAlleleCode;
            edescrips_write_iter[1] = kMissingAlleleCode - 1;
            // mom is always rendered as diploid on chrX; try to render as
            // haploid on chrM.
            if (is_mt && (mom_ac0 == mom_ac1)) {
              edescrips_write_iter[3] = kMissingAlleleCode;
            }
            if (child_ac0 == child_ac1) {
              edescrips_write_iter[5] = kMissingAlleleCode;
            }
            break;
          case 11:
          case 12:
            // chrY non-female: render dad and child as haploid, don't render
            // mom at all (don't even include preceding "x")
            edescrips_write_iter[1] = kMissingAlleleCode;
            edescrips_write_iter[2] = kMissingAlleleCode;
            edescrips_write_iter[3] = 0;
            edescrips_write_iter[5] = kMissingAlleleCode;
            break;
          }
          edescrips_write_iter = &(edescrips_write_iter[6]);
        }
      }
      variant_error_cts[cur_variant_bidx] = variant_error_ct;
      variant_obs_cts[cur_variant_bidx] = variant_error_ct_denom;
      if ((!var_first) || (variant_error_ct <= S_CAST(uint32_t, S_CAST(int32_t, u31tod(variant_error_ct_denom) * max_var_error)))) {
        // This variant counts toward per-trio stats.
        if (variant_error_ct) {
          if (var_first) {
            U32CastVecAdd(trio_error_ct_cur_buf, trio_vec_ct, trio_error_ct_acc_buf);
          }
          ++error_ct_fill;
          if (error_ct_fill == 255) {
            unsigned char* read_iter = R_CAST(unsigned char*, trio_error_ct_acc_buf);
            uint32_t* trio_error_cts_write_iter = trio_error_cts;
            for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
              trio_error_cts_write_iter[0] += read_iter[0];
              trio_error_cts_write_iter[1] += read_iter[1];
              trio_error_cts_write_iter[2] += read_iter[2];
              trio_error_cts_write_iter = &(trio_error_cts_write_iter[3]);
              read_iter = &(read_iter[4]);
            }
            error_ct_fill = 0;
            ZeroU32Arr(trio_ct, trio_error_ct_acc_buf);
            ZeroU32Arr(trio_ct, trio_missing_ct_acc_buf);
          }
        }
        if (variant_error_ct_denom != trio_ct) {
          if (var_first) {
            U32CastVecAdd(trio_missing_ct_cur_buf, trio_vec_ct, trio_missing_ct_acc_buf);
          }
          missing_ct_fill += exclude_missing_from_denom;
          if (missing_ct_fill == 255) {
            unsigned char* read_iter = R_CAST(unsigned char*, trio_missing_ct_acc_buf);
            uint32_t* trio_missing_cts_write_iter = trio_missing_cts;
            for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
              trio_missing_cts_write_iter[0] += read_iter[0];
              trio_missing_cts_write_iter[1] += read_iter[1];
              trio_missing_cts_write_iter[2] += read_iter[2];
              trio_missing_cts_write_iter = &(trio_missing_cts_write_iter[3]);
              read_iter = &(read_iter[4]);
            }
          }
        }
      }
      if (cur_variant_errors) {
        cur_variant_errors = &(cur_variant_errors[trio_ctl]);
        if (cur_edescrips) {
          cur_edescrips = &(cur_edescrips[edescrip_stride]);
        }
      }
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  if (error_ct_fill) {
    unsigned char* read_iter = R_CAST(unsigned char*, trio_error_ct_acc_buf);
    uint32_t* trio_error_cts_write_iter = trio_error_cts;
    for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
      trio_error_cts_write_iter[0] += read_iter[0];
      trio_error_cts_write_iter[1] += read_iter[1];
      trio_error_cts_write_iter[2] += read_iter[2];
      trio_error_cts_write_iter = &(trio_error_cts_write_iter[3]);
      read_iter = &(read_iter[4]);
    }
  }
  if (missing_ct_fill) {
    unsigned char* read_iter = R_CAST(unsigned char*, trio_missing_ct_acc_buf);
    uint32_t* trio_missing_cts_write_iter = trio_missing_cts;
    for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
      trio_missing_cts_write_iter[0] += read_iter[0];
      trio_missing_cts_write_iter[1] += read_iter[1];
      trio_missing_cts_write_iter[2] += read_iter[2];
      trio_missing_cts_write_iter = &(trio_missing_cts_write_iter[3]);
      read_iter = &(read_iter[4]);
    }
  }
  while (0) {
  MendelErrorScanThread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}

PglErr MendelErrorScan(const PedigreeIdInfo* piip, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const MendelInfo* mip, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t generate_reports, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uintptr_t* outer_sample_include, uintptr_t* variant_include, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* l_cswritep = nullptr;
  char* big_cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState l_css;
  CompressStreamState big_css;
  ThreadGroup tg;
  PreinitCstream(&l_css);
  PreinitCstream(&big_css);
  PreinitThreads(&tg);
  MendelErrorScanCtx ctx;
  {
    if (unlikely(IsSet(cip->haploid_mask, 0))) {
      logerrputs("Error: --me/--mendel cannot be used on haploid genomes.\n");
      goto MendelErrorScan_ret_INCONSISTENT_INPUT;
    }
    const MendelFlags flags = mip->flags;
    TrioFlags trio_flags = (flags & kfMendelDuos)? (kfTrioDuos | kfTrioPopulateIds) : kfTrioPopulateIds;
    uint32_t sid_col = 0;
    if (generate_reports) {
      sid_col = SidColIsRequired(piip->sii.sids, flags / kfMendelRptColMaybesid);
      if (sid_col) {
        trio_flags |= kfTrioPopulateSids;
      }
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* trio_sample_include;
    uint32_t* trio_sample_include_cumulative_popcounts;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &trio_sample_include) ||
                 bigstack_alloc_u32(raw_sample_ctl, &trio_sample_include_cumulative_popcounts))) {
      goto MendelErrorScan_ret_NOMEM;
    }
    FamilyInfo family_info;
    PreinitFamilyInfo(&family_info);
    reterr = GetTriosAndFamilies(outer_sample_include, piip, founder_info, sex_nm, sex_male, raw_sample_ct, trio_flags, &sample_ct, trio_sample_include, &family_info);
    if (unlikely(reterr)) {
      goto MendelErrorScan_ret_1;
    }
    if (!sample_ct) {
      logerrprintf("Warning: Skipping --me/--mendel since there are no %strios.\n", (flags & kfMendelDuos)? "duos or " : "");
      goto MendelErrorScan_ret_1;
    }
    FillCumulativePopcounts(trio_sample_include, raw_sample_ctl, trio_sample_include_cumulative_popcounts);

    char* chr_buf = nullptr;
    const uint32_t output_zst = (flags / kfMendelRptZs) & 1;
    const uint32_t big_report = generate_reports && (!(flags & kfMendelRptSummariesOnly));
    const uint32_t fid_col = FidColIsRequired(&(piip->sii), flags / kfMendelRptColMaybefid);
    const uint32_t chrom_col = (flags / kfMendelRptColChrom) & 1;
    const uint32_t pos_col = (flags / kfMendelRptColPos) & 1;
    const uint32_t ref_col = (flags / kfMendelRptColRef) & 1;
    const uint32_t alt_col = (flags / kfMendelRptColAlt) & 1;
    const uint32_t ecode_col = big_report && (flags & kfMendelRptColCode);
    const uint32_t edescrip_col = big_report && (flags & kfMendelRptColError);
    const uint32_t nl_col = (flags / kfMendelRptLcolN) & 1;
    const uint32_t nobsl_col = (flags / kfMendelRptLcolNobs) & 1;
    const uint32_t fracl_col = (flags / kfMendelRptLcolFrac) & 1;
    if (generate_reports) {
      uint32_t max_chr_blen = 0;
      if (chrom_col) {
        max_chr_blen = GetMaxChrSlen(cip) + 1;
        if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
          goto MendelErrorScan_ret_NOMEM;
        }
      }
      const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 128 + (ref_col | alt_col) * max_allele_slen;
      OutnameZstSet(".lmendel", output_zst, outname_end);
      reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &l_css, &l_cswritep);
      if (unlikely(reterr)) {
        goto MendelErrorScan_ret_1;
      }
      *l_cswritep++ = '#';
      if (chrom_col) {
        l_cswritep = strcpya_k(l_cswritep, "CHROM\t");
      }
      if (pos_col) {
        l_cswritep = strcpya_k(l_cswritep, "POS\t");
      }
      l_cswritep = strcpya_k(l_cswritep, "ID");
      if (ref_col) {
        l_cswritep = strcpya_k(l_cswritep, "\tREF");
      }
      if (alt_col) {
        l_cswritep = strcpya_k(l_cswritep, "\tALT");
      }
      if (nl_col) {
        l_cswritep = strcpya_k(l_cswritep, "\tN");
      }
      if (nobsl_col) {
        l_cswritep = strcpya_k(l_cswritep, "\tOBS_CT");
      }
      if (fracl_col) {
        l_cswritep = strcpya_k(l_cswritep, "\tF_ERR");
      }
      AppendBinaryEoln(&l_cswritep);

      if (big_report) {
        OutnameZstSet(".mendel", output_zst, outname_end);
        reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &big_css, &big_cswritep);
        if (unlikely(reterr)) {
          goto MendelErrorScan_ret_1;
        }
        *big_cswritep++ = '#';
        if (fid_col) {
          big_cswritep = strcpya_k(big_cswritep, "FID\t");
        }
        big_cswritep = strcpya_k(big_cswritep, "KID\t");
        if (sid_col) {
          big_cswritep = strcpya_k(big_cswritep, "KID_SID\t");
        }
        if (chrom_col) {
          big_cswritep = strcpya_k(big_cswritep, "CHROM\t");
        }
        if (pos_col) {
          big_cswritep = strcpya_k(big_cswritep, "POS\t");
        }
        big_cswritep = strcpya_k(big_cswritep, "ID");
        if (ref_col) {
          l_cswritep = strcpya_k(l_cswritep, "\tREF");
        }
        if (alt_col) {
          l_cswritep = strcpya_k(l_cswritep, "\tALT");
        }
        if (ecode_col) {
          big_cswritep = strcpya_k(big_cswritep, "\tCODE");
        }
        if (edescrip_col) {
          big_cswritep = strcpya_k(big_cswritep, "\tERROR");
        }
        AppendBinaryEoln(&big_cswritep);
      }
    }

    ctx.fip = &family_info;
    ctx.variant_include = variant_include;
    ctx.cip = cip;
    const uint32_t thread_mhc = (max_allele_ct > 2);
    ctx.allele_idx_offsets = thread_mhc? allele_idx_offsets : nullptr;
    ctx.sample_include = trio_sample_include;
    ctx.sample_include_cumulative_popcounts = trio_sample_include_cumulative_popcounts;
    const uintptr_t sample_ctl = BitCtToWordCt(sample_ct);
    {
      const uint32_t sample_ctv2 = NypCtToVecCt(sample_ct);
      uintptr_t* sex_male_collapsed;
      uintptr_t* sex_male_collapsed_interleaved;
      uintptr_t* sex_female_collapsed;
      uintptr_t* sex_female_collapsed_interleaved;
      if (unlikely(bigstack_alloc_w(sample_ctl, &sex_male_collapsed) ||
                   bigstack_alloc_w(sample_ctv2 * kWordsPerVec, &sex_male_collapsed_interleaved) ||
                   bigstack_alloc_w(sample_ctl, &sex_female_collapsed) ||
                   bigstack_alloc_w(sample_ctv2 * kWordsPerVec, &sex_female_collapsed_interleaved))) {
        goto MendelErrorScan_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, trio_sample_include, sample_ct, sex_male_collapsed);
      FillInterleavedMaskVec(sex_male_collapsed, sample_ctv2, sex_male_collapsed_interleaved);
      CopyBitarrSubset(sex_nm, trio_sample_include, sample_ct, sex_female_collapsed);
      BitvecInvmask(sex_male, sample_ctl, sex_female_collapsed);
      FillInterleavedMaskVec(sex_female_collapsed, sample_ctv2, sex_female_collapsed_interleaved);
      ctx.sex_male_collapsed = sex_male_collapsed;
      ctx.sex_male_collapsed_interleaved = sex_male_collapsed_interleaved;
      ctx.sex_female_collapsed = sex_female_collapsed;
      ctx.sex_female_collapsed_interleaved = sex_female_collapsed_interleaved;
    }
    ctx.flags = flags;
    ctx.sample_ct = sample_ct;
    ctx.male_ct = PopcountWords(ctx.sex_male_collapsed, sample_ctl);
    ctx.female_ct = PopcountWords(ctx.sex_female_collapsed, sample_ctl);
    const uint32_t max_difflist_len = sample_ct / kPglMaxDifflistLenDivisor;
    ctx.max_difflist_len = max_difflist_len;
    const uintptr_t raregeno_vec_ct = NypCtToVecCt(max_difflist_len);
    const uintptr_t difflist_sample_id_vec_ct = Int32CtToVecCt(max_difflist_len);

    const double max_var_error = mip->max_var_error * (1 + kSmallEpsilon);
    ctx.max_var_error = max_var_error;
    // not strictly needed, but it removes a potential source of confusion, and
    // memory isn't that tight here
    uintptr_t* new_variant_include = nullptr;
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    if (max_var_error < 1.0) {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, &new_variant_include))) {
        goto MendelErrorScan_ret_NOMEM;
      }
      memcpy(new_variant_include, variant_include, raw_variant_ctl * sizeof(intptr_t));
    }

    uint32_t calc_thread_ct = max_thread_ct - (max_thread_ct > 4);
    if (unlikely(bigstack_alloc_wp(calc_thread_ct, &ctx.raregenos) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.difflist_sample_id_bufs) ||
                 bigstack_alloc_wp(calc_thread_ct, &ctx.thread_trio_include_bufs) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_trio_error_ct_acc_bufs) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_trio_missing_ct_acc_bufs) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_trio_error_cts) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_trio_missing_cts))) {
      goto MendelErrorScan_ret_NOMEM;
    }
    const uintptr_t trio_ct = family_info.trio_ct;
    const uintptr_t trio_include_vec_ct = BitCtToVecCt(trio_ct);
    const uintptr_t trio_error_acc_vec_ct = Int32CtToVecCt(trio_ct);
    const uintptr_t trio_error_main_cts_vec_ct = Int32CtToVecCt(trio_ct * 3);
    // additional per-thread vec requirements:
    //   raregenos: raregeno_vec_ct
    //   difflist_sample_id_bufs: difflist_sample_id_vec_ct
    //   trio_include_bufs: trio_include_vec_ct
    //   trio_error_ct_acc_bufs: trio_error_acc_vec_ct
    //   trio_error_ct_missing_bufs: trio_error_acc_vec_ct
    //   if var_first, 2 * trio_error_acc_vec_ct for cur_bufs
    //   if thread_mhc, wide_codes_vec_ct
    //   trio_error_cts: trio_error_main_cts_vec_ct
    //   trio_missing_cts: trio_error_main_cts_vec_ct
    uintptr_t thread_xalloc_vec_ct = raregeno_vec_ct + difflist_sample_id_vec_ct + trio_include_vec_ct + 2 * (trio_error_acc_vec_ct + trio_error_main_cts_vec_ct);
    const uint32_t var_first = (flags / kfMendelFilterVarFirst) & 1;
    ctx.thread_trio_error_ct_cur_bufs = nullptr;
    ctx.thread_trio_missing_ct_cur_bufs = nullptr;
    if (var_first) {
      if (unlikely(bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_trio_error_ct_cur_bufs) ||
                   bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_trio_missing_ct_cur_bufs))) {
        goto MendelErrorScan_ret_NOMEM;
      }
      thread_xalloc_vec_ct += 2 * trio_error_acc_vec_ct;
    }
    const uint32_t duo_exists = family_info.duo_exists;
    const uintptr_t wide_codes_vec_ct = AlleleCodeCtToVecCt(2 * (sample_ct + duo_exists));
    ctx.thread_wide_code_bufs = nullptr;
    if (thread_mhc) {
      if (unlikely(bigstack_alloc_acp(calc_thread_ct, &ctx.thread_wide_code_bufs))) {
        goto MendelErrorScan_ret_NOMEM;
      }
      thread_xalloc_vec_ct += wide_codes_vec_ct;
    }
    const uintptr_t thread_xalloc_cacheline_ct = VecCtToCachelineCt(thread_xalloc_vec_ct);

    // additional per-variant requirements:
    //   variant_error_cts: 2 * sizeof(int32_t)
    //   variant_obs_cts: 2 * sizeof(int32_t)
    //   if big_report:
    //     variant_errors: 2 * trio_ctl * sizeof(intptr_t)
    //     edescrips: 2 * edescrip_stride * sizeof(AlleleCode)
    ctx.variant_errors[0] = nullptr;
    ctx.variant_errors[1] = nullptr;
    ctx.edescrips[0] = nullptr;
    ctx.edescrips[1] = nullptr;
    ctx.err_info = (~0LLU) << 32;
    const uintptr_t trio_ctl = BitCtToWordCt(trio_ct);
    const uintptr_t edescrip_stride = trio_ct * (ecode_col + 6 * edescrip_col);
    const uintptr_t per_variant_xalloc_byte_ct = 4 * sizeof(int32_t) + 2 * (big_report * trio_ctl * sizeof(intptr_t) + edescrip_stride * sizeof(AlleleCode));

    // defend against adverse rounding for per-variant allocations
    uintptr_t bytes_avail = bigstack_left();
    if (unlikely(bytes_avail < 7 * kCacheline)) {
      goto MendelErrorScan_ret_NOMEM;
    }
    bytes_avail -= 7 * kCacheline;

    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    ctx.thread_read_mhc = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct + duo_exists, variant_ct, bytes_avail, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, 0, pgfip, &calc_thread_ct, &ctx.genovecs, thread_mhc? (&ctx.thread_read_mhc) : nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto MendelErrorScan_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto MendelErrorScan_ret_NOMEM;
    }
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      unsigned char* cur_alloc = S_CAST(unsigned char*, bigstack_alloc_raw(thread_xalloc_cacheline_ct * kCacheline));
      ctx.raregenos[tidx] = R_CAST(uintptr_t*, cur_alloc);
      cur_alloc = &(cur_alloc[raregeno_vec_ct * kBytesPerVec]);
      ctx.difflist_sample_id_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[difflist_sample_id_vec_ct * kBytesPerVec]);
      ctx.thread_trio_include_bufs[tidx] = R_CAST(uintptr_t*, cur_alloc);
      cur_alloc = &(cur_alloc[trio_include_vec_ct * kBytesPerVec]);
      ctx.thread_trio_error_ct_acc_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[trio_error_acc_vec_ct * kBytesPerVec]);
      ctx.thread_trio_missing_ct_acc_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[trio_error_acc_vec_ct * kBytesPerVec]);
      if (var_first) {
        ctx.thread_trio_error_ct_cur_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[trio_error_acc_vec_ct * kBytesPerVec]);
        ctx.thread_trio_missing_ct_cur_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[trio_error_acc_vec_ct * kBytesPerVec]);
      }
      if (thread_mhc) {
        ctx.thread_wide_code_bufs[tidx] = R_CAST(AlleleCode*, cur_alloc);
        cur_alloc = &(cur_alloc[wide_codes_vec_ct * kBytesPerVec]);
      }
      ctx.thread_trio_error_cts[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[trio_error_main_cts_vec_ct * kBytesPerVec]);
      ctx.thread_trio_missing_cts[tidx] = R_CAST(uint32_t*, cur_alloc);
      assert(&(cur_alloc[trio_error_main_cts_vec_ct * kBytesPerVec]) <= g_bigstack_base);
    }
    if (unlikely(bigstack_alloc_u32(read_block_size, &(ctx.variant_error_cts[0])) ||
                 bigstack_alloc_u32(read_block_size, &(ctx.variant_error_cts[1])) ||
                 bigstack_alloc_u32(read_block_size, &(ctx.variant_obs_cts[0])) ||
                 bigstack_alloc_u32(read_block_size, &(ctx.variant_obs_cts[1])))) {
      assert(0);
      goto MendelErrorScan_ret_NOMEM;
    }
    if (big_report) {
      if (unlikely(bigstack_alloc_w(trio_ctl * read_block_size, &(ctx.variant_errors[0])) ||
                   bigstack_alloc_w(trio_ctl * read_block_size, &(ctx.variant_errors[1])))) {
        assert(0);
        goto MendelErrorScan_ret_NOMEM;
      }
      if (edescrip_stride) {
        if (unlikely(bigstack_alloc_ac(edescrip_stride * read_block_size, &(ctx.edescrips[0])) ||
                     bigstack_alloc_ac(edescrip_stride * read_block_size, &(ctx.edescrips[1])))) {
          assert(0);
          goto MendelErrorScan_ret_NOMEM;
        }
      }
    }

    SetThreadFuncAndData(MendelErrorScanThread, &ctx, &tg);
    uint64_t tot_error_ct = 0;
    fputs("--me/--mendel: 0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    const uint64_t* trio_list = family_info.trio_list;
    const char* trio_fids = family_info.trio_fids;
    const char* iids = family_info.iids;
    const char* sids = family_info.sids;
    const uintptr_t max_fid_blen = family_info.max_fid_blen;
    const uintptr_t max_iid_blen = family_info.max_iid_blen;
    const uintptr_t max_sid_blen = family_info.max_sid_blen;
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t write_variant_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t prev_block_variant_ct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto MendelErrorScan_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, ctx.err_info >> 32);
          goto MendelErrorScan_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_size = cur_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_size, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_size == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto MendelErrorScan_ret_THREAD_CREATE_FAIL;
        }
      }

      parity = 1 - parity;
      if (variant_idx) {
        // process *previous* block results
        const uint32_t* block_variant_error_cts = ctx.variant_error_cts[parity];
        const uint32_t* block_variant_obs_cts = ctx.variant_obs_cts[parity];
        const uintptr_t* block_variant_errors = ctx.variant_errors[parity];
        const AlleleCode* block_edescrips = ctx.edescrips[parity];
        uint32_t cur_allele_ct = 2;
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &write_variant_bits);
          if (chr_buf && (write_variant_uidx >= chr_end)) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
            *chr_name_end = '\t';
            chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
          }
          const uint32_t cur_error_ct = block_variant_error_cts[variant_bidx];
          const uint32_t cur_obs_ct = block_variant_obs_cts[variant_bidx];
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            cur_allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offset_base;
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          if (generate_reports) {
            if (chrom_col) {
              l_cswritep = memcpya(l_cswritep, chr_buf, chr_buf_blen);
            }
            if (pos_col) {
              l_cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', l_cswritep);
            }
            l_cswritep = strcpya(l_cswritep, variant_ids[write_variant_uidx]);
            if (ref_col) {
              *l_cswritep++ = '\t';
              l_cswritep = strcpya(l_cswritep, cur_alleles[0]);
            }
            if (alt_col) {
              *l_cswritep++ = '\t';
              for (uint32_t allele_idx = 1; allele_idx != cur_allele_ct; ++allele_idx) {
                if (unlikely(Cswrite(&l_css, &l_cswritep))) {
                  goto MendelErrorScan_ret_WRITE_FAIL;
                }
                l_cswritep = strcpyax(l_cswritep, cur_alleles[allele_idx], ',');
              }
              --l_cswritep;
            }
            if (nl_col) {
              *l_cswritep++ = '\t';
              l_cswritep = u32toa(cur_error_ct, l_cswritep);
            }
            if (nobsl_col) {
              *l_cswritep++ = '\t';
              l_cswritep = u32toa(cur_obs_ct, l_cswritep);
            }
            if (fracl_col) {
              *l_cswritep++ = '\t';
              l_cswritep = dtoa_g(u31tod(cur_error_ct) / u31tod(cur_obs_ct), l_cswritep);
            }
            AppendBinaryEoln(&l_cswritep);
            if (unlikely(Cswrite(&l_css, &l_cswritep))) {
              goto MendelErrorScan_ret_WRITE_FAIL;
            }
          }
          if (!cur_error_ct) {
            continue;
          }
          tot_error_ct += cur_error_ct;
          if (max_var_error * u31tod(cur_obs_ct) < u31tod(cur_error_ct)) {
            ClearBit(write_variant_uidx, new_variant_include);
          }
          if (!big_report) {
            continue;
          }
          const uintptr_t* cur_variant_errors = &(block_variant_errors[variant_bidx * trio_ctl]);
          const AlleleCode* edescrips_read_iter = &(block_edescrips[variant_bidx * edescrip_stride]);
          uintptr_t trio_idx_base = 0;
          uintptr_t cur_bits = cur_variant_errors[0];
          for (uintptr_t error_idx = 0; error_idx != cur_error_ct; ++error_idx) {
            const uintptr_t trio_idx = BitIter1(cur_variant_errors, &trio_idx_base, &cur_bits);
            if (fid_col) {
              big_cswritep = strcpyax(big_cswritep, &(trio_fids[trio_idx * max_fid_blen]), '\t');
            }
            const uint32_t child_idx = S_CAST(uint32_t, trio_list[trio_idx]);
            big_cswritep = strcpyax(big_cswritep, &(iids[child_idx * max_iid_blen]), '\t');
            if (sid_col) {
              big_cswritep = strcpyax(big_cswritep, &(sids[child_idx * max_sid_blen]), '\t');
            }
            if (chrom_col) {
              big_cswritep = memcpya(big_cswritep, chr_buf, chr_buf_blen);
            }
            if (pos_col) {
              big_cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', big_cswritep);
            }
            big_cswritep = strcpya(big_cswritep, variant_ids[write_variant_uidx]);
            if (ref_col) {
              *big_cswritep++ = '\t';
              big_cswritep = strcpya(big_cswritep, cur_alleles[0]);
            }
            if (alt_col) {
              *big_cswritep++ = '\t';
              for (uint32_t allele_idx = 1; allele_idx != cur_allele_ct; ++allele_idx) {
                if (unlikely(Cswrite(&big_css, &big_cswritep))) {
                  goto MendelErrorScan_ret_WRITE_FAIL;
                }
                big_cswritep = strcpyax(big_cswritep, cur_alleles[allele_idx], ',');
              }
              --big_cswritep;
            }
            if (ecode_col) {
              *big_cswritep++ = '\t';
              const uint32_t ecode = *edescrips_read_iter++;
              big_cswritep = u32toa(ecode, big_cswritep);
            }
            if (edescrip_col) {
              *big_cswritep++ = '\t';
              const uint32_t dad_ac0 = edescrips_read_iter[0];
              const uint32_t dad_ac1 = edescrips_read_iter[1];
              if (dad_ac0 == kMissingAlleleCode) {
                *big_cswritep++ = '*';
                if (dad_ac1 == kMissingAlleleCode) {
                  big_cswritep = strcpya_k(big_cswritep, "/*");
                }
              } else {
                big_cswritep = u32toa(dad_ac0, big_cswritep);
                if (dad_ac1 != kMissingAlleleCode) {
                  *big_cswritep++ = '/';
                  big_cswritep = u32toa(dad_ac1, big_cswritep);
                }
              }
              const uint32_t mom_ac0 = edescrips_read_iter[2];
              const uint32_t mom_ac1 = edescrips_read_iter[3];
              if (mom_ac0 == kMissingAlleleCode) {
                if (mom_ac1 != 0) {
                  big_cswritep = strcpya_k(big_cswritep, "x*");
                  if (mom_ac1 == kMissingAlleleCode) {
                    big_cswritep = strcpya_k(big_cswritep, "/*");
                  }
                }
              } else {
                *big_cswritep++ = 'x';
                big_cswritep = u32toa(mom_ac0, big_cswritep);
                if (mom_ac1 != kMissingAlleleCode) {
                  *big_cswritep++ = '/';
                  big_cswritep = u32toa(mom_ac1, big_cswritep);
                }
              }
              big_cswritep = strcpya_k(big_cswritep, "->");
              const uint32_t child_ac0 = edescrips_read_iter[4];
              const uint32_t child_ac1 = edescrips_read_iter[5];
              big_cswritep = u32toa(child_ac0, big_cswritep);
              if (child_ac1 != kMissingAlleleCode) {
                *big_cswritep++ = '/';
                big_cswritep = u32toa(child_ac1, big_cswritep);
              }
              edescrips_read_iter = &(edescrips_read_iter[6]);
            }
            AppendBinaryEoln(&big_cswritep);
            if (unlikely(Cswrite(&big_css, &big_cswritep))) {
              goto MendelErrorScan_ret_WRITE_FAIL;
            }
          }
        }
      }
      if (variant_idx == variant_ct) {
        break;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }
      ++read_block_idx;
      prev_block_variant_ct = cur_block_size;
      variant_idx += cur_block_size;
      pgfip->block_base = main_loadbufs[parity];
    }

    putc_unlocked('\r', stdout);
    logprintf("--me/--mendel: %" PRIu64 " Mendel error%s detected.\n", tot_error_ct, (tot_error_ct == 1)? "" : "s");
    uint32_t* trio_error_cts = ctx.thread_trio_error_cts[0];
    uint32_t* trio_missing_cts = ctx.thread_trio_missing_cts[0];
    {
      for (uint32_t thread_idx = 1; thread_idx != calc_thread_ct; ++thread_idx) {
        U32CastVecAdd(ctx.thread_trio_error_cts[thread_idx], trio_error_main_cts_vec_ct, trio_error_cts);
        U32CastVecAdd(ctx.thread_trio_missing_cts[thread_idx], trio_error_main_cts_vec_ct, trio_missing_cts);
      }
    }
    const uint64_t* family_list = family_info.family_list;
    if (generate_reports) {
      if (unlikely(CswriteCloseNull(&l_css, l_cswritep))) {
        goto MendelErrorScan_ret_WRITE_FAIL;
      }
      snprintf(outname_end, kMaxOutfnameExtBlen, ".fmendel");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
        goto MendelErrorScan_ret_OPEN_FAIL;
      }
      char* write_iter = g_textbuf;
      char* textbuf_flush = &(write_iter[kMaxMediumLine]);
      *write_iter++ = '#';
      if (fid_col) {
        write_iter = strcpya_k(write_iter, "FID\t");
      }
      write_iter = strcpya_k(write_iter, "PAT\tMAT\tCHLD\tN" EOLN_STR);
      const uint32_t family_ct = family_info.family_ct;
      uintptr_t trio_idx = 0;
      for (uint32_t family_idx = 0; family_idx != family_ct; ++family_idx) {
        assert((trio_list[trio_idx] >> 32) == family_idx);
        if (fid_col) {
          write_iter = strcpyax(write_iter, &(trio_fids[trio_idx * max_fid_blen]), '\t');
        }
        const uint64_t dad_and_mom_idxs = family_list[family_idx];
        const uint32_t dad_idx = S_CAST(uint32_t, dad_and_mom_idxs);
        const uint32_t mom_idx = dad_and_mom_idxs >> 32;
        write_iter = strcpyax(write_iter, &(iids[dad_idx * max_iid_blen]), '\t');
        write_iter = strcpyax(write_iter, &(iids[mom_idx * max_iid_blen]), '\t');
        uint64_t family_error_ct = 0;
        const uintptr_t trio_idx_start = trio_idx;
        do {
          family_error_ct += trio_error_cts[3 * trio_idx];
          if (++trio_idx == trio_ct) {
            break;
          }
        } while ((trio_list[trio_idx] >> 32) == family_idx);
        write_iter = u32toa_x(trio_idx - trio_idx_start, '\t', write_iter);
        write_iter = i64toa(family_error_ct, write_iter);
        AppendBinaryEoln(&write_iter);
        if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
          goto MendelErrorScan_ret_WRITE_FAIL;
        }
      }
      if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
        goto MendelErrorScan_ret_WRITE_FAIL;
      }

      outname_end[1] = 'i';
      if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
        goto MendelErrorScan_ret_OPEN_FAIL;
      }
      write_iter = g_textbuf;
      *write_iter++ = '#';
      const uint32_t trionum_col = (flags / kfMendelRptIcolTrionum) & 1;
      const uint32_t ni_col = (flags / kfMendelRptIcolN) & 1;
      const uint32_t nobsi_col = (flags / kfMendelRptIcolNobs) & 1;
      const uint32_t fraci_col = (flags / kfMendelRptIcolFrac) & 1;
      if (trionum_col) {
        write_iter = strcpya_k(write_iter, "TRIO_NUM\t");
      }
      if (fid_col) {
        write_iter = strcpya_k(write_iter, "FID\t");
      }
      write_iter = strcpya_k(write_iter, "IID");
      if (sid_col) {
        write_iter = strcpya_k(write_iter, "\tSID");
      }
      if (ni_col) {
        write_iter = strcpya_k(write_iter, "\tN");
      }
      if (nobsi_col) {
        write_iter = strcpya_k(write_iter, "\tOBS_CT");
      }
      if (fraci_col) {
        write_iter = strcpya_k(write_iter, "\tF_ERR");
      }
      AppendBinaryEoln(&write_iter);
      for (trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
        const char* trio_fid = &(trio_fids[trio_idx * max_fid_blen]);
        const uint64_t child_and_family_idxs = trio_list[trio_idx];
        const uint32_t child_idx = S_CAST(uint32_t, child_and_family_idxs);
        const uint32_t family_idx = child_and_family_idxs >> 32;
        const uint64_t dad_and_mom_idxs = family_list[family_idx];
        const uint32_t dad_idx = S_CAST(uint32_t, dad_and_mom_idxs);
        const uint32_t mom_idx = dad_and_mom_idxs >> 32;
        uint32_t sample_idxs[3];
        uint32_t cur_error_cts[3];
        uint32_t cur_denoms[3];
        uint32_t trio_member_ct;
        {
          uint32_t trio_member_idx = 0;
          if (dad_idx != sample_ct) {
            sample_idxs[0] = dad_idx;
            cur_error_cts[0] = trio_error_cts[3 * trio_idx + 1];
            cur_denoms[0] = variant_ct - trio_missing_cts[3 * trio_idx + 1];
            ++trio_member_idx;
          }
          if (mom_idx != sample_ct) {
            sample_idxs[trio_member_idx] = mom_idx;
            cur_error_cts[trio_member_idx] = trio_error_cts[3 * trio_idx + 2];
            cur_denoms[trio_member_idx] = variant_ct - trio_missing_cts[3 * trio_idx + 2];
            ++trio_member_idx;
          }
          sample_idxs[trio_member_idx] = child_idx;
          cur_error_cts[trio_member_idx] = trio_error_cts[3 * trio_idx];
          cur_denoms[trio_member_idx] = variant_ct - trio_missing_cts[3 * trio_idx];
          trio_member_ct = trio_member_idx + 1;
        }
        for (uint32_t trio_member_idx = 0; trio_member_idx != trio_member_ct; ++trio_member_idx) {
          if (trionum_col) {
            write_iter = u32toa_x(trio_idx + 1, '\t', write_iter);
          }
          if (fid_col) {
            write_iter = strcpyax(write_iter, trio_fid, '\t');
          }
          const uintptr_t sample_idx = sample_idxs[trio_member_idx];
          write_iter = strcpya(write_iter, &(iids[sample_idx * max_iid_blen]));
          if (sid_col) {
            *write_iter++ = '\t';
            write_iter = strcpya(write_iter, &(sids[sample_idx * max_sid_blen]));
          }
          const uint32_t cur_error_ct = cur_error_cts[trio_member_idx];
          if (ni_col) {
            *write_iter++ = '\t';
            write_iter = u32toa(cur_error_cts[trio_member_idx], write_iter);
          }
          const uint32_t cur_denom = cur_denoms[trio_member_idx];
          if (nobsi_col) {
            *write_iter++ = '\t';
            write_iter = u32toa(cur_denom, write_iter);
          }
          if (fraci_col) {
            *write_iter++ = '\t';
            write_iter = dtoa_g(u31tod(cur_error_ct) / u31tod(cur_denom), write_iter);
          }
          AppendBinaryEoln(&write_iter);
          if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
            goto MendelErrorScan_ret_WRITE_FAIL;
          }
        }
      }
      if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
        goto MendelErrorScan_ret_WRITE_FAIL;
      }

      *outname_end = '\0';
      if (big_report) {
        if (unlikely(CswriteCloseNull(&big_css, big_cswritep))) {
          goto MendelErrorScan_ret_WRITE_FAIL;
        }
        logprintfww("--mendel: Reports written to %s.mendel%s + %s.imendel + %s.fmendel + %s.lmendel%s .\n", outname, output_zst? ".zst" : "", outname, outname, outname, output_zst? ".zst" : "");
      } else {
        logprintfww("--mendel: Reports written to %s.imendel + %s.fmendel + %s.lmendel%s .\n", outname, outname, outname, output_zst? ".zst" : "");
      }
    }
    const double max_trio_error = mip->max_trio_error * (1 + kSmallEpsilon);
    if (max_trio_error < 1.0) {
      const double exclude_one_ratio = mip->exclude_one_ratio * (1 + kSmallEpsilon);
      for (uintptr_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
        if (u31tod(trio_error_cts[3 * trio_idx]) <= u31tod(variant_ct - trio_missing_cts[3 * trio_idx]) * max_trio_error) {
          continue;
        }
        const uint64_t child_and_family_idxs = trio_list[trio_idx];
        const uint32_t child_idx = S_CAST(uint32_t, child_and_family_idxs);
        const uint32_t family_idx = child_and_family_idxs >> 32;
        const uint64_t dad_and_mom_idxs = family_list[family_idx];
        const uint32_t dad_idx = S_CAST(uint32_t, dad_and_mom_idxs);
        const uint32_t mom_idx = dad_and_mom_idxs >> 32;
        if ((exclude_one_ratio <= 0.0) || (dad_idx == sample_ct) || (mom_idx == sample_ct)) {
          const uint32_t child_uidx = IdxToUidx(trio_sample_include, trio_sample_include_cumulative_popcounts, 0, raw_sample_ctl, child_idx);
          ClearBit(child_uidx, outer_sample_include);
          if (exclude_one_ratio == 0.0) {
            if (dad_idx != sample_ct) {
              const uint32_t dad_uidx = IdxToUidx(trio_sample_include, trio_sample_include_cumulative_popcounts, 0, raw_sample_ctl, dad_idx);
              ClearBit(dad_uidx, outer_sample_include);
            }
            if (mom_idx != sample_ct) {
              const uint32_t mom_uidx = IdxToUidx(trio_sample_include, trio_sample_include_cumulative_popcounts, 0, raw_sample_ctl, mom_idx);
              ClearBit(mom_uidx, outer_sample_include);
            }
          }
        } else {
          const uintptr_t dad_numer = trio_error_cts[3 * trio_idx + 1];
          const uint64_t dad_denom = variant_ct - trio_missing_cts[3 * trio_idx + 1];
          const uintptr_t mom_numer = trio_error_cts[3 * trio_idx + 2];
          const uint64_t mom_denom = variant_ct - trio_missing_cts[3 * trio_idx + 2];
          // both rates multiplied by dad_denom * mom_denom to avoid nan
          const double dad_scaled_error = u63tod(dad_numer * mom_denom);
          const double mom_scaled_error = u63tod(mom_numer * dad_denom);
          uint32_t clear_idx = child_idx;
          if (dad_scaled_error > mom_scaled_error * exclude_one_ratio) {
            clear_idx = dad_idx;
          } else if (mom_scaled_error > dad_scaled_error * exclude_one_ratio) {
            clear_idx = mom_idx;
          }
          const uint32_t clear_uidx = IdxToUidx(trio_sample_include, trio_sample_include_cumulative_popcounts, 0, raw_sample_ctl, clear_idx);
          ClearBit(clear_uidx, outer_sample_include);
        }
      }
    }
    if (new_variant_include) {
      memcpy(variant_include, new_variant_include, raw_variant_ctl * sizeof(intptr_t));
    }
    // caller prints --me summary
  }
  while (0) {
  MendelErrorScan_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MendelErrorScan_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  MendelErrorScan_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  MendelErrorScan_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  MendelErrorScan_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  MendelErrorScan_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 MendelErrorScan_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&big_css, big_cswritep);
  CswriteCloseCond(&l_css, l_cswritep);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

uint32_t EraseMendelErrors(const FamilyInfo* fip, const uintptr_t* sex_male_collapsed, const uintptr_t* sex_male_collapsed_interleaved, const uintptr_t* sex_female_collapsed, const uintptr_t* sex_female_collapsed_interleaved, uint32_t sample_ct, uint32_t male_ct, uint32_t female_ct, uint32_t is_x, uint32_t is_y, uint32_t is_mt, uintptr_t* genoarr, uint32_t* patch_01_ctp, uintptr_t* patch_01_set, AlleleCode* patch_01_vals, uint32_t* patch_10_ctp, uintptr_t* patch_10_set, AlleleCode* patch_10_vals, uintptr_t* erase_map, uintptr_t* genovec_buf, AlleleCode* wide_codes_buf) {
  const uint32_t* trio_lookup_iter = fip->trio_lookup;
  const uint32_t patch_01_ct = *patch_01_ctp;
  const uint32_t patch_10_ct = *patch_10_ctp;
  const uintptr_t trio_ct = fip->trio_ct;
  const uint32_t duo_exists = fip->duo_exists;
  uint32_t variant_error_ct = 0;
  if ((!patch_01_ct) && (!patch_10_ct)) {
    // ok for wide_codes_buf to be nullptr in this case
    uintptr_t* genoarr_read = genoarr;
    const uint32_t* biallelic_mendel_error_table = kBiallelicMendelErrorTableAutosomalOrX;
    if ((!fip->ordered_erase_ok) || duo_exists || is_x || is_y) {
      memcpy(genovec_buf, genoarr, NypCtToWordCt(sample_ct) * sizeof(intptr_t));
      if (is_x) {
        SetMaleHetMissing(sex_male_collapsed_interleaved, NypCtToVecCt(sample_ct), genovec_buf);
      } else if (is_y) {
        InterleavedSetMissing(sex_female_collapsed_interleaved, NypCtToVecCt(sample_ct), genovec_buf);
        biallelic_mendel_error_table = kBiallelicMendelErrorTableChrY;
      }
      if (duo_exists) {
        SetNyparrEntryTo3(sample_ct, genovec_buf);
      }
      genoarr_read = genovec_buf;
    }
    if (is_mt) {
      biallelic_mendel_error_table = kBiallelicMendelErrorTableChrM;
    }
    for (uintptr_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
      const uint32_t child_idx = trio_lookup_iter[0];
      trio_lookup_iter = &(trio_lookup_iter[3]);
      const uint32_t child_geno = GetNyparrEntry(genoarr_read, child_idx);
      if (child_geno == 3) {
        continue;
      }
      const uint32_t dad_idx = trio_lookup_iter[-2];
      const uint32_t mom_idx = trio_lookup_iter[-1];
      uint32_t dad_geno;
      if (is_x && IsSet(sex_male_collapsed, child_idx)) {
        dad_geno = 3;
      } else {
        dad_geno = GetNyparrEntry(genoarr_read, dad_idx);
      }
      const uint32_t mom_geno = GetNyparrEntry(genoarr_read, mom_idx);
      const uint32_t error_result = biallelic_mendel_error_table[dad_geno + mom_geno * 4 + child_geno * 16];
      if (!error_result) {
        continue;
      }
      SetNyparrEntryTo3(child_idx, genoarr);
      if (erase_map) {
        SetBit(child_idx, erase_map);
      }
      if (error_result & 0x100) {
        SetNyparrEntryTo3(dad_idx, genoarr);
        if (erase_map) {
          SetBit(dad_idx, erase_map);
        }
      }
      if (error_result & 0x10000) {
        SetNyparrEntryTo3(mom_idx, genoarr);
        if (erase_map) {
          SetBit(mom_idx, erase_map);
        }
      }
      ++variant_error_ct;
    }
    return variant_error_ct;
  }
  const uint16_t* hc_to_allele_codes = is_y? kHetMissingToAlleleCodes : kHcToAlleleCodes;
  InitMultiallelicWideCodes(hc_to_allele_codes, sex_male_collapsed, sex_female_collapsed, genoarr, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, sample_ct, male_ct, female_ct, patch_01_ct, patch_10_ct, is_x, is_y, wide_codes_buf);
  if (duo_exists) {
    wide_codes_buf[2 * sample_ct] = kMissingAlleleCode;
    wide_codes_buf[2 * sample_ct + 1] = kMissingAlleleCode;
  }
  const uint32_t* multiallelic_mendel_error_table = kMultiallelicMendelErrorTableAutosomalOrX;
  if (is_y) {
    multiallelic_mendel_error_table = kMultiallelicMendelErrorTableChrY;
  } else if (is_mt) {
    multiallelic_mendel_error_table = kMultiallelicMendelErrorTableChrM;
  }
  for (uintptr_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
    const uint32_t child_idx = trio_lookup_iter[0];
    trio_lookup_iter = &(trio_lookup_iter[3]);
    const AlleleCode child_ac0 = wide_codes_buf[child_idx * 2];
    if (child_ac0 == kMissingAlleleCode) {
      continue;
    }
    const AlleleCode child_ac1 = wide_codes_buf[child_idx * 2 + 1];
    const uint32_t dad_idx = trio_lookup_iter[-2];
    const uint32_t mom_idx = trio_lookup_iter[-1];
    const AlleleCode dad_ac0 = wide_codes_buf[dad_idx * 2];
    uint32_t d0;
    uint32_t d1;
    if ((is_x && IsSet(sex_male_collapsed, child_idx)) || (dad_ac0 == kMissingAlleleCode)) {
      d0 = 1;
      d1 = 1;
    } else {
      const AlleleCode dad_ac1 = wide_codes_buf[dad_idx * 2 + 1];
      d0 = (dad_ac0 == child_ac0) || (dad_ac1 == child_ac0);
      d1 = (dad_ac0 == child_ac1) || (dad_ac1 == child_ac1);
    }
    const AlleleCode mom_ac0 = wide_codes_buf[mom_idx * 2];
    const AlleleCode mom_ac1 = wide_codes_buf[mom_idx * 2 + 1];
    const uint32_t m0 = (mom_ac0 == child_ac0) || (mom_ac1 == child_ac0) || (mom_ac0 == kMissingAlleleCode);
    const uint32_t m1 = (mom_ac0 == child_ac1) || (mom_ac1 == child_ac1) || (mom_ac0 == kMissingAlleleCode);
    // child_ref_ct only affects error code numbers, we don't need it here.
    const uint32_t error_result = multiallelic_mendel_error_table[d0 + d1 * 2 + m0 * 4 + m1 * 8];
    if (!error_result) {
      continue;
    }
    SetNyparrEntryTo3(child_idx, genoarr);
    if (erase_map) {
      SetBit(child_idx, erase_map);
    }
    if (error_result & 0x100) {
      SetNyparrEntryTo3(dad_idx, genoarr);
      if (erase_map) {
        SetBit(dad_idx, erase_map);
      }
    }
    if (error_result & 0x10000) {
      SetNyparrEntryTo3(mom_idx, genoarr);
      if (erase_map) {
        SetBit(mom_idx, erase_map);
      }
    }
    ++variant_error_ct;
  }
  if (!variant_error_ct) {
    return 0;
  }
  if (patch_01_ct) {
    ClearGenoarrMissing1bit8Unsafe(genoarr, patch_01_ctp, patch_01_set, patch_01_vals);
  }
  if (patch_10_ct) {
    ClearGenoarrMissing1bit16Unsafe(genoarr, patch_10_ctp, patch_10_set, patch_10_vals);
  }
  return variant_error_ct;
}

// QFAM family-based association tests for quantitative traits, PLINK 1.9's
// --qfam / --qfam-parents / --qfam-between / --qfam-total.
//
// Each family, sibship or unrelated singleton gets a between-family genotype
// score B (its mean ALT dosage, centered at 1), and each sample gets a
// within-family deviate W = (own ALT dosage - 1) - B.  The four tests regress
// the phenotype on W (--qfam and --qfam-parents), on B (--qfam-between), or on
// B + W (--qfam-total).  Only the empirical p-value from the permutation is
// interpretable: the asymptotic one ignores the family structure, which is why
// PLINK 1.9 names that column RAW_P, and this port keeps the name.
//
// The permutation permutes B across families and flips the sign of W within
// them, so the null it tests is "this variant's transmission is unrelated to
// the phenotype" rather than "the genotypes are exchangeable".

static const double kQfamEpsilon = 0.000000000931322574615478515625;
static const double kQfamSmallEpsilon = 0.00000000000005684341886080801486968994140625;

typedef struct QfamGroupsStruct {
  // fss = families, then sibships, then singletons.
  // fss_starts[k] is where group k's members begin in fss_contents; a family's
  // first two members are its parents.
  uint32_t* fss_starts;
  uint32_t* fss_contents;
  // For each of the lm_ct samples in the regression, in increasing sample_idx
  // order, the group it belongs to.
  uint32_t* sample_lm_to_fss_idx;
  uintptr_t* lm_eligible;
  // --qfam-parents only: founder parents whose partner's genotype has to be
  // present too.  nullptr otherwise.
  uintptr_t* lm_within2_founder;
  uint32_t family_ct;
  uint32_t fs_ct;
  uint32_t fss_ct;
  uint32_t lm_ct;
} QfamGroups;

// B and W for one variant, plus the phenotype sums with the samples whose
// genotype is missing taken back out.
void QfamComputeBw(const uintptr_t* genovec, const QfamGroups* qgp, const double* pheno_d2, double qt_sum_all, double qt_ssq_all, uintptr_t* nm_fss, uintptr_t* nm_lm, double* qfam_b, double* qfam_w, double* qt_sum_ptr, double* qt_ssq_ptr) {
  const uint32_t* fss_starts = qgp->fss_starts;
  const uint32_t* fss_contents = qgp->fss_contents;
  const uint32_t family_ct = qgp->family_ct;
  const uint32_t fs_ct = qgp->fs_ct;
  const uint32_t fss_ct = qgp->fss_ct;
  const uint32_t lm_ct = qgp->lm_ct;
  double qt_sum = qt_sum_all;
  double qt_ssq = qt_ssq_all;
  SetAllBits(fss_ct, nm_fss);
  for (uint32_t fss_idx = 0; fss_idx != family_ct; ++fss_idx) {
    const uint32_t cur_start = fss_starts[fss_idx];
    const uint32_t cur_end = fss_starts[fss_idx + 1];
    const uint32_t dad_geno = GetNyparrEntry(genovec, fss_contents[cur_start]);
    const uint32_t mom_geno = GetNyparrEntry(genovec, fss_contents[cur_start + 1]);
    if ((dad_geno != 3) && (mom_geno != 3)) {
      qfam_b[fss_idx] = 0.5 * S_CAST(double, S_CAST(int32_t, dad_geno + mom_geno)) - 1.0;
      continue;
    }
    // Fall back on the children's mean when a parent's genotype is missing.
    uint32_t sib_ct = 0;
    uint32_t alt_ct = 0;
    for (uint32_t uii = cur_start + 2; uii != cur_end; ++uii) {
      const uint32_t cur_geno = GetNyparrEntry(genovec, fss_contents[uii]);
      if (cur_geno != 3) {
        ++sib_ct;
        alt_ct += cur_geno;
      }
    }
    if (sib_ct) {
      qfam_b[fss_idx] = S_CAST(double, S_CAST(int32_t, alt_ct)) / u31tod(sib_ct) - 1.0;
    } else {
      ClearBit(fss_idx, nm_fss);
    }
  }
  for (uint32_t fss_idx = family_ct; fss_idx != fs_ct; ++fss_idx) {
    const uint32_t cur_start = fss_starts[fss_idx];
    const uint32_t cur_end = fss_starts[fss_idx + 1];
    uint32_t sib_ct = 0;
    uint32_t alt_ct = 0;
    for (uint32_t uii = cur_start; uii != cur_end; ++uii) {
      const uint32_t cur_geno = GetNyparrEntry(genovec, fss_contents[uii]);
      if (cur_geno != 3) {
        ++sib_ct;
        alt_ct += cur_geno;
      }
    }
    if (sib_ct) {
      qfam_b[fss_idx] = S_CAST(double, S_CAST(int32_t, alt_ct)) / u31tod(sib_ct) - 1.0;
    } else {
      ClearBit(fss_idx, nm_fss);
    }
  }
  for (uint32_t fss_idx = fs_ct; fss_idx != fss_ct; ++fss_idx) {
    const uint32_t cur_geno = GetNyparrEntry(genovec, fss_contents[fss_starts[fss_idx]]);
    if (cur_geno != 3) {
      qfam_b[fss_idx] = S_CAST(double, S_CAST(int32_t, cur_geno)) - 1.0;
    } else {
      ClearBit(fss_idx, nm_fss);
    }
  }
  const uintptr_t* lm_eligible = qgp->lm_eligible;
  const uintptr_t* lm_within2_founder = qgp->lm_within2_founder;
  const uint32_t* sample_lm_to_fss_idx = qgp->sample_lm_to_fss_idx;
  SetAllBits(lm_ct, nm_lm);
  uintptr_t sample_idx_base = 0;
  uintptr_t cur_bits = lm_eligible[0];
  for (uint32_t lm_idx = 0; lm_idx != lm_ct; ++lm_idx) {
    const uint32_t sample_idx = BitIter1(lm_eligible, &sample_idx_base, &cur_bits);
    const uint32_t cur_geno = GetNyparrEntry(genovec, sample_idx);
    if (cur_geno != 3) {
      const uint32_t fss_idx = sample_lm_to_fss_idx[lm_idx];
      if (IsSet(nm_fss, fss_idx)) {
        if ((!lm_within2_founder) || (!IsSet(lm_within2_founder, sample_idx))) {
          qfam_w[lm_idx] = S_CAST(double, S_CAST(int32_t, cur_geno)) - 1.0 - qfam_b[fss_idx];
          continue;
        }
        // --qfam-parents: a founder parent's deviate is only meaningful when
        // the other parent is genotyped too.
        const uint32_t cur_start = fss_starts[fss_idx];
        const uint32_t partner_idx = (fss_contents[cur_start] == sample_idx)? fss_contents[cur_start + 1] : fss_contents[cur_start];
        if (GetNyparrEntry(genovec, partner_idx) != 3) {
          qfam_w[lm_idx] = S_CAST(double, S_CAST(int32_t, cur_geno)) - 1.0 - qfam_b[fss_idx];
          continue;
        }
      }
    }
    const double cur_pheno = pheno_d2[lm_idx];
    qt_sum -= cur_pheno;
    qt_ssq -= cur_pheno * cur_pheno;
    ClearBit(lm_idx, nm_lm);
  }
  *qt_sum_ptr = qt_sum;
  *qt_ssq_ptr = qt_ssq;
}

// --qfam/--qfam-parents only: the genotype sum of squares doesn't depend on
// the flips, and samples with W == 0 can be skipped entirely.
void QfamFlipPrecalc(uint32_t lm_ct, const double* qfam_w, const double* pheno_d2, uintptr_t* nm_lm, double* geno_sum_ptr, double* geno_ssq_ptr, double* qt_g_prod_ptr) {
  double geno_sum = 0.0;
  double geno_ssq = 0.0;
  double qt_g_prod = 0.0;
  uintptr_t lm_idx_base = 0;
  uintptr_t cur_bits = nm_lm[0];
  const uint32_t nm_ct = PopcountWords(nm_lm, BitCtToWordCt(lm_ct));
  for (uint32_t uii = 0; uii != nm_ct; ++uii) {
    const uint32_t lm_idx = BitIter1(nm_lm, &lm_idx_base, &cur_bits);
    const double cur_geno = qfam_w[lm_idx];
    if (fabs(cur_geno) < kQfamSmallEpsilon) {
      ClearBit(lm_idx, nm_lm);
    } else {
      geno_sum += cur_geno;
      geno_ssq += cur_geno * cur_geno;
      qt_g_prod += cur_geno * pheno_d2[lm_idx];
    }
  }
  *geno_sum_ptr = geno_sum * 0.5;
  *geno_ssq_ptr = geno_ssq;
  *qt_g_prod_ptr = qt_g_prod * 0.5;
}

// Univariate OLS of the phenotype on the (permuted, flipped) genotype score.
// Returns 1 when there is nothing to report.
BoolErr QfamRegress(uint32_t test_type, uint32_t nind, uint32_t lm_ct, const uint32_t* sample_lm_to_fss_idx, const uintptr_t* nm_lm, const double* pheno_d2, const double* qfam_b, const double* qfam_w, const uint32_t* qfam_permute, const uintptr_t* qfam_flip, double nind_recip, double qt_sum, double qt_ssq, double geno_sum, double geno_ssq, double qt_g_prod, double* beta_ptr, double* tstat_ptr) {
  if (nind < 3) {
    return 1;
  }
  const uint32_t lm_ctl = BitCtToWordCt(lm_ct);
  if (test_type & (kQfamWithin1 | kQfamWithin2)) {
    for (uint32_t widx = 0; widx != lm_ctl; ++widx) {
      uintptr_t cur_word = nm_lm[widx] & qfam_flip[widx];
      while (cur_word) {
        const uint32_t lm_idx = widx * kBitsPerWord + ctzw(cur_word);
        const double neg_w = -qfam_w[lm_idx];
        geno_sum += neg_w;
        qt_g_prod += neg_w * pheno_d2[lm_idx];
        cur_word &= cur_word - 1;
      }
    }
    geno_sum *= 2;
    qt_g_prod *= 2;
  } else {
    uintptr_t lm_idx_base = 0;
    uintptr_t cur_bits = nm_lm[0];
    for (uint32_t uii = 0; uii != nind; ++uii) {
      const uint32_t lm_idx = BitIter1(nm_lm, &lm_idx_base, &cur_bits);
      const uint32_t fss_idx = qfam_permute[sample_lm_to_fss_idx[lm_idx]];
      double cur_geno = qfam_b[fss_idx];
      if (test_type == kQfamTotal) {
        const double cur_w = qfam_w[lm_idx];
        if (IsSet(qfam_flip, fss_idx)) {
          cur_geno -= cur_w;
        } else {
          cur_geno += cur_w;
        }
      }
      geno_sum += cur_geno;
      geno_ssq += cur_geno * cur_geno;
      qt_g_prod += cur_geno * pheno_d2[lm_idx];
    }
  }
  const double qt_mean = qt_sum * nind_recip;
  const double geno_mean = geno_sum * nind_recip;
  const double dof_recip = 1.0 / u31tod(nind - 1);
  const double qt_var = (qt_ssq - qt_sum * qt_mean) * dof_recip;
  const double geno_var = (geno_ssq - geno_sum * geno_mean) * dof_recip;
  if (geno_var == 0.0) {
    return 1;
  }
  const double qt_g_covar = (qt_g_prod - qt_sum * geno_mean) * dof_recip;
  const double geno_var_recip = 1.0 / geno_var;
  const double beta = qt_g_covar * geno_var_recip;
  const double resid = qt_var * geno_var_recip - beta * beta;
  *beta_ptr = beta;
  *tstat_ptr = beta * sqrt(u31tod(nind - 2) / resid);
  return 0;
}

// Families (both parents present), then sibships (samples sharing FID/PAT/MAT
// whose parents aren't both present), then singletons.  Parents who appear in
// more than one family are dropped from the regression, since otherwise the
// result would depend on which of their families came first in the file.
PglErr QfamBuildGroups(const PedigreeIdInfo* piip, const uintptr_t* founder_collapsed, const FamilyInfo* fip, const uint32_t* idx_to_uidx, uint32_t sample_ct, uint32_t test_type, QfamGroups* qgp) {
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t family_ct = fip->family_ct;
    const uint32_t trio_ct = fip->trio_ct;
    const uint64_t* family_list = fip->family_list;
    const uint64_t* trio_list = fip->trio_list;
    uint32_t* family_child_offsets;
    uint32_t* family_children;
    uintptr_t* not_in_family;
    uintptr_t* is_child;
    uintptr_t* double_parent;
    uintptr_t* is_parent;
    uint32_t* sample_to_fss_idx;
    if (unlikely(bigstack_calloc_u32(family_ct + 1, &family_child_offsets) ||
                 bigstack_alloc_u32(trio_ct, &family_children) ||
                 bigstack_alloc_w(sample_ctl, &not_in_family) ||
                 bigstack_calloc_w(sample_ctl, &is_child) ||
                 bigstack_calloc_w(sample_ctl, &double_parent) ||
                 bigstack_calloc_w(sample_ctl, &is_parent) ||
                 bigstack_alloc_u32(sample_ct, &sample_to_fss_idx))) {
      goto QfamBuildGroups_ret_NOMEM;
    }
    SetAllBits(sample_ct, not_in_family);
    SetAllU32Arr(sample_ct, sample_to_fss_idx);
    for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
      family_child_offsets[1 + (trio_list[trio_idx] >> 32)] += 1;
    }
    for (uint32_t family_idx = 0; family_idx != family_ct; ++family_idx) {
      family_child_offsets[family_idx + 1] += family_child_offsets[family_idx];
    }
    {
      uint32_t* write_idxs;
      if (unlikely(bigstack_end_alloc_u32(family_ct + 1, &write_idxs))) {
        goto QfamBuildGroups_ret_NOMEM;
      }
      memcpy(write_idxs, family_child_offsets, (family_ct + 1) * sizeof(int32_t));
      for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
        const uint64_t trio_code = trio_list[trio_idx];
        const uint32_t child_idx = S_CAST(uint32_t, trio_code);
        family_children[write_idxs[trio_code >> 32]++] = child_idx;
        SetBit(child_idx, is_child);
        sample_to_fss_idx[child_idx] = trio_code >> 32;
      }
      BigstackEndReset(write_idxs);
    }
    for (uint32_t family_idx = 0; family_idx != family_ct; ++family_idx) {
      const uint64_t family_code = family_list[family_idx];
      const uint32_t parent_idxs[2] = {S_CAST(uint32_t, family_code), S_CAST(uint32_t, family_code >> 32)};
      for (uint32_t uii = 0; uii != 2; ++uii) {
        const uint32_t parent_idx = parent_idxs[uii];
        SetBit(parent_idx, is_parent);
        if (IsSet(not_in_family, parent_idx)) {
          if (sample_to_fss_idx[parent_idx] == UINT32_MAX) {
            sample_to_fss_idx[parent_idx] = family_idx;
          }
          ClearBit(parent_idx, not_in_family);
        } else {
          SetBit(parent_idx, double_parent);
        }
      }
    }
    BitvecInvmask(is_child, sample_ctl, not_in_family);
    BitvecInvmask(is_child, sample_ctl, double_parent);
    // Children have their family's index no matter what, and a parent keeps
    // the index of the first family it appeared in.
    for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
      const uint64_t trio_code = trio_list[trio_idx];
      sample_to_fss_idx[S_CAST(uint32_t, trio_code)] = trio_code >> 32;
    }

    // Upper bound: every family contributes 2 parents plus its children, and
    // every other sample appears in exactly one sibship or singleton entry.
    const uintptr_t fss_contents_capacity = 2 * S_CAST(uintptr_t, family_ct) + trio_ct + sample_ct;
    const uintptr_t fss_capacity = family_ct + sample_ct + 2;
    uint32_t* fss_starts;
    uint32_t* fss_contents;
    if (unlikely(bigstack_alloc_u32(fss_capacity, &fss_starts) ||
                 bigstack_alloc_u32(fss_contents_capacity, &fss_contents))) {
      goto QfamBuildGroups_ret_NOMEM;
    }
    uint32_t* contents_iter = fss_contents;
    for (uint32_t family_idx = 0; family_idx != family_ct; ++family_idx) {
      fss_starts[family_idx] = contents_iter - fss_contents;
      const uint64_t family_code = family_list[family_idx];
      *contents_iter++ = S_CAST(uint32_t, family_code);
      *contents_iter++ = family_code >> 32;
      const uint32_t child_end = family_child_offsets[family_idx + 1];
      for (uint32_t uii = family_child_offsets[family_idx]; uii != child_end; ++uii) {
        *contents_iter++ = family_children[uii];
      }
    }
    uint32_t fs_ct = family_ct;
    {
      // Sibships, keyed by FID/PAT/MAT.
      uint32_t candidate_ct = 0;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        if (IsSet(not_in_family, sample_idx) && (!IsSet(founder_collapsed, sample_idx))) {
          ++candidate_ct;
        }
      }
      if (candidate_ct > 1) {
        const SampleIdInfo* siip = &(piip->sii);
        const ParentalIdInfo* parental_id_infop = &(piip->parental_id_info);
        const uintptr_t max_key_blen = siip->max_sample_id_blen + parental_id_infop->max_paternal_id_blen + parental_id_infop->max_maternal_id_blen;
        char* key_strbox;
        uint32_t* key_id_map;
        if (unlikely(bigstack_end_alloc_c(max_key_blen * candidate_ct, &key_strbox) ||
                     bigstack_end_alloc_u32(candidate_ct, &key_id_map))) {
          goto QfamBuildGroups_ret_NOMEM;
        }
        uint32_t key_idx = 0;
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          if ((!IsSet(not_in_family, sample_idx)) || IsSet(founder_collapsed, sample_idx)) {
            continue;
          }
          const uint32_t sample_uidx = idx_to_uidx[sample_idx];
          const char* cur_sample_id = &(siip->sample_ids[sample_uidx * siip->max_sample_id_blen]);
          const char* fid_end = AdvToDelim(cur_sample_id, '\t');
          char* write_iter = &(key_strbox[key_idx * max_key_blen]);
          write_iter = memcpyax(write_iter, cur_sample_id, fid_end - cur_sample_id, '\t');
          write_iter = strcpyax(write_iter, &(parental_id_infop->paternal_ids[sample_uidx * parental_id_infop->max_paternal_id_blen]), '\t');
          // strcpya() doesn't write the terminator, and the keys are compared
          // as C strings below.
          *strcpya(write_iter, &(parental_id_infop->maternal_ids[sample_uidx * parental_id_infop->max_maternal_id_blen])) = '\0';
          key_id_map[key_idx] = sample_idx;
          ++key_idx;
        }
        if (unlikely(SortStrboxIndexed(candidate_ct, max_key_blen, 0, key_strbox, key_id_map))) {
          goto QfamBuildGroups_ret_NOMEM;
        }
        uint32_t run_start = 0;
        while (run_start != candidate_ct) {
          const char* cur_key = &(key_strbox[run_start * max_key_blen]);
          uint32_t run_end = run_start + 1;
          while ((run_end != candidate_ct) && (!strcmp(cur_key, &(key_strbox[run_end * max_key_blen])))) {
            ++run_end;
          }
          if (run_end - run_start >= 2) {
            fss_starts[fs_ct] = contents_iter - fss_contents;
            for (uint32_t uii = run_start; uii != run_end; ++uii) {
              const uint32_t sample_idx = key_id_map[uii];
              sample_to_fss_idx[sample_idx] = fs_ct;
              ClearBit(sample_idx, not_in_family);
              *contents_iter++ = sample_idx;
            }
            ++fs_ct;
          }
          run_start = run_end;
        }
        BigstackEndReset(key_strbox);
      }
    }
    uint32_t fss_ct = fs_ct;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      if (!IsSet(not_in_family, sample_idx)) {
        continue;
      }
      fss_starts[fss_ct] = contents_iter - fss_contents;
      sample_to_fss_idx[sample_idx] = fss_ct;
      *contents_iter++ = sample_idx;
      ++fss_ct;
    }
    fss_starts[fss_ct] = contents_iter - fss_contents;

    uintptr_t* lm_eligible;
    if (unlikely(bigstack_alloc_w(sample_ctl, &lm_eligible))) {
      goto QfamBuildGroups_ret_NOMEM;
    }
    SetAllBits(sample_ct, lm_eligible);
    BitvecInvmask(double_parent, sample_ctl, lm_eligible);
    if (test_type == kQfamWithin1) {
      BitvecInvmask(founder_collapsed, sample_ctl, lm_eligible);
    }
    const uint32_t lm_ct = PopcountWords(lm_eligible, sample_ctl);
    uint32_t* sample_lm_to_fss_idx;
    if (unlikely(bigstack_alloc_u32(lm_ct, &sample_lm_to_fss_idx))) {
      goto QfamBuildGroups_ret_NOMEM;
    }
    {
      uintptr_t sample_idx_base = 0;
      uintptr_t cur_bits = lm_eligible[0];
      for (uint32_t lm_idx = 0; lm_idx != lm_ct; ++lm_idx) {
        sample_lm_to_fss_idx[lm_idx] = sample_to_fss_idx[BitIter1(lm_eligible, &sample_idx_base, &cur_bits)];
      }
    }
    uintptr_t* lm_within2_founder = nullptr;
    if (test_type == kQfamWithin2) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &lm_within2_founder))) {
        goto QfamBuildGroups_ret_NOMEM;
      }
      memcpy(lm_within2_founder, is_parent, sample_ctl * sizeof(intptr_t));
      BitvecInvmask(double_parent, sample_ctl, lm_within2_founder);
      BitvecAnd(founder_collapsed, sample_ctl, lm_within2_founder);
    }
    qgp->fss_starts = fss_starts;
    qgp->fss_contents = fss_contents;
    qgp->sample_lm_to_fss_idx = sample_lm_to_fss_idx;
    qgp->lm_eligible = lm_eligible;
    qgp->lm_within2_founder = lm_within2_founder;
    qgp->family_ct = family_ct;
    qgp->fs_ct = fs_ct;
    qgp->fss_ct = fss_ct;
    qgp->lm_ct = lm_ct;
  }
  while (0) {
  QfamBuildGroups_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
  return reterr;
}

PglErr QfamReport(const uintptr_t* orig_sample_include, const PedigreeIdInfo* piip, const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const PermConfig* perm_config_ptr, uint32_t raw_sample_ct, uint32_t pheno_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, PgenGlobalFlags gflags, QfamFlags flags, uint32_t mperm_ct, uint32_t max_thread_ct, PgenReader* simple_pgrp, sfmt_t* sfmtp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PreinitCstream(&css);
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t test_type = flags & kfQfamTestMask;
    const uint32_t only_within = (test_type & (kQfamWithin1 | kQfamWithin2))? 1 : 0;
    const char* flag_suffix = (test_type == kQfamWithin1)? "within" : ((test_type == kQfamWithin2)? "parents" : ((test_type == kQfamTotal)? "total" : "between"));
    const char* test_str = only_within? "WITHIN" : ((test_type == kQfamTotal)? "TOTAL" : "BETWEEN");
    const PhenoCol* qt_pheno_col = nullptr;
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      if (pheno_cols[pheno_idx].type_code == kPhenoDtypeQt) {
        qt_pheno_col = &(pheno_cols[pheno_idx]);
        break;
      }
    }
    if (unlikely(!qt_pheno_col)) {
      logerrprintfww("Error: --qfam-%s requires a quantitative phenotype.\n", flag_suffix);
      goto QfamReport_ret_INCONSISTENT_INPUT;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sample_include))) {
      goto QfamReport_ret_NOMEM;
    }
    BitvecAndCopy(orig_sample_include, qt_pheno_col->nonmiss, raw_sample_ctl, sample_include);
    uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    if (unlikely(!sample_ct)) {
      logerrprintfww("Error: --qfam-%s requires at least one sample with a nonmissing quantitative phenotype.\n", flag_suffix);
      goto QfamReport_ret_INCONSISTENT_INPUT;
    }

    // As in PLINK 1.x, the haploid genome and chrMT are out of scope.
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* autosomal_variant_include;
    if (unlikely(bigstack_alloc_w(raw_variant_ctl, &autosomal_variant_include))) {
      goto QfamReport_ret_NOMEM;
    }
    memcpy(autosomal_variant_include, variant_include, raw_variant_ctl * sizeof(intptr_t));
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    uint32_t skipped_variant_ct = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != cip->chr_ct; ++chr_fo_idx) {
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if ((!IsSet(cip->haploid_mask, chr_idx)) && (chr_idx != mt_code)) {
        continue;
      }
      const uint32_t chr_vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      const uint32_t chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      skipped_variant_ct += PopcountBitRange(autosomal_variant_include, chr_vidx_start, chr_vidx_end);
      ClearBitsNz(chr_vidx_start, chr_vidx_end, autosomal_variant_include);
    }
    if (skipped_variant_ct) {
      logprintf("--qfam-%s: Excluding %u haploid/MT variant%s.\n", flag_suffix, skipped_variant_ct, (skipped_variant_ct == 1)? "" : "s");
    }
    const uint32_t autosomal_variant_ct = variant_ct - skipped_variant_ct;
    if (unlikely(!autosomal_variant_ct)) {
      logerrprintfww("Error: No variants remaining for --qfam-%s.\n", flag_suffix);
      goto QfamReport_ret_INCONSISTENT_INPUT;
    }

    FamilyInfo family_info;
    PreinitFamilyInfo(&family_info);
    reterr = GetTriosAndFamilies(sample_include, piip, founder_info, sex_nm, sex_male, raw_sample_ct, kfTrio0, &sample_ct, nullptr, &family_info);
    if (unlikely(reterr)) {
      goto QfamReport_ret_1;
    }
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* founder_collapsed;
    uint32_t* idx_to_uidx;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(sample_ctl, &founder_collapsed) ||
                 bigstack_alloc_u32(sample_ct, &idx_to_uidx))) {
      goto QfamReport_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    CopyBitarrSubset(founder_info, sample_include, sample_ct, founder_collapsed);
    ZeroTrailingBits(sample_ct, founder_collapsed);
    {
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        idx_to_uidx[sample_idx] = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      }
    }
    QfamGroups qgroups;
    reterr = QfamBuildGroups(piip, founder_collapsed, &family_info, idx_to_uidx, sample_ct, test_type, &qgroups);
    if (unlikely(reterr)) {
      goto QfamReport_ret_1;
    }
    const uint32_t fss_ct = qgroups.fss_ct;
    const uint32_t lm_ct = qgroups.lm_ct;
    if (unlikely(fss_ct < 2)) {
      logerrprintfww("Error: --qfam-%s requires at least two families.\n", flag_suffix);
      goto QfamReport_ret_INCONSISTENT_INPUT;
    }
    if (unlikely(lm_ct < 3)) {
      logerrprintfww("Error: Fewer than three eligible %ss for --qfam-%s.\n", (test_type == kQfamWithin1)? "nonfounder" : "sample", flag_suffix);
      goto QfamReport_ret_INCONSISTENT_INPUT;
    }
    const uint32_t fss_ctl = BitCtToWordCt(fss_ct);
    const uint32_t lm_ctl = BitCtToWordCt(lm_ct);
    const uint32_t flip_ctl = only_within? lm_ctl : fss_ctl;

    double* pheno_d2;
    if (unlikely(bigstack_alloc_d(lm_ct, &pheno_d2))) {
      goto QfamReport_ret_NOMEM;
    }
    double qt_sum_all = 0.0;
    double qt_ssq_all = 0.0;
    {
      const double* pheno_qt = qt_pheno_col->data.qt;
      uintptr_t sample_idx_base = 0;
      uintptr_t cur_bits = qgroups.lm_eligible[0];
      for (uint32_t lm_idx = 0; lm_idx != lm_ct; ++lm_idx) {
        const uint32_t sample_idx = BitIter1(qgroups.lm_eligible, &sample_idx_base, &cur_bits);
        const double cur_pheno = pheno_qt[idx_to_uidx[sample_idx]];
        pheno_d2[lm_idx] = cur_pheno;
        qt_sum_all += cur_pheno;
        qt_ssq_all += cur_pheno * cur_pheno;
      }
    }

    const uint32_t perm_adapt = !mperm_ct;
    const uint32_t perms_total = perm_adapt? perm_config_ptr->aperm_max : mperm_ct;
    const double aperm_alpha = perm_config_ptr->aperm_alpha;
    const double adaptive_intercept = perm_config_ptr->aperm_init_interval;
    const double adaptive_slope = perm_config_ptr->aperm_interval_slope;
    const double adaptive_ci_zt = perm_adapt? QuantileToZscore(1.0 - perm_config_ptr->aperm_beta / (2.0 * u31tod(autosomal_variant_ct))) : 0.0;
    uint32_t first_adapt_check = perm_config_ptr->aperm_min;
    if (u31tod(first_adapt_check) < adaptive_intercept) {
      first_adapt_check = S_CAST(uint32_t, adaptive_intercept);
    }

    const uint32_t emp_se = (flags / kfQfamEmpSe) & 1;
    const uint32_t perm_count = (flags / kfQfamPermCount) & 1;
    const uint32_t output_zst = (flags / kfQfamZs) & 1;
    const uint32_t chr_col = flags & kfQfamColChrom;
    const uint32_t ref_col = flags & kfQfamColRef;
    const uint32_t alt1_col = flags & kfQfamColAlt1;
    const uint32_t alt_col = flags & kfQfamColAlt;
    const uint32_t a1_col = flags & kfQfamColA1;
    const uint32_t all_nonref = (gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    uint32_t provref_col = 0;
    if (ref_col) {
      if (flags & kfQfamColProvref) {
        provref_col = 1;
      } else if (flags & kfQfamColMaybeprovref) {
        provref_col = all_nonref || (nonref_flags && (!IntersectionRangeIsEmpty(autosomal_variant_include, nonref_flags, 0, raw_variant_ct)));
      }
    }
    const uint32_t test_col = flags & kfQfamColTest;
    const uint32_t nind_col = flags & kfQfamColNind;
    const uint32_t beta_col = flags & kfQfamColBeta;
    const uint32_t stat_col = flags & kfQfamColStat;
    const uint32_t rawp_col = flags & kfQfamColRawp;
    const uint32_t emp1_col = flags & kfQfamColEmp1;
    const uint32_t np_col = flags & kfQfamColNp;

    uintptr_t* nm_fss;
    uintptr_t* nm_lm;
    uintptr_t* qfam_flip;
    uintptr_t* fss_flip;
    double* qfam_b;
    double* qfam_w;
    uint32_t* dummy_perm;
    uint32_t* cur_perm;
    uint32_t* permute_edit_buf;
    uintptr_t* genovec;
    uintptr_t* genovec_buf;
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    if (unlikely(bigstack_alloc_w(fss_ctl, &nm_fss) ||
                 bigstack_alloc_w(lm_ctl, &nm_lm) ||
                 bigstack_calloc_w(flip_ctl, &qfam_flip) ||
                 bigstack_alloc_w(fss_ctl, &fss_flip) ||
                 bigstack_alloc_d(fss_ct, &qfam_b) ||
                 bigstack_alloc_d(lm_ct, &qfam_w) ||
                 bigstack_alloc_u32(fss_ct, &dummy_perm) ||
                 bigstack_alloc_u32(fss_ct, &cur_perm) ||
                 bigstack_alloc_u32(fss_ct, &permute_edit_buf) ||
                 bigstack_alloc_w(sample_ctl2, &genovec) ||
                 bigstack_alloc_w(sample_ctl2, &genovec_buf))) {
      goto QfamReport_ret_NOMEM;
    }
    for (uint32_t fss_idx = 0; fss_idx != fss_ct; ++fss_idx) {
      dummy_perm[fss_idx] = fss_idx;
    }
    uintptr_t* dummy_flip;
    if (unlikely(bigstack_calloc_w(flip_ctl, &dummy_flip))) {
      goto QfamReport_ret_NOMEM;
    }

    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    char* chr_buf;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto QfamReport_ret_NOMEM;
    }
    char* outname_end2 = strcpya_k(outname_end, ".qfam.");
    outname_end2 = strcpya(outname_end2, flag_suffix);
    OutnameZstSet("", output_zst, outname_end2);
    const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 512 + 2 * max_allele_slen + max_chr_blen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, MAXV(max_thread_ct - 1, 1), overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto QfamReport_ret_1;
    }
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (flags & kfQfamColPos) {
      cswritep = strcpya_k(cswritep, "POS\t");
    } else {
      variant_bps = nullptr;
    }
    cswritep = strcpya_k(cswritep, "ID");
    if (ref_col) {
      cswritep = strcpya_k(cswritep, "\tREF");
    }
    if (alt1_col) {
      cswritep = strcpya_k(cswritep, "\tALT1");
    }
    if (alt_col) {
      cswritep = strcpya_k(cswritep, "\tALT");
    }
    if (provref_col) {
      cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
    }
    if (a1_col) {
      cswritep = strcpya_k(cswritep, "\tA1");
    }
    if (test_col) {
      cswritep = strcpya_k(cswritep, "\tTEST");
    }
    if (nind_col) {
      cswritep = strcpya_k(cswritep, "\tNIND");
    }
    if (beta_col) {
      cswritep = strcpya_k(cswritep, "\tBETA");
    }
    if (stat_col) {
      cswritep = strcpya_k(cswritep, "\tSTAT");
    }
    if (rawp_col) {
      cswritep = strcpya_k(cswritep, "\tRAW_P");
    }
    if (emp_se) {
      cswritep = strcpya_k(cswritep, "\tEMP_BETA\tEMP_SE");
    }
    if (emp1_col) {
      cswritep = strcpya_k(cswritep, "\tEMP1");
    }
    if (np_col) {
      cswritep = strcpya_k(cswritep, "\tNP");
    }
    AppendBinaryEoln(&cswritep);
    logprintfww("--qfam-%s: Permuting %u families/sibships/singletons, %u sample%s in the regression, up to %u permutation%s.\n", flag_suffix, fss_ct, lm_ct, (lm_ct == 1)? "" : "s", perms_total, (perms_total == 1)? "" : "s");

    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    const uint32_t* trio_lookup = family_info.trio_lookup;
    const uint32_t trio_ct = family_info.trio_ct;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = autosomal_variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_blen = 0;
    uint32_t regress_fail_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (autosomal_variant_ct + 99) / 100;
    printf("--qfam-%s: 0%%", flag_suffix);
    fflush(stdout);
    for (uint32_t variant_idx = 0; variant_idx != autosomal_variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(autosomal_variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        char* chr_name_end = chrtoa(cip, cip->chr_file_order[chr_fo_idx], chr_buf);
        *chr_name_end++ = '\t';
        chr_blen = chr_name_end - chr_buf;
      }
      reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, genovec);
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, variant_uidx);
        goto QfamReport_ret_1;
      }
      ZeroTrailingNyps(sample_ct, genovec);
      // Mendel errors are set to missing first, as in PLINK 1.x, reading every
      // trio from the unmodified genotypes.
      {
        memcpy(genovec_buf, genovec, sample_ctl2 * sizeof(intptr_t));
        const uint32_t* trio_lookup_iter = trio_lookup;
        for (uint32_t trio_idx = 0; trio_idx != trio_ct; ++trio_idx) {
          const uint32_t child_idx = trio_lookup_iter[0];
          const uint32_t dad_idx = trio_lookup_iter[1];
          const uint32_t mom_idx = trio_lookup_iter[2];
          trio_lookup_iter = &(trio_lookup_iter[3]);
          const uint32_t child_geno = GetNyparrEntry(genovec_buf, child_idx);
          if (child_geno == 3) {
            continue;
          }
          const uint32_t error_result = kBiallelicMendelErrorTableAutosomalOrX[GetNyparrEntry(genovec_buf, dad_idx) + GetNyparrEntry(genovec_buf, mom_idx) * 4 + child_geno * 16];
          if (!error_result) {
            continue;
          }
          SetNyparrEntryTo3(child_idx, genovec);
          if (error_result & 0x100) {
            SetNyparrEntryTo3(dad_idx, genovec);
          }
          if (error_result & 0x10000) {
            SetNyparrEntryTo3(mom_idx, genovec);
          }
        }
      }
      double qt_sum;
      double qt_ssq;
      QfamComputeBw(genovec, &qgroups, pheno_d2, qt_sum_all, qt_ssq_all, nm_fss, nm_lm, qfam_b, qfam_w, &qt_sum, &qt_ssq);
      const uint32_t cur_fss_ct = PopcountWords(nm_fss, fss_ctl);
      const uint32_t nind = PopcountWords(nm_lm, lm_ctl);
      const double nind_recip = 1.0 / u31tod(nind);
      double geno_sum = 0.0;
      double geno_ssq = 0.0;
      double qt_g_prod = 0.0;
      if (only_within) {
        QfamFlipPrecalc(lm_ct, qfam_w, pheno_d2, nm_lm, &geno_sum, &geno_ssq, &qt_g_prod);
      }
      double beta = 0.0;
      double tstat = 0.0;
      const uint32_t regress_ok = !QfamRegress(test_type, nind, lm_ct, qgroups.sample_lm_to_fss_idx, nm_lm, pheno_d2, qfam_b, qfam_w, dummy_perm, dummy_flip, nind_recip, qt_sum, qt_ssq, geno_sum, geno_ssq, qt_g_prod, &beta, &tstat);
      uint32_t success_2 = 0;
      uint32_t perm_attempt_ct = perms_total;
      double beta_sum = 0.0;
      double beta_ssq = 0.0;
      uint32_t beta_fail_ct = 0;
      if (!regress_ok) {
        ++regress_fail_ct;
      } else {
        const double stat_high = fabs(tstat) + kQfamEpsilon;
        const double stat_low = fabs(tstat) - kQfamEpsilon;
        uint32_t next_adapt_check = first_adapt_check;
        for (uint32_t perm_idx = 0; perm_idx != perms_total; ++perm_idx) {
          // One random sign flip per family, and (for the tests that use B) a
          // random permutation of the families.
          for (uint32_t widx = 0; widx != fss_ctl; ++widx) {
            uintptr_t cur_word = sfmt_genrand_uint32(sfmtp);
#ifdef __LP64__
            cur_word |= S_CAST(uintptr_t, sfmt_genrand_uint32(sfmtp)) << 32;
#endif
            fss_flip[widx] = cur_word;
          }
          const uint32_t* perm_ptr = dummy_perm;
          if (only_within) {
            ZeroWArr(lm_ctl, qfam_flip);
            for (uint32_t lm_idx = 0; lm_idx != lm_ct; ++lm_idx) {
              if (IsSet(fss_flip, qgroups.sample_lm_to_fss_idx[lm_idx])) {
                SetBit(lm_idx, qfam_flip);
              }
            }
          } else {
            memcpy(qfam_flip, fss_flip, fss_ctl * sizeof(intptr_t));
            for (uint32_t fss_idx = 0; fss_idx != fss_ct; ++fss_idx) {
              cur_perm[fss_idx] = fss_idx;
            }
            for (uint32_t uii = fss_ct - 1; uii; --uii) {
              const uint32_t ujj = RandU32(uii + 1, sfmtp);
              const uint32_t tmp = cur_perm[uii];
              cur_perm[uii] = cur_perm[ujj];
              cur_perm[ujj] = tmp;
            }
            perm_ptr = cur_perm;
            if (cur_fss_ct != fss_ct) {
              // A family with no genotype must not be permuted in; walk the
              // permutation cycle until a genotyped one is found, as PLINK
              // 1.07's qfam.cpp does.
              memcpy(permute_edit_buf, cur_perm, fss_ct * sizeof(int32_t));
              uintptr_t fss_idx_base = 0;
              uintptr_t nm_bits = nm_fss[0];
              for (uint32_t uii = 0; uii != cur_fss_ct; ++uii) {
                const uint32_t orig_fss_idx = BitIter1(nm_fss, &fss_idx_base, &nm_bits);
                uint32_t new_fss_idx = permute_edit_buf[orig_fss_idx];
                if (IsSet(nm_fss, new_fss_idx)) {
                  continue;
                }
                uint32_t ujj;
                while (1) {
                  ujj = permute_edit_buf[new_fss_idx];
                  permute_edit_buf[new_fss_idx] = new_fss_idx;
                  if (IsSet(nm_fss, ujj)) {
                    break;
                  }
                  new_fss_idx = ujj;
                }
                permute_edit_buf[orig_fss_idx] = ujj;
              }
              perm_ptr = permute_edit_buf;
            }
          }
          double perm_beta;
          double perm_tstat;
          if (!QfamRegress(test_type, nind, lm_ct, qgroups.sample_lm_to_fss_idx, nm_lm, pheno_d2, qfam_b, qfam_w, perm_ptr, qfam_flip, nind_recip, qt_sum, qt_ssq, geno_sum, geno_ssq, qt_g_prod, &perm_beta, &perm_tstat)) {
            beta_sum += perm_beta;
            beta_ssq += perm_beta * perm_beta;
            const double perm_abs_tstat = fabs(perm_tstat);
            if (perm_abs_tstat > stat_high) {
              success_2 += 2;
            } else if (perm_abs_tstat > stat_low) {
              success_2 += 1;
            }
          } else {
            // Conservative handling of a permutation that can't be fitted.
            success_2 += 2;
            ++beta_fail_ct;
          }
          if (perm_adapt && (perm_idx + 1 == next_adapt_check)) {
            if (success_2) {
              const double pval = S_CAST(double, S_CAST(int32_t, success_2 + 2)) / (2 * u31tod(next_adapt_check + 1));
              const double ci_halfwidth = adaptive_ci_zt * sqrt(pval * (1 - pval) / u31tod(next_adapt_check));
              if ((pval - ci_halfwidth > aperm_alpha) || (pval + ci_halfwidth < aperm_alpha)) {
                perm_attempt_ct = next_adapt_check;
                break;
              }
            }
            next_adapt_check += S_CAST(int32_t, adaptive_intercept + u31tod(next_adapt_check) * adaptive_slope);
          }
        }
      }

      if (chr_col) {
        cswritep = memcpya(cswritep, chr_buf, chr_blen);
      }
      if (variant_bps) {
        cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      }
      cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
      const uintptr_t allele_idx_offset_base = allele_idx_offsets? allele_idx_offsets[variant_uidx] : (2 * variant_uidx);
      const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base) : 2;
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (ref_col) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[0]);
      }
      if (alt1_col) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[1]);
      }
      if (alt_col) {
        *cswritep++ = '\t';
        for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto QfamReport_ret_WRITE_FAIL;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
        }
        --cswritep;
      }
      if (provref_col) {
        *cswritep++ = '\t';
        *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
      }
      if (a1_col) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[1]);
      }
      if (test_col) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, test_str);
      }
      if (nind_col) {
        *cswritep++ = '\t';
        cswritep = u32toa(nind, cswritep);
      }
      if (!regress_ok) {
        if (beta_col) {
          cswritep = strcpya_k(cswritep, "\tNA");
        }
        if (stat_col) {
          cswritep = strcpya_k(cswritep, "\tNA");
        }
        if (rawp_col) {
          cswritep = strcpya_k(cswritep, "\tNA");
        }
        if (emp_se) {
          cswritep = strcpya_k(cswritep, "\tNA\tNA");
        }
        if (emp1_col) {
          cswritep = strcpya_k(cswritep, "\tNA");
        }
        if (np_col) {
          cswritep = strcpya_k(cswritep, "\tNA");
        }
      } else {
        if (beta_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(beta, cswritep);
        }
        if (stat_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(tstat, cswritep);
        }
        if (rawp_col) {
          // Deliberately not passed through --output-min-p: only the empirical
          // p-value is meant to be interpreted here.
          *cswritep++ = '\t';
          cswritep = lntoa_g(TstatToLnP(tstat, nind - 2), cswritep);
        }
        if (emp_se) {
          const uint32_t beta_ct = perm_attempt_ct - beta_fail_ct;
          if (beta_ct <= 1) {
            cswritep = strcpya_k(cswritep, "\tNA\tNA");
          } else {
            const double beta_mean = beta_sum / u31tod(beta_ct);
            *cswritep++ = '\t';
            cswritep = dtoa_g(beta_mean, cswritep);
            *cswritep++ = '\t';
            cswritep = dtoa_g(sqrt((beta_ssq - beta_sum * beta_mean) / u31tod(beta_ct - 1)), cswritep);
          }
        }
        if (emp1_col) {
          *cswritep++ = '\t';
          if (perm_count) {
            cswritep = dtoa_g(S_CAST(double, S_CAST(int32_t, success_2)) * 0.5, cswritep);
          } else {
            cswritep = dtoa_g(S_CAST(double, S_CAST(int32_t, success_2 + 2)) / (2 * u31tod(perm_attempt_ct + 1)), cswritep);
          }
        }
        if (np_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(perm_attempt_ct, cswritep);
        }
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto QfamReport_ret_WRITE_FAIL;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / autosomal_variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, autosomal_variant_ct)) / 100;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto QfamReport_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    if (regress_fail_ct) {
      logprintf("--qfam-%s: %u variant%s had no usable regression.\n", flag_suffix, regress_fail_ct, (regress_fail_ct == 1)? "" : "s");
    }
    logprintfww("--qfam-%s report written to %s .\n", flag_suffix, outname);
  }
  while (0) {
  QfamReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  QfamReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  QfamReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 QfamReport_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
