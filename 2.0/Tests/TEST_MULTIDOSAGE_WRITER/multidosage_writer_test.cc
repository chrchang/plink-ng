// Unit test for SpgwAppendMultiallelicGenovecDphase16(): writes .pgen files
// with multiallelic dosage records (rare-alt dosage tracks #5-6, #9-10) from
// random but valid data, then checks each record's type byte and that the
// record ends with independently computed bytes for tracks #5-10.  The
// header is parsed by hand rather than with pgenlib's reader.

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "pgenlib_write.h"

using namespace plink2;

static int g_fail_ct = 0;
#define CHECK(cond) do { if (!(cond)) { fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #cond); ++g_fail_ct; } } while (0)

static uint64_t g_rng = 0x2545f4914f6cdd1dLLU;
static uint64_t NextRand() {
  uint64_t z = (g_rng += 0x9e3779b97f4a7c15LLU);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9LLU;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebLLU;
  return z ^ (z >> 31);
}
static uint32_t RandBelow(uint32_t n) { return NextRand() % n; }

// Zeroed, cacheline-aligned word buffer with room for vector-width reads.
static uintptr_t* AllocWords(uint32_t bit_ct) {
  const uintptr_t byte_ct = RoundUpPow2(DivUp(bit_ct, 8) + 64, 64);
  void* ptr = aligned_alloc(64, byte_ct);
  memset(ptr, 0, byte_ct);
  return static_cast<uintptr_t*>(ptr);
}

struct VariantData {
  uint32_t allele_ct;
  std::vector<uintptr_t*> bufs;  // owned
  PgenVariant pgv;
  // storage
  std::vector<AlleleCode> patch_01_vals, patch_10_vals, md_codes, mdp_codes;
  std::vector<uint16_t> dosage_main, md_vals;
  std::vector<int16_t> dphase_delta, mdp_delta;
  std::vector<unsigned char> md_cts, mdp_cts;
  // expected track bytes
  std::vector<unsigned char> expected_suffix;
  unsigned char expected_vrtype_bits;  // bits 3, 5-7 (bit 4 checked separately)
  uint32_t hphase;
};

static uintptr_t* NewBits(VariantData& vd, uint32_t bit_ct) {
  uintptr_t* ptr = AllocWords(bit_ct);
  vd.bufs.push_back(ptr);
  return ptr;
}

// 64-bit deltalist plus 16-bit values, encoded independently of the writer's
// multidosage code path (Deltalist64Append() has its own unit test).
static void AppendTrackPair(const std::vector<uint64_t>& idxs, const std::vector<uint16_t>& vals, uint64_t idx_end, std::vector<unsigned char>& out) {
  std::vector<unsigned char> buf(Deltalist64MaxByteCt(idxs.size(), Deltalist64IdxByteCt(idx_end)) + 16);
  unsigned char* end = Deltalist64Append(idxs.data(), idxs.size(), Deltalist64IdxByteCt(idx_end), buf.data());
  out.insert(out.end(), buf.data(), end);
  for (uint16_t v : vals) {
    out.push_back(v & 255);
    out.push_back(v >> 8);
  }
}

// dosage_mode: 0 = few dosages (list), 1 = about half (bitarray),
// 2 = all (unconditional, some missing).  rare_density: 0 = none,
// 1 = sparse, 2 = dense.
static void MakeVariant(uint32_t sample_ct, uint32_t allele_ct, uint32_t dosage_mode, uint32_t dphase, uint32_t hphase, uint32_t rare_density, VariantData& vd) {
  vd.allele_ct = allele_ct;
  vd.hphase = hphase;
  memset(&vd.pgv, 0, sizeof(PgenVariant));
  PgenVariant& pgv = vd.pgv;
  const uint32_t rare_ct = allele_ct - 2;
  pgv.genovec = NewBits(vd, 2 * sample_ct);
  pgv.patch_01_set = NewBits(vd, sample_ct);
  pgv.patch_10_set = NewBits(vd, sample_ct);
  pgv.dosage_present = NewBits(vd, sample_ct);
  pgv.dphase_present = NewBits(vd, sample_ct);
  pgv.multidosage_present = NewBits(vd, sample_ct);
  pgv.multidphase_present = NewBits(vd, sample_ct);
  pgv.phaseinfo = NewBits(vd, sample_ct);
  std::vector<uint32_t> cats(sample_ct);
  uint32_t het_ct = 0;
  for (uint32_t s = 0; s != sample_ct; ++s) {
    cats[s] = RandBelow(4);
  }
  // Unconditional mode: missing dosages need missing hardcalls.
  std::vector<uint32_t> dosage_missing(sample_ct, 0);
  std::vector<uint32_t> has_dosage(sample_ct, 0);
  for (uint32_t s = 0; s != sample_ct; ++s) {
    if (dosage_mode == 2) {
      has_dosage[s] = 1;
      if (RandBelow(10) == 0) {
        dosage_missing[s] = 1;
        cats[s] = 3;
      }
    } else if (dosage_mode == 1) {
      has_dosage[s] = RandBelow(2);
    } else {
      has_dosage[s] = (RandBelow(40) == 0);
    }
  }
  for (uint32_t s = 0; s != sample_ct; ++s) {
    pgv.genovec[s / kBitsPerWordD2] |= S_CAST(uintptr_t, cats[s]) << (2 * (s % kBitsPerWordD2));
    if ((cats[s] == 1) && RandBelow(3) == 0) {
      SetBit(s, pgv.patch_01_set);
      vd.patch_01_vals.push_back(2 + RandBelow(rare_ct));
    } else if ((cats[s] == 2) && RandBelow(3) == 0) {
      SetBit(s, pgv.patch_10_set);
      uint32_t lo = 1 + RandBelow(allele_ct - 1);
      uint32_t hi = 1 + RandBelow(allele_ct - 1);
      if (lo > hi) {
        std::swap(lo, hi);
      }
      if (hi == 1) {
        hi = 2;
      }
      vd.patch_10_vals.push_back(lo);
      vd.patch_10_vals.push_back(hi);
      het_ct += (lo != hi);
    }
    het_ct += (cats[s] == 1);
    if (RandBelow(2)) {
      SetBit(s, pgv.phaseinfo);
    }
  }
  pgv.patch_01_vals = vd.patch_01_vals.data();
  pgv.patch_10_vals = vd.patch_10_vals.data();
  pgv.patch_01_ct = vd.patch_01_vals.size();
  pgv.patch_10_ct = vd.patch_10_vals.size() / 2;
  if (hphase && het_ct) {
    pgv.phasepresent = nullptr;
    pgv.phasepresent_ct = het_ct;
  }
  // Dosages and rare-alt dosages.
  const uint64_t idx_end = S_CAST(uint64_t, sample_ct) << RareAltIdxBitCt(allele_ct);
  const uint32_t rare_alt_bit_ct = RareAltIdxBitCt(allele_ct);
  std::vector<uint64_t> md_idxs, mdp_idxs;
  std::vector<uint16_t> mdp_vals_u;
  std::vector<uint32_t> sample_dosage(sample_ct, 0);
  std::vector<uint16_t> dphase_bits;  // per dosage entry: 1 if explicit dphase
  for (uint32_t s = 0; s != sample_ct; ++s) {
    if (!has_dosage[s]) {
      continue;
    }
    SetBit(s, pgv.dosage_present);
    if (dosage_missing[s]) {
      vd.dosage_main.push_back(65535);
      if (dphase && (dosage_mode == 2)) {
        SetBit(s, pgv.dphase_present);
        vd.dphase_delta.push_back(-32768);
      }
      dphase_bits.push_back(dphase && (dosage_mode == 2));
      continue;
    }
    const uint32_t v = RandBelow(32769);
    vd.dosage_main.push_back(v);
    // rare-alt dosages
    std::vector<uint32_t> codes;
    std::vector<uint32_t> vals;
    if (rare_density) {
      const uint32_t max_entries = (rare_density == 1)? (RandBelow(4) == 0) : MINV(rare_ct, 255);
      uint32_t remaining = v;
      for (uint32_t code = 2; (code < allele_ct) && (codes.size() < max_entries) && remaining; ++code) {
        if ((rare_density == 2) || RandBelow(rare_ct) == 0) {
          const uint32_t cur = 1 + RandBelow(MINV(remaining, 1 + (v / 4)));
          codes.push_back(code);
          vals.push_back(cur);
          remaining -= cur;
        }
      }
    }
    if (!codes.empty()) {
      SetBit(s, pgv.multidosage_present);
      vd.md_cts.push_back(codes.size());
      for (size_t i = 0; i != codes.size(); ++i) {
        vd.md_codes.push_back(codes[i]);
        vd.md_vals.push_back(vals[i]);
        md_idxs.push_back((S_CAST(uint64_t, s) << rare_alt_bit_ct) | (codes[i] - 2));
      }
    }
    // dosage phase
    uint32_t explicit_dphase = 0;
    if (dphase && ((dosage_mode == 2) || RandBelow(2))) {
      explicit_dphase = 1;
      const int32_t lim = MINV(v, 32768 - v);
      const int32_t delta = lim? (S_CAST(int32_t, RandBelow(2 * lim + 1)) - lim) : 0;
      SetBit(s, pgv.dphase_present);
      vd.dphase_delta.push_back(delta);
      // rare-alt deltas, keeping each haplotype's alt1 share nonnegative
      int64_t sum_plus = 0;
      int64_t sum_minus = 0;
      std::vector<uint32_t> pcodes;
      std::vector<int32_t> pdeltas;
      for (size_t i = 0; i != codes.size(); ++i) {
        if (RandBelow(2)) {
          continue;
        }
        const int32_t vr = vals[i];
        const int32_t dr = S_CAST(int32_t, RandBelow(2 * vr + 1)) - vr;
        if ((!dr) || (sum_plus + vr + dr > v + delta) || (sum_minus + vr - dr > S_CAST(int32_t, v) - delta)) {
          sum_plus += vr;
          sum_minus += vr;
          continue;
        }
        sum_plus += vr + dr;
        sum_minus += vr - dr;
        pcodes.push_back(codes[i]);
        pdeltas.push_back(dr);
      }
      if (!pcodes.empty()) {
        SetBit(s, pgv.multidphase_present);
        vd.mdp_cts.push_back(pcodes.size());
        for (size_t i = 0; i != pcodes.size(); ++i) {
          vd.mdp_codes.push_back(pcodes[i]);
          vd.mdp_delta.push_back(pdeltas[i]);
          mdp_idxs.push_back((S_CAST(uint64_t, s) << rare_alt_bit_ct) | (pcodes[i] - 2));
          mdp_vals_u.push_back(S_CAST(uint16_t, pdeltas[i]));
        }
      }
    }
    dphase_bits.push_back(explicit_dphase);
  }
  pgv.dosage_main = vd.dosage_main.data();
  pgv.dosage_ct = vd.dosage_main.size();
  pgv.dphase_delta = vd.dphase_delta.data();
  pgv.dphase_ct = vd.dphase_delta.size();
  pgv.multidosage_cts = vd.md_cts.data();
  pgv.multidosage_codes = vd.md_codes.data();
  pgv.multidosage_vals = vd.md_vals.data();
  pgv.multidosage_sample_ct = vd.md_cts.size();
  pgv.multidphase_cts = vd.mdp_cts.data();
  pgv.multidphase_codes = vd.mdp_codes.data();
  pgv.multidphase_delta = vd.mdp_delta.data();
  pgv.multidphase_sample_ct = vd.mdp_cts.size();

  // Expected record suffix (tracks #5-10) and type bits.
  vd.expected_vrtype_bits = (pgv.patch_01_ct || pgv.patch_10_ct)? 8 : 0;
  vd.expected_suffix.clear();
  if (pgv.dosage_ct) {
    const uint32_t dosage_ct = pgv.dosage_ct;
    const uint32_t dphase_ct = pgv.dphase_ct;
    if (dosage_ct < sample_ct / kPglMaxDeltalistLenDivisor) {
      vd.expected_vrtype_bits |= 0x20;
    } else if ((dosage_ct == sample_ct) && ((!dphase_ct) || (dphase_ct == sample_ct))) {
      vd.expected_vrtype_bits |= 0x40;
    } else {
      vd.expected_vrtype_bits |= 0x60;
    }
    AppendTrackPair(md_idxs, vd.md_vals, idx_end, vd.expected_suffix);
    if (dphase_ct) {
      vd.expected_vrtype_bits |= 0x80;
      if (dphase_ct != sample_ct) {
        // track #7: bitarray over dosage entries
        for (uint32_t i = 0; i < dosage_ct; i += 8) {
          uint32_t byte = 0;
          for (uint32_t j = 0; (j != 8) && (i + j < dosage_ct); ++j) {
            byte |= dphase_bits[i + j] << j;
          }
          vd.expected_suffix.push_back(byte);
        }
      }
      for (int16_t d : vd.dphase_delta) {
        vd.expected_suffix.push_back(S_CAST(uint16_t, d) & 255);
        vd.expected_suffix.push_back(S_CAST(uint16_t, d) >> 8);
      }
      AppendTrackPair(mdp_idxs, mdp_vals_u, idx_end, vd.expected_suffix);
    }
  }
}

static void FreeVariant(VariantData& vd) {
  for (uintptr_t* p : vd.bufs) {
    free(p);
  }
  vd.bufs.clear();
}

static std::vector<unsigned char> ReadFile(const char* fname) {
  std::vector<unsigned char> out;
  FILE* ff = fopen(fname, "rb");
  if (!ff) {
    return out;
  }
  unsigned char buf[65536];
  size_t n;
  while ((n = fread(buf, 1, sizeof(buf), ff)) > 0) {
    out.insert(out.end(), buf, buf + n);
  }
  fclose(ff);
  return out;
}

static PglErr OpenWriter(const char* fname, const uintptr_t* allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags gflags, STPgenWriter* spgwp, unsigned char** alloc_ptr) {
  PreinitSpgw(spgwp);
  uintptr_t alloc_cacheline_ct;
  uint32_t max_vrec_len;
  PglErr reterr = SpgwInitPhase1(fname, allele_idx_offsets, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, gflags, 1, spgwp, &alloc_cacheline_ct, &max_vrec_len);
  if (reterr) {
    return reterr;
  }
  // kCacheline is 128 on e.g. Apple silicon
  *alloc_ptr = static_cast<unsigned char*>(aligned_alloc(kCacheline, alloc_cacheline_ct * kCacheline));
  SpgwInitPhase2(max_vrec_len, spgwp, *alloc_ptr);
  return kPglRetSuccess;
}

// Writes a file mixing biallelic records (existing API) and multiallelic
// dosage records, then checks it.
// combo_offset selects which (dosage form, dphase, hphase, rare density)
// combinations the multiallelic records cycle through.
static void RunFile(const char* fname, uint32_t sample_ct, uint32_t variant_ct, uint32_t combo_offset) {
  const uint32_t allele_cts[] = {3, 4, 5, 20, 255};
  std::vector<VariantData> vds(variant_ct);
  std::vector<uintptr_t> allele_idx_offsets(variant_ct + 1);
  std::vector<uint32_t> is_multi(variant_ct);
  // Each group of 6 variants: one per multiallelic allele count, then one
  // biallelic.
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    is_multi[vidx] = (vidx % 6 != 5);
    const uint32_t allele_ct = is_multi[vidx]? allele_cts[vidx % 6] : 2;
    allele_idx_offsets[vidx + 1] = allele_idx_offsets[vidx] + allele_ct;
    const uint32_t combo = combo_offset + vidx / 6;
    MakeVariant(sample_ct, MAXV(allele_ct, 3), combo % 3, (combo / 3) % 2, (combo / 6) % 2, (combo / 12) % 3, vds[vidx]);
    vds[vidx].allele_ct = allele_ct;
  }
  STPgenWriter spgw;
  unsigned char* alloc;
  const PgenGlobalFlags gflags = kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent | kfPgenGlobalMultiallelicDosagePresent;
  CHECK(!OpenWriter(fname, allele_idx_offsets.data(), variant_ct, sample_ct, gflags, &spgw, &alloc));
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    VariantData& vd = vds[vidx];
    if (is_multi[vidx]) {
      CHECK(!SpgwAppendMultiallelicGenovecDphase16(&vd.pgv, vd.allele_ct, &spgw));
    } else {
      // biallelic: existing API, plain genotypes only (no patches)
      PgenVariant& pgv = vd.pgv;
      ZeroWArr(BitCtToWordCt(sample_ct), pgv.patch_01_set);
      CHECK(!SpgwAppendBiallelicGenovecDosage16(pgv.genovec, pgv.dosage_present, pgv.dosage_main, pgv.dosage_ct, &spgw));
    }
  }
  PglErr reterr = SpgwFinish(&spgw);
  CHECK(!reterr);
  CleanupSpgw(&spgw, &reterr);
  free(alloc);

  // Hand-parse the header: mode 0x10, single vblock, 8-bit vrtypes, then
  // record lengths.  (Allele counts aren't stored in the header; as usual,
  // they'd come from the .pvar.)
  std::vector<unsigned char> file = ReadFile(fname);
  CHECK(file.size() > 12);
  CHECK(file[2] == 0x10);
  uint32_t file_variant_ct, file_sample_ct;
  memcpy(&file_variant_ct, &file[3], 4);
  memcpy(&file_sample_ct, &file[7], 4);
  CHECK(file_variant_ct == variant_ct);
  CHECK(file_sample_ct == sample_ct);
  const uint32_t ctrl = file[11];
  CHECK((ctrl & 12) == 4);  // 8-bit vrtypes
  CHECK(((ctrl >> 4) & 3) == 0);  // no allele counts
  const uint32_t len_byte_ct = 1 + (ctrl & 3);
  uint64_t fpos;
  memcpy(&fpos, &file[12], 8);
  const unsigned char* vrtypes = &file[20];
  const unsigned char* lens = &vrtypes[variant_ct];
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    uint32_t len = 0;
    memcpy(&len, &lens[vidx * len_byte_ct], len_byte_ct);
    CHECK(fpos + len <= file.size());
    if (is_multi[vidx]) {
      const VariantData& vd = vds[vidx];
      const uint32_t vrtype = vrtypes[vidx];
      CHECK((vrtype & 0xe8) == vd.expected_vrtype_bits);
      CHECK(((vrtype >> 4) & 1) == (vd.pgv.phasepresent_ct != 0));
      const std::vector<unsigned char>& suffix = vd.expected_suffix;
      CHECK(len >= suffix.size());
      if (len >= suffix.size()) {
        const int ok = !memcmp(&file[fpos + len - suffix.size()], suffix.data(), suffix.size());
        if (!ok) {
          fprintf(stderr, "variant %u (allele_ct %u, dosage_ct %u, dphase_ct %u, md %u, mdp %u): record suffix mismatch\n", vidx, vd.allele_ct, vd.pgv.dosage_ct, vd.pgv.dphase_ct, vd.pgv.multidosage_sample_ct, vd.pgv.multidphase_sample_ct);
        }
        CHECK(ok);
      }
    }
    fpos += len;
  }
  CHECK(fpos == file.size());
  uint64_t md_entry_ct = 0;
  uint64_t mdp_entry_ct = 0;
  uint32_t uncond_ct = 0;
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    if (is_multi[vidx]) {
      md_entry_ct += vds[vidx].md_codes.size();
      mdp_entry_ct += vds[vidx].mdp_codes.size();
      uncond_ct += ((vds[vidx].expected_vrtype_bits & 0x60) == 0x40);
    }
  }
  printf("%s: %u variants, %" PRIu64 " rare-alt dosage entries, %" PRIu64 " rare-alt dphase entries, %u unconditional-dosage records\n", fname, variant_ct, md_entry_ct, mdp_entry_ct, uncond_ct);
  for (VariantData& vd : vds) {
    FreeVariant(vd);
  }
}

int main() {
  // 300 samples: 2-byte group starts for most allele counts.  70000 samples:
  // 3-byte group starts, and sample IDs past 2^16.
  RunFile("tmp_small.pgen", 300, 6 * 36, 0);
  // combos 24-35: dense rare-alt dosages
  RunFile("tmp_big.pgen", 70000, 6 * 12, 24);

  // Error paths.
  {
    VariantData vd;
    MakeVariant(300, 5, 1, 0, 0, 2, vd);
    CHECK(vd.pgv.multidosage_sample_ct);
    uintptr_t allele_idx_offsets[2] = {0, 5};
    STPgenWriter spgw;
    unsigned char* alloc;
    // Nonempty rare-alt dosages without the opt-in flag.
    CHECK(!OpenWriter("tmp_noflag.pgen", allele_idx_offsets, 1, 300, kfPgenGlobalDosagePresent, &spgw, &alloc));
    CHECK(SpgwAppendMultiallelicGenovecDphase16(&vd.pgv, 5, &spgw) == kPglRetVarRecordTooLarge);
    PglErr reterr = kPglRetSuccess;
    CleanupSpgw(&spgw, &reterr);
    free(alloc);
    // Decreasing codes within a sample.
    CHECK(!OpenWriter("tmp_badcode.pgen", allele_idx_offsets, 1, 300, kfPgenGlobalDosagePresent | kfPgenGlobalMultiallelicDosagePresent, &spgw, &alloc));
    size_t i = 0;
    for (; i != vd.md_cts.size(); ++i) {
      if (vd.md_cts[i] >= 2) {
        break;
      }
    }
    CHECK(i != vd.md_cts.size());
    size_t offset = 0;
    for (size_t j = 0; j != i; ++j) {
      offset += vd.md_cts[j];
    }
    std::swap(vd.md_codes[offset], vd.md_codes[offset + 1]);
    CHECK(SpgwAppendMultiallelicGenovecDphase16(&vd.pgv, 5, &spgw) == kPglRetVarRecordTooLarge);
    std::swap(vd.md_codes[offset], vd.md_codes[offset + 1]);
    // Code out of range.
    const AlleleCode saved = vd.md_codes[0];
    vd.md_codes[0] = 5;
    CHECK(SpgwAppendMultiallelicGenovecDphase16(&vd.pgv, 5, &spgw) == kPglRetVarRecordTooLarge);
    vd.md_codes[0] = saved;
    CleanupSpgw(&spgw, &reterr);
    free(alloc);
    // Empty rare-alt lists are fine without the flag.
    VariantData vd2;
    MakeVariant(300, 5, 1, 1, 1, 0, vd2);
    CHECK(!vd2.pgv.multidosage_sample_ct);
    CHECK(!OpenWriter("tmp_empty.pgen", allele_idx_offsets, 1, 300, kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent, &spgw, &alloc));
    CHECK(!SpgwAppendMultiallelicGenovecDphase16(&vd2.pgv, 5, &spgw));
    reterr = SpgwFinish(&spgw);
    CHECK(!reterr);
    CleanupSpgw(&spgw, &reterr);
    free(alloc);
    FreeVariant(vd);
    FreeVariant(vd2);
  }

  if (g_fail_ct) {
    fprintf(stderr, "%d check(s) failed.\n", g_fail_ct);
    return 1;
  }
  printf("multidosage writer tests passed.\n");
  return 0;
}
