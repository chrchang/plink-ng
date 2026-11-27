// Round-trip test for multiallelic dosages (rare-alt dosage tracks #5-6,
// #9-10).  Writes .pgen files with SpgwAppendMultiallelicGenovecDphase16()
// from random but valid data, and checks each record's type byte and that the
// record ends with independently computed bytes for tracks #5-10 (parsing the
// header by hand).  Then reads the files back with PgrGetMDp(), PgrGetMD(),
// PgrGetDp(), PgrGet1D(), PgrGetInv1D(), and PgrGetInv1Dp(), with and without
// sample subsetting, and compares against the source data.

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "pgenlib_read.h"
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
  // offset of track #6 within expected_suffix
  size_t md_vals_suffix_offset;
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
  std::vector<uint32_t> is_het(sample_ct, 0);
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
      is_het[s] = (lo != hi);
    }
    het_ct += (cats[s] == 1);
    is_het[s] |= (cats[s] == 1);
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
  // No implicit dosage-phase for multiallelic variants: a phased het with a
  // dosage needs an explicit dosage-phase value.  (All hets are phased when
  // hphase is set.)
  const uint32_t het_dosage_needs_dphase = hphase && het_ct;
  for (uint32_t s = 0; s != sample_ct; ++s) {
    if (het_dosage_needs_dphase && (!dphase) && (is_het[s] || dosage_missing[s])) {
      // Dropping the hets' dosages means the record can't use unconditional
      // dosage, so missing (65535) entries must go too.
      has_dosage[s] = 0;
    }
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
    if (dphase && ((dosage_mode == 2) || (het_dosage_needs_dphase && is_het[s]) || RandBelow(2))) {
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
    vd.md_vals_suffix_offset = vd.expected_suffix.size() - 2 * vd.md_vals.size();
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

// ***** Read-back *****

// One sample's data, independent of PgenVariant encoding details.
struct SampleData {
  uint32_t a1 = 0, a2 = 0;  // hardcall alleles, a1 <= a2; missing = UINT32_MAX
  uint32_t phased = 0, phaseinfo = 0;
  uint32_t has_dosage = 0, dosage = 0;
  uint32_t has_dphase = 0;
  int32_t dphase = 0;
  std::vector<std::pair<uint32_t, uint32_t>> md;  // (allele, value)
  std::vector<std::pair<uint32_t, int32_t>> mdp;
  bool operator==(const SampleData& o) const {
    return (a1 == o.a1) && (a2 == o.a2) && (phased == o.phased) && (phaseinfo == o.phaseinfo) && (has_dosage == o.has_dosage) && (dosage == o.dosage) && (has_dphase == o.has_dphase) && (dphase == o.dphase) && (md == o.md) && (mdp == o.mdp);
  }
};

// Decodes a PgenVariant.  As in the writer's input, phasepresent == nullptr
// with phasepresent_ct > 0 means every het is phased, and 65535 dosages /
// -32768 deltas are missing.  Fields not selected by load_* are left zero.
static std::vector<SampleData> Decode(const PgenVariant& pgv, uint32_t n, uint32_t load_phase, uint32_t load_dphase, uint32_t load_multi) {
  std::vector<SampleData> out(n);
  const AlleleCode* p01 = pgv.patch_01_vals;
  const AlleleCode* p10 = pgv.patch_10_vals;
  uint32_t dosage_idx = 0;
  uint32_t dphase_idx = 0;
  uint32_t md_sample_idx = 0;
  uint32_t md_entry_idx = 0;
  uint32_t mdp_sample_idx = 0;
  uint32_t mdp_entry_idx = 0;
  for (uint32_t s = 0; s != n; ++s) {
    SampleData& sd = out[s];
    const uint32_t cat = GetNyparrEntry(pgv.genovec, s);
    if (cat == 0) {
      sd.a1 = 0;
      sd.a2 = 0;
    } else if (cat == 1) {
      sd.a1 = 0;
      sd.a2 = (pgv.patch_01_ct && IsSet(pgv.patch_01_set, s))? (*p01++) : 1;
    } else if (cat == 2) {
      if (pgv.patch_10_ct && IsSet(pgv.patch_10_set, s)) {
        sd.a1 = p10[0];
        sd.a2 = p10[1];
        p10 += 2;
      } else {
        sd.a1 = 1;
        sd.a2 = 1;
      }
    } else {
      sd.a1 = UINT32_MAX;
      sd.a2 = UINT32_MAX;
    }
    if (load_phase && pgv.phasepresent_ct) {
      const uint32_t het = (cat != 3) && (sd.a1 != sd.a2);
      sd.phased = pgv.phasepresent? IsSet(pgv.phasepresent, s) : het;
      if (sd.phased) {
        sd.phaseinfo = IsSet(pgv.phaseinfo, s);
      }
    }
    if (pgv.dosage_ct && IsSet(pgv.dosage_present, s)) {
      const uint32_t v = pgv.dosage_main[dosage_idx++];
      if (v != 65535) {
        sd.has_dosage = 1;
        sd.dosage = v;
      }
    }
    if (load_dphase && pgv.dphase_ct && IsSet(pgv.dphase_present, s)) {
      const int32_t d = pgv.dphase_delta[dphase_idx++];
      if (d != -32768) {
        sd.has_dphase = 1;
        sd.dphase = d;
      }
    }
    if (load_multi && pgv.multidosage_sample_ct && IsSet(pgv.multidosage_present, s)) {
      const uint32_t ct = pgv.multidosage_cts[md_sample_idx++];
      for (uint32_t j = 0; j != ct; ++j, ++md_entry_idx) {
        sd.md.emplace_back(pgv.multidosage_codes[md_entry_idx], pgv.multidosage_vals[md_entry_idx]);
      }
    }
    if (load_multi && load_dphase && pgv.multidphase_sample_ct && IsSet(pgv.multidphase_present, s)) {
      const uint32_t ct = pgv.multidphase_cts[mdp_sample_idx++];
      for (uint32_t j = 0; j != ct; ++j, ++mdp_entry_idx) {
        sd.mdp.emplace_back(pgv.multidphase_codes[mdp_entry_idx], pgv.multidphase_delta[mdp_entry_idx]);
      }
    }
  }
  return out;
}

static int g_mismatch_msg_ct = 0;

static void CheckSamples(const std::vector<SampleData>& expected, const std::vector<SampleData>& actual, const char* what, uint32_t vidx) {
  CHECK(expected.size() == actual.size());
  for (size_t i = 0; i != expected.size(); ++i) {
    if (!(expected[i] == actual[i])) {
      if (g_mismatch_msg_ct++ < 10) {
        const SampleData& e = expected[i];
        const SampleData& a = actual[i];
        fprintf(stderr, "%s mismatch: variant %u, sample idx %zu: expected (%u/%u ph %u/%u d %u:%u dp %u:%d md %zu mdp %zu), got (%u/%u ph %u/%u d %u:%u dp %u:%d md %zu mdp %zu)\n", what, vidx, i, e.a1, e.a2, e.phased, e.phaseinfo, e.has_dosage, e.dosage, e.has_dphase, e.dphase, e.md.size(), e.mdp.size(), a.a1, a.a2, a.phased, a.phaseinfo, a.has_dosage, a.dosage, a.has_dphase, a.dphase, a.md.size(), a.mdp.size());
      }
      ++g_fail_ct;
      return;
    }
  }
}

struct ReaderBufs {
  std::vector<uintptr_t*> bufs;
  PgenVariant pgv;
};

static uintptr_t* ReaderAlloc(ReaderBufs& rb, uint64_t bit_ct) {
  uintptr_t* ptr = AllocWords(bit_ct);
  rb.bufs.push_back(ptr);
  return ptr;
}

static void InitReaderBufs(uint32_t sample_ct, uint32_t max_allele_ct, ReaderBufs& rb) {
  PreinitPgv(&rb.pgv);
  PgenVariant& pgv = rb.pgv;
  const uint64_t max_entry_ct = S_CAST(uint64_t, sample_ct) * MINV(max_allele_ct - 2, 255);
  pgv.genovec = ReaderAlloc(rb, 2 * sample_ct);
  pgv.patch_01_set = ReaderAlloc(rb, sample_ct);
  pgv.patch_01_vals = R_CAST(AlleleCode*, ReaderAlloc(rb, 8 * sample_ct));
  pgv.patch_10_set = ReaderAlloc(rb, sample_ct);
  pgv.patch_10_vals = R_CAST(AlleleCode*, ReaderAlloc(rb, 16 * sample_ct));
  pgv.phasepresent = ReaderAlloc(rb, sample_ct);
  pgv.phaseinfo = ReaderAlloc(rb, sample_ct);
  pgv.dosage_present = ReaderAlloc(rb, sample_ct);
  pgv.dosage_main = R_CAST(uint16_t*, ReaderAlloc(rb, 16 * sample_ct));
  pgv.dphase_present = ReaderAlloc(rb, sample_ct);
  pgv.dphase_delta = R_CAST(int16_t*, ReaderAlloc(rb, 16 * sample_ct));
  pgv.multidosage_present = ReaderAlloc(rb, sample_ct);
  pgv.multidosage_cts = R_CAST(unsigned char*, ReaderAlloc(rb, 8 * sample_ct));
  pgv.multidosage_codes = R_CAST(AlleleCode*, ReaderAlloc(rb, 8 * max_entry_ct));
  pgv.multidosage_vals = R_CAST(uint16_t*, ReaderAlloc(rb, 16 * max_entry_ct));
  pgv.multidphase_present = ReaderAlloc(rb, sample_ct);
  pgv.multidphase_cts = R_CAST(unsigned char*, ReaderAlloc(rb, 8 * sample_ct));
  pgv.multidphase_codes = R_CAST(AlleleCode*, ReaderAlloc(rb, 8 * max_entry_ct));
  pgv.multidphase_delta = R_CAST(int16_t*, ReaderAlloc(rb, 16 * max_entry_ct));
}

static void FreeReaderBufs(ReaderBufs& rb) {
  for (uintptr_t* p : rb.bufs) {
    free(p);
  }
  rb.bufs.clear();
}

struct Reader {
  PgenFileInfo pgfi;
  PgenReader pgr;
  unsigned char* pgfi_alloc;
  unsigned char* pgr_alloc;
};

static PglErr OpenReader(const char* fname, const uintptr_t* allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, Reader& r) {
  char errstr[kPglErrstrBufBlen];
  PreinitPgfi(&r.pgfi);
  PreinitPgr(&r.pgr);
  r.pgfi_alloc = nullptr;
  r.pgr_alloc = nullptr;
  PgenHeaderCtrl header_ctrl;
  uintptr_t pgfi_alloc_cacheline_ct;
  PglErr reterr = PgfiInitPhase1(fname, nullptr, variant_ct, sample_ct, &header_ctrl, &r.pgfi, &pgfi_alloc_cacheline_ct, errstr);
  if (reterr) {
    fputs(errstr, stderr);
    return reterr;
  }
  // the reader doesn't modify these when allele_cts_already_loaded is set
  r.pgfi.allele_idx_offsets = const_cast<uintptr_t*>(allele_idx_offsets);
  uint32_t max_allele_ct = 2;
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    max_allele_ct = MAXV(max_allele_ct, allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]);
  }
  r.pgfi.max_allele_ct = max_allele_ct;
  r.pgfi_alloc = static_cast<unsigned char*>(aligned_alloc(kCacheline, (pgfi_alloc_cacheline_ct + 1) * kCacheline));
  uint32_t max_vrec_width;
  uintptr_t pgr_alloc_cacheline_ct;
  reterr = PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, variant_ct, &max_vrec_width, &r.pgfi, r.pgfi_alloc, &pgr_alloc_cacheline_ct, errstr);
  if (reterr) {
    fputs(errstr, stderr);
    return reterr;
  }
  r.pgr_alloc = static_cast<unsigned char*>(aligned_alloc(kCacheline, (pgr_alloc_cacheline_ct + 1) * kCacheline));
  return PgrInit(fname, max_vrec_width, &r.pgfi, &r.pgr, r.pgr_alloc);
}

static void CloseReader(Reader& r) {
  PglErr reterr = kPglRetSuccess;
  CleanupPgr(&r.pgr, &reterr);
  CleanupPgfi(&r.pgfi, &reterr);
  free(r.pgr_alloc);
  free(r.pgfi_alloc);
}

// Expected single-allele (inverse-)dosage data.
static void ExpectedAllele(const std::vector<SampleData>& full, uint32_t allele_idx, uint32_t invert, std::vector<SampleData>& out) {
  out.assign(full.size(), SampleData());
  for (size_t i = 0; i != full.size(); ++i) {
    const SampleData& f = full[i];
    SampleData& o = out[i];
    // a1 = allele count (3 = missing)
    if (f.a1 == UINT32_MAX) {
      o.a1 = 3;
    } else {
      o.a1 = (f.a1 == allele_idx) + (f.a2 == allele_idx);
      if (invert) {
        o.a1 = 2 - o.a1;
      }
    }
    if (f.has_dosage) {
      int32_t v;
      if (!allele_idx) {
        v = 32768 - f.dosage;
      } else if (allele_idx == 1) {
        v = f.dosage;
        for (const auto& pr : f.md) {
          v -= pr.second;
        }
      } else {
        v = 0;
        for (const auto& pr : f.md) {
          if (pr.first == allele_idx) {
            v = pr.second;
          }
        }
      }
      o.has_dosage = 1;
      o.dosage = invert? (32768 - v) : v;
    }
    if (f.has_dphase && allele_idx) {
      int32_t d;
      if (allele_idx == 1) {
        d = f.dphase;
        for (const auto& pr : f.mdp) {
          d -= pr.second;
        }
      } else {
        d = 0;
        for (const auto& pr : f.mdp) {
          if (pr.first == allele_idx) {
            d = pr.second;
          }
        }
      }
      o.has_dphase = 1;
      o.dphase = invert? -d : d;
    }
  }
}

// Decodes single-allele output (count in a1).
static void DecodeAllele(const uintptr_t* countvec, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, const uintptr_t* dphase_present, const int16_t* dphase_delta, uint32_t dphase_ct, uint32_t n, std::vector<SampleData>& out) {
  out.assign(n, SampleData());
  uint32_t dosage_idx = 0;
  uint32_t dphase_idx = 0;
  for (uint32_t s = 0; s != n; ++s) {
    out[s].a1 = GetNyparrEntry(countvec, s);
    if (dosage_ct && IsSet(dosage_present, s)) {
      out[s].has_dosage = 1;
      out[s].dosage = dosage_main[dosage_idx++];
    }
    if (dphase_ct && IsSet(dphase_present, s)) {
      out[s].has_dphase = 1;
      out[s].dphase = dphase_delta[dphase_idx++];
    }
  }
  CHECK(dosage_idx == dosage_ct);
  CHECK(dphase_idx == dphase_ct);
}

static std::vector<SampleData> SubsetSamples(const std::vector<SampleData>& full, const uintptr_t* sample_include) {
  std::vector<SampleData> out;
  for (size_t i = 0; i != full.size(); ++i) {
    if (!sample_include || IsSet(sample_include, i)) {
      out.push_back(full[i]);
    }
  }
  return out;
}

// Reads every variant back, with sample_include == nullptr (all samples) or a
// subset, and compares against the source data.
static void ReadBackPass(Reader& r, const std::vector<VariantData>& vds, const std::vector<uint32_t>& is_multi, uint32_t sample_ct, const uintptr_t* sample_include, uint32_t subset_ct, const uint32_t* cumulative_popcounts, ReaderBufs& rb) {
  PgrSampleSubsetIndex pssi;
  if (sample_include) {
    PgrSetSampleSubsetIndex(cumulative_popcounts, &r.pgr, &pssi);
  } else {
    PgrClearSampleSubsetIndex(&r.pgr, &pssi);
  }
  PgenVariant& out = rb.pgv;
  const uint32_t variant_ct = vds.size();
  std::vector<SampleData> expected;
  std::vector<SampleData> actual;
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    const VariantData& vd = vds[vidx];
    const uint32_t allele_ct = vd.allele_ct;
    if (!is_multi[vidx]) {
      // biallelic records were written with genotypes and unphased dosages
      // only
      PgenVariant src = vd.pgv;
      src.patch_01_ct = 0;
      src.patch_10_ct = 0;
      src.phasepresent_ct = 0;
      src.dphase_ct = 0;
      src.multidosage_sample_ct = 0;
      src.multidphase_sample_ct = 0;
      CHECK(!PgrGetMDp(sample_include, pssi, subset_ct, vidx, &r.pgr, &out));
      CheckSamples(SubsetSamples(Decode(src, sample_ct, 1, 1, 1), sample_include), Decode(out, subset_ct, 1, 1, 1), "biallelic PgrGetMDp", vidx);
      continue;
    }
    const std::vector<SampleData> full = SubsetSamples(Decode(vd.pgv, sample_ct, 1, 1, 1), sample_include);
    // PgrGetMDp(): everything
    PglErr reterr = PgrGetMDp(sample_include, pssi, subset_ct, vidx, &r.pgr, &out);
    CHECK(!reterr);
    if (!reterr) {
      CheckSamples(full, Decode(out, subset_ct, 1, 1, 1), "PgrGetMDp", vidx);
    }
    // PgrGetMD(): no phase
    reterr = PgrGetMD(sample_include, pssi, subset_ct, vidx, &r.pgr, &out);
    CHECK(!reterr);
    if (!reterr) {
      expected = full;
      for (SampleData& sd : expected) {
        sd.phased = 0;
        sd.phaseinfo = 0;
        sd.has_dphase = 0;
        sd.dphase = 0;
        sd.mdp.clear();
      }
      CheckSamples(expected, Decode(out, subset_ct, 0, 0, 1), "PgrGetMD", vidx);
    }
    // PgrGetDp(): collapsed to ref/nonref
    reterr = PgrGetDp(sample_include, pssi, subset_ct, vidx, &r.pgr, &out);
    CHECK(!reterr);
    if (!reterr) {
      expected = full;
      for (SampleData& sd : expected) {
        if (sd.a1 != UINT32_MAX) {
          // alt/alt hets aren't hets after collapsing
          if (sd.a1) {
            sd.phased = 0;
            sd.phaseinfo = 0;
          }
          sd.a1 = (sd.a1 != 0);
          sd.a2 = (sd.a2 != 0);
        }
        sd.md.clear();
        sd.mdp.clear();
      }
      out.patch_01_ct = 0;
      out.patch_10_ct = 0;
      out.multidosage_sample_ct = 0;
      out.multidphase_sample_ct = 0;
      CheckSamples(expected, Decode(out, subset_ct, 1, 1, 0), "PgrGetDp", vidx);
    }
    // Single-allele functions, for a few alleles.
    const uint32_t allele_idxs[] = {0, 1, 2, allele_ct / 2, allele_ct - 1};
    for (uint32_t allele_idx : allele_idxs) {
      reterr = PgrGet1D(sample_include, pssi, subset_ct, vidx, allele_idx, &r.pgr, out.genovec, out.dosage_present, out.dosage_main, &out.dosage_ct);
      CHECK(!reterr);
      if (!reterr) {
        ExpectedAllele(full, allele_idx, 0, expected);
        for (SampleData& sd : expected) {
          sd.has_dphase = 0;
          sd.dphase = 0;
        }
        DecodeAllele(out.genovec, out.dosage_present, out.dosage_main, out.dosage_ct, nullptr, nullptr, 0, subset_ct, actual);
        CheckSamples(expected, actual, "PgrGet1D", vidx);
      }
      reterr = PgrGetInv1D(sample_include, pssi, subset_ct, vidx, allele_idx, &r.pgr, out.genovec, out.dosage_present, out.dosage_main, &out.dosage_ct);
      CHECK(!reterr);
      if (!reterr) {
        ExpectedAllele(full, allele_idx, 1, expected);
        for (SampleData& sd : expected) {
          sd.has_dphase = 0;
          sd.dphase = 0;
        }
        DecodeAllele(out.genovec, out.dosage_present, out.dosage_main, out.dosage_ct, nullptr, nullptr, 0, subset_ct, actual);
        CheckSamples(expected, actual, "PgrGetInv1D", vidx);
      }
      if (allele_idx) {
        reterr = PgrGetInv1Dp(sample_include, pssi, subset_ct, vidx, allele_idx, &r.pgr, &out);
        CHECK(!reterr);
        if (!reterr) {
          ExpectedAllele(full, allele_idx, 1, expected);
          DecodeAllele(out.genovec, out.dosage_present, out.dosage_main, out.dosage_ct, out.dphase_present, out.dphase_delta, out.dphase_ct, subset_ct, actual);
          CheckSamples(expected, actual, "PgrGetInv1Dp", vidx);
        }
      }
    }
  }
}

static void ReadBack(const char* fname, const std::vector<VariantData>& vds, const std::vector<uint32_t>& is_multi, const uintptr_t* allele_idx_offsets, uint32_t sample_ct) {
  const uint32_t variant_ct = vds.size();
  Reader r;
  PglErr reterr = OpenReader(fname, allele_idx_offsets, variant_ct, sample_ct, r);
  CHECK(!reterr);
  if (reterr) {
    return;
  }
  ReaderBufs rb;
  InitReaderBufs(sample_ct, 255, rb);
  ReadBackPass(r, vds, is_multi, sample_ct, nullptr, sample_ct, nullptr, rb);
  // Random subset of about 3/4 of the samples (so that subsetted dosage
  // lists stay nonempty), with the first and last samples excluded.
  uintptr_t* sample_include = AllocWords(sample_ct);
  uint32_t subset_ct = 0;
  for (uint32_t s = 1; s + 1 < sample_ct; ++s) {
    if (RandBelow(4)) {
      SetBit(s, sample_include);
      ++subset_ct;
    }
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  std::vector<uint32_t> cumulative_popcounts(sample_ctl);
  FillCumulativePopcounts(sample_include, sample_ctl, cumulative_popcounts.data());
  ReadBackPass(r, vds, is_multi, sample_ct, sample_include, subset_ct, cumulative_popcounts.data(), rb);

  // Rare-alt dosages without the opt-in buffers.
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    if (is_multi[vidx] && vds[vidx].pgv.dosage_ct) {
      PgrSampleSubsetIndex pssi;
      PgrClearSampleSubsetIndex(&r.pgr, &pssi);
      PgenVariant pgv = rb.pgv;
      pgv.multidosage_present = nullptr;
      CHECK(PgrGetMD(nullptr, pssi, sample_ct, vidx, &r.pgr, &pgv) == kPglRetNotYetSupported);
      CHECK(PgrGetMDp(nullptr, pssi, sample_ct, vidx, &r.pgr, &pgv) == kPglRetNotYetSupported);
      break;
    }
  }
  free(sample_include);
  FreeReaderBufs(rb);
  CloseReader(r);
}

// Corrupts one track #6 value of a copy of fname, and checks that the readers
// reject it.
static void CheckCorruptValue(const char* fname, const std::vector<VariantData>& vds, const std::vector<uint32_t>& is_multi, const std::vector<uint64_t>& record_ends, const uintptr_t* allele_idx_offsets, uint32_t sample_ct, uint32_t bad_val) {
  const uint32_t variant_ct = vds.size();
  uint32_t vidx = 0;
  for (; vidx != variant_ct; ++vidx) {
    if (is_multi[vidx] && (!vds[vidx].md_vals.empty())) {
      break;
    }
  }
  CHECK(vidx != variant_ct);
  if (vidx == variant_ct) {
    return;
  }
  const VariantData& vd = vds[vidx];
  std::vector<unsigned char> file = ReadFile(fname);
  const uint64_t pos = record_ends[vidx] - vd.expected_suffix.size() + vd.md_vals_suffix_offset;
  file[pos] = bad_val & 255;
  file[pos + 1] = bad_val >> 8;
  const char* corrupt_fname = "tmp_corrupt.pgen";
  FILE* ff = fopen(corrupt_fname, "wb");
  CHECK(ff && (fwrite(file.data(), 1, file.size(), ff) == file.size()));
  fclose(ff);
  Reader r;
  CHECK(!OpenReader(corrupt_fname, allele_idx_offsets, variant_ct, sample_ct, r));
  ReaderBufs rb;
  InitReaderBufs(sample_ct, 255, rb);
  PgrSampleSubsetIndex pssi;
  PgrClearSampleSubsetIndex(&r.pgr, &pssi);
  CHECK(PgrGetMDp(nullptr, pssi, sample_ct, vidx, &r.pgr, &rb.pgv) == kPglRetMalformedInput);
  CHECK(PgrGet1D(nullptr, pssi, sample_ct, vidx, 1, &r.pgr, rb.pgv.genovec, rb.pgv.dosage_present, rb.pgv.dosage_main, &rb.pgv.dosage_ct) == kPglRetMalformedInput);
  // the corrupted entry's own allele
  CHECK(PgrGet1D(nullptr, pssi, sample_ct, vidx, vd.md_codes[0], &r.pgr, rb.pgv.genovec, rb.pgv.dosage_present, rb.pgv.dosage_main, &rb.pgv.dosage_ct) == kPglRetMalformedInput);
  FreeReaderBufs(rb);
  CloseReader(r);
}

// Writes a file mixing biallelic records (existing API) and multiallelic
// dosage records, then checks it.
// combo_offset selects which (dosage form, dphase, hphase, rare density)
// combinations the multiallelic records cycle through.
static void RunFile(const char* fname, uint32_t sample_ct, uint32_t variant_ct, uint32_t combo_offset, uint32_t check_corrupt) {
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
  std::vector<uint64_t> record_ends(variant_ct);
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
    record_ends[vidx] = fpos;
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
  ReadBack(fname, vds, is_multi, allele_idx_offsets.data(), sample_ct);
  if (check_corrupt) {
    CheckCorruptValue(fname, vds, is_multi, record_ends, allele_idx_offsets.data(), sample_ct, 0);
    CheckCorruptValue(fname, vds, is_multi, record_ends, allele_idx_offsets.data(), sample_ct, 32769);
  }
  for (VariantData& vd : vds) {
    FreeVariant(vd);
  }
}

int main() {
  // 300 samples: 2-byte group starts for most allele counts.  70000 samples:
  // 3-byte group starts, and sample IDs past 2^16.
  RunFile("tmp_small.pgen", 300, 6 * 36, 0, 1);
  // combos 24-35: dense rare-alt dosages
  RunFile("tmp_big.pgen", 70000, 6 * 12, 24, 0);

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
    // Phased het with a dosage but no dosage-phase value.
    {
      VariantData vd3;
      MakeVariant(300, 5, 1, 1, 1, 1, vd3);
      CHECK(vd3.pgv.phasepresent_ct && vd3.pgv.dphase_ct);
      // up to 3 records are appended below (the failed attempts don't count)
      const uintptr_t badphase_allele_idx_offsets[4] = {0, 5, 10, 15};
      CHECK(!OpenWriter("tmp_badphase.pgen", badphase_allele_idx_offsets, 3, 300, kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent | kfPgenGlobalMultiallelicDosagePresent, &spgw, &alloc));
      PgenVariant pgv = vd3.pgv;
      pgv.dphase_ct = 0;
      pgv.multidphase_sample_ct = 0;
      CHECK(SpgwAppendMultiallelicGenovecDphase16(&pgv, 5, &spgw) == kPglRetVarRecordTooLarge);
      // ...but it's fine once hardcall phase is dropped.
      pgv.phasepresent_ct = 0;
      CHECK(!SpgwAppendMultiallelicGenovecDphase16(&pgv, 5, &spgw));
      // Explicit phasepresent: phasing only the hets without dosages is fine;
      // adding one het with a dosage isn't.
      const std::vector<SampleData> sds = Decode(vd3.pgv, 300, 0, 0, 0);
      uintptr_t* phasepresent = AllocWords(300);
      uint32_t het_without_dosage_ct = 0;
      uint32_t het_with_dosage = UINT32_MAX;
      for (uint32_t s = 0; s != 300; ++s) {
        if ((sds[s].a1 == UINT32_MAX) || (sds[s].a1 == sds[s].a2)) {
          continue;
        }
        if (!IsSet(pgv.dosage_present, s)) {
          SetBit(s, phasepresent);
          ++het_without_dosage_ct;
        } else {
          het_with_dosage = s;
        }
      }
      CHECK(het_without_dosage_ct && (het_with_dosage != UINT32_MAX));
      pgv.phasepresent = phasepresent;
      pgv.phasepresent_ct = het_without_dosage_ct;
      CHECK(!SpgwAppendMultiallelicGenovecDphase16(&pgv, 5, &spgw));
      SetBit(het_with_dosage, phasepresent);
      pgv.phasepresent_ct = het_without_dosage_ct + 1;
      CHECK(SpgwAppendMultiallelicGenovecDphase16(&pgv, 5, &spgw) == kPglRetVarRecordTooLarge);
      free(phasepresent);
      CleanupSpgw(&spgw, &reterr);
      free(alloc);
      FreeVariant(vd3);
    }
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
