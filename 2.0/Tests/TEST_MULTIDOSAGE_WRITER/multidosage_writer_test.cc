// Round-trip test for multiallelic dosages (rare-alt dosage tracks #5-6,
// #9-10).  Writes .pgen files with SpgwAppendMultiallelicGenovecDphase16()
// from random but valid data, and checks each record's type byte and that the
// record ends with independently computed bytes for tracks #5-10 (parsing the
// header by hand).  Then reads the files back with PgrGetMDp(), PgrGetMD(),
// PgrGetDp(), PgrGet1D(), PgrGetInv1D(), PgrGetInv1Dp(), PgrGetMDCounts(), and
// PgrGetDCounts(),
// with and without sample subsetting, and with PgrGetRaw(), and compares
// against the source data, and checks that
// PgrValidate() accepts the files.  Finally, PgrValidate() is run on small
// hand-built records that each violate one multiallelic-dosage rule.

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <set>
#include <string>
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
  // hardcall alleles (UINT32_MAX = missing)
  std::vector<uint32_t> hc_a1(sample_ct);
  std::vector<uint32_t> hc_a2(sample_ct);
  for (uint32_t s = 0; s != sample_ct; ++s) {
    pgv.genovec[s / kBitsPerWordD2] |= S_CAST(uintptr_t, cats[s]) << (2 * (s % kBitsPerWordD2));
    hc_a1[s] = (cats[s] == 3)? UINT32_MAX : (cats[s] == 2);
    hc_a2[s] = (cats[s] == 3)? UINT32_MAX : (cats[s] != 0);
    // No rare-alt hardcalls when there are to be no rare-alt dosages, since
    // they'd be inconsistent with each other.
    if (rare_density && (cats[s] == 1) && RandBelow(3) == 0) {
      SetBit(s, pgv.patch_01_set);
      vd.patch_01_vals.push_back(2 + RandBelow(rare_ct));
      hc_a2[s] = vd.patch_01_vals.back();
    } else if (rare_density && (cats[s] == 2) && RandBelow(3) == 0) {
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
      hc_a1[s] = lo;
      hc_a2[s] = hi;
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
  std::vector<uint32_t> u(allele_ct, 0);
  std::vector<uint32_t> touched;
  auto add_dosage = [&](uint32_t allele, uint32_t amt) {
    u[allele] += amt;
    touched.push_back(allele);
  };
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
    // Per-allele dosages u[], summing to 32768.  With a non-missing hardcall,
    // start from the hardcall and move at most half a copy's worth (8192)
    // away, so that sum_i |16384 * h_i - u_i| <= 16384.
    touched.clear();
    uint32_t v;
    if (cats[s] == 3) {
      v = RandBelow(32769);
      add_dosage(0, 32768 - v);
      uint32_t remaining = v;
      if (rare_density) {
        const uint32_t max_entries = (rare_density == 1)? (RandBelow(4) == 0) : MINV(rare_ct, 255);
        uint32_t entry_ct = 0;
        for (uint32_t code = 2; (code < allele_ct) && (entry_ct < max_entries) && remaining; ++code) {
          if ((rare_density == 2) || RandBelow(rare_ct) == 0) {
            const uint32_t cur = 1 + RandBelow(MINV(remaining, 1 + (v / 4)));
            add_dosage(code, cur);
            remaining -= cur;
            ++entry_ct;
          }
        }
      }
      add_dosage(1, remaining);
    } else {
      const uint32_t a1 = hc_a1[s];
      const uint32_t a2 = hc_a2[s];
      add_dosage(a1, 16384);
      add_dosage(a2, 16384);
      uint32_t budget = RandBelow(8193);
      const uint32_t transfer_ct = (rare_density == 0)? 1 : ((rare_density == 1)? (1 + (RandBelow(4) == 0)) : MINV(rare_ct, 255));
      for (uint32_t t = 0; (t != transfer_ct) && budget; ++t) {
        uint32_t src = RandBelow(2)? a1 : a2;
        if (!u[src]) {
          src = (src == a1)? a2 : a1;
        }
        uint32_t dst;
        if (rare_density == 0) {
          dst = (src == 0)? 1 : 0;
        } else if (rare_density == 2) {
          dst = 2 + t;
        } else {
          dst = RandBelow(allele_ct);
        }
        if ((dst == src) || (!u[src])) {
          continue;
        }
        const uint32_t amt = 1 + RandBelow(MINV(u[src], budget));
        u[src] -= amt;
        add_dosage(dst, amt);
        budget -= amt;
      }
      v = 32768 - u[0];
    }
    std::sort(touched.begin(), touched.end());
    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
    std::vector<uint32_t> nz_alleles;
    for (uint32_t a : touched) {
      if (u[a]) {
        nz_alleles.push_back(a);
      }
    }
    vd.dosage_main.push_back(v);
    // rare-alt dosages
    std::vector<uint32_t> codes;
    std::vector<uint32_t> vals;
    for (uint32_t a : nz_alleles) {
      if (a >= 2) {
        codes.push_back(a);
        vals.push_back(u[a]);
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
      // Split u[] into two haplotypes h1[], h2[], each summing to 16384 with
      // nonnegative entries: h1_i in [max(0, u_i - 16384), min(u_i, 16384)].
      std::vector<uint32_t> h1(nz_alleles.size());
      std::vector<uint32_t> caps(nz_alleles.size());
      uint32_t rem = 16384;
      for (size_t i = 0; i != nz_alleles.size(); ++i) {
        const uint32_t ui = u[nz_alleles[i]];
        h1[i] = (ui > 16384)? (ui - 16384) : 0;
        caps[i] = MINV(ui, 16384) - h1[i];
        rem -= h1[i];
      }
      for (uint32_t pass = 0; pass != 2; ++pass) {
        for (size_t i = 0; (i != nz_alleles.size()) && rem; ++i) {
          const uint32_t limit = MINV(caps[i], rem);
          const uint32_t give = pass? limit : RandBelow(limit + 1);
          h1[i] += give;
          caps[i] -= give;
          rem -= give;
        }
      }
      // delta_i = h1_i - h2_i = 2 * h1_i - u_i; the alt-sum delta is minus
      // the ref delta.
      int32_t delta = 0;
      std::vector<uint32_t> pcodes;
      std::vector<int32_t> pdeltas;
      for (size_t i = 0; i != nz_alleles.size(); ++i) {
        const uint32_t a = nz_alleles[i];
        const int32_t cur_delta = 2 * S_CAST(int32_t, h1[i]) - S_CAST(int32_t, u[a]);
        if (!a) {
          delta = -cur_delta;
        } else if ((a >= 2) && cur_delta) {
          pcodes.push_back(a);
          pdeltas.push_back(cur_delta);
        }
      }
      SetBit(s, pgv.dphase_present);
      vd.dphase_delta.push_back(delta);
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
    for (uint32_t a : touched) {
      u[a] = 0;
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

// Checks PgrGetMDCounts() against values computed from the source data.
static void CheckMDCounts(Reader& r, const uintptr_t* sample_include, const uintptr_t* sample_include_interleaved_vec, PgrSampleSubsetIndex pssi, uint32_t subset_ct, uint32_t vidx, uint32_t allele_ct, const std::vector<SampleData>& full) {
  std::vector<uint64_t> expected_sums(allele_ct, 0);
  std::vector<uint64_t> ssqs[2] = {std::vector<uint64_t>(allele_ct, 0), std::vector<uint64_t>(allele_ct, 0)};
  STD_ARRAY_DECL(uint32_t, 4, expected_genocounts);
  STD_ARRAY_FILL0(expected_genocounts);
  uint32_t expected_het_ct = 0;
  uint32_t nm_ct = 0;
  uint32_t extra_phased_het_ct = 0;
  for (const SampleData& sd : full) {
    const uint32_t hc_missing = (sd.a1 == UINT32_MAX);
    if (hc_missing) {
      expected_genocounts[3] += 1;
    } else {
      expected_genocounts[(sd.a1 != 0) + (sd.a2 != 0)] += 1;
      expected_het_ct += (sd.a1 != sd.a2);
    }
    if (sd.has_dosage) {
      std::vector<int64_t> u(allele_ct, 0);
      std::vector<int64_t> d(allele_ct, 0);
      u[0] = 32768 - sd.dosage;
      u[1] = sd.dosage;
      for (const auto& pr : sd.md) {
        u[pr.first] = pr.second;
        u[1] -= pr.second;
      }
      if (sd.has_dphase) {
        d[0] = -sd.dphase;
        d[1] = sd.dphase;
        for (const auto& pr : sd.mdp) {
          d[pr.first] = pr.second;
          d[1] -= pr.second;
        }
      }
      for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
        expected_sums[aidx] += u[aidx];
        ssqs[0][aidx] += u[aidx] * u[aidx];
        ssqs[1][aidx] += u[aidx] * u[aidx] + d[aidx] * d[aidx];
      }
      ++nm_ct;
    } else if (!hc_missing) {
      for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
        const uint64_t cur_dosage = 16384 * ((sd.a1 == aidx) + (sd.a2 == aidx));
        expected_sums[aidx] += cur_dosage;
        ssqs[0][aidx] += cur_dosage * cur_dosage;
        ssqs[1][aidx] += cur_dosage * cur_dosage;
      }
      extra_phased_het_ct += sd.phased;
      ++nm_ct;
    }
  }
  for (uint32_t is_minimac3_r2 = 0; is_minimac3_r2 != 2; ++is_minimac3_r2) {
    STD_ARRAY_DECL(uint32_t, 4, genocounts);
    std::vector<uint64_t> all_dosages(allele_ct);
    uint32_t het_ct = 0;
    double imp_r2 = 0.0;
    double* imp_r2_ptr = (allele_ct > 2)? (&imp_r2) : nullptr;
    const PglErr reterr = PgrGetMDCounts(sample_include, sample_include_interleaved_vec, pssi, subset_ct, vidx, is_minimac3_r2, &r.pgr, imp_r2_ptr, &het_ct, genocounts, all_dosages.data());
    CHECK(!reterr);
    if (reterr) {
      return;
    }
    int ok = (all_dosages == expected_sums) && (het_ct == expected_het_ct);
    for (uint32_t uii = 0; uii != 4; ++uii) {
      ok = ok && (genocounts[uii] == expected_genocounts[uii]);
    }
    if (imp_r2_ptr) {
      double expected_r2 = MultiallelicDiploidMinimac3R2(expected_sums.data(), ssqs[is_minimac3_r2].data(), nm_ct, allele_ct, is_minimac3_r2? extra_phased_het_ct : 0);
      if (!is_minimac3_r2) {
        expected_r2 *= 2;
      }
      // NaN when there are no nonmissing samples
      ok = ok && ((imp_r2 == expected_r2) || ((imp_r2 != imp_r2) && (expected_r2 != expected_r2)));
    }
    if ((!ok) && (g_mismatch_msg_ct++ < 10)) {
      fprintf(stderr, "PgrGetMDCounts mismatch: variant %u (allele_ct %u, minimac3 %u): dosage[0] %" PRIu64 " vs %" PRIu64 ", het %u vs %u, r2 %g vs %g\n", vidx, allele_ct, is_minimac3_r2, all_dosages[0], expected_sums[0], het_ct, expected_het_ct, imp_r2, imp_r2_ptr? MultiallelicDiploidMinimac3R2(expected_sums.data(), ssqs[is_minimac3_r2].data(), nm_ct, allele_ct, is_minimac3_r2? extra_phased_het_ct : 0) : 0.0);
    }
    CHECK(ok);
  }
}

// Checks PgrGetDCounts()'s ref vs. nonref collapse of a multiallelic
// variant, and that it refuses to compute imputation-r2 there.
static void CheckDCountsCollapse(Reader& r, const uintptr_t* sample_include, const uintptr_t* sample_include_interleaved_vec, PgrSampleSubsetIndex pssi, uint32_t subset_ct, uint32_t vidx, const std::vector<SampleData>& full) {
  uint64_t expected_dosages[2] = {0, 0};
  STD_ARRAY_DECL(uint32_t, 4, expected_genocounts);
  STD_ARRAY_FILL0(expected_genocounts);
  for (const SampleData& sd : full) {
    const uint32_t hc_missing = (sd.a1 == UINT32_MAX);
    const uint32_t nonref_ct = hc_missing? 0 : ((sd.a1 != 0) + (sd.a2 != 0));
    expected_genocounts[hc_missing? 3 : nonref_ct] += 1;
    if (sd.has_dosage) {
      expected_dosages[0] += 32768 - sd.dosage;
      expected_dosages[1] += sd.dosage;
    } else if (!hc_missing) {
      expected_dosages[0] += 16384 * (2 - nonref_ct);
      expected_dosages[1] += 16384 * nonref_ct;
    }
  }
  for (uint32_t is_minimac3_r2 = 0; is_minimac3_r2 != 2; ++is_minimac3_r2) {
    STD_ARRAY_DECL(uint32_t, 4, genocounts);
    uint64_t dosages[2];
    const PglErr reterr = PgrGetDCounts(sample_include, sample_include_interleaved_vec, pssi, subset_ct, vidx, is_minimac3_r2, &r.pgr, nullptr, genocounts, dosages);
    CHECK(!reterr);
    if (reterr) {
      return;
    }
    int ok = (dosages[0] == expected_dosages[0]) && (dosages[1] == expected_dosages[1]);
    for (uint32_t uii = 0; uii != 4; ++uii) {
      ok = ok && (genocounts[uii] == expected_genocounts[uii]);
    }
    if ((!ok) && (g_mismatch_msg_ct++ < 10)) {
      fprintf(stderr, "PgrGetDCounts mismatch: variant %u: dosages %" PRIu64 "/%" PRIu64 " vs %" PRIu64 "/%" PRIu64 "\n", vidx, dosages[0], dosages[1], expected_dosages[0], expected_dosages[1]);
    }
    CHECK(ok);
    double imp_r2;
    CHECK(PgrGetDCounts(sample_include, sample_include_interleaved_vec, pssi, subset_ct, vidx, is_minimac3_r2, &r.pgr, &imp_r2, genocounts, dosages) == kPglRetImproperFunctionCall);
  }
}

// Decodes one rare-alt track pair in PgrGetRaw() format; returns the end.
static const uintptr_t* DecodeRawRareAlt(const uintptr_t* iter, uint32_t raw_sample_ct, uintptr_t** present_ptr, unsigned char** cts_ptr, AlleleCode** codes_ptr, uint16_t** vals_ptr, uint32_t* sample_ct_ptr) {
  const uintptr_t entry_ct = iter[0];
  iter = &(iter[kWordsPerVec]);
  *present_ptr = const_cast<uintptr_t*>(iter);
  *sample_ct_ptr = PopcountWords(iter, BitCtToWordCt(raw_sample_ct));
  iter = &(iter[BitCtToAlignedWordCt(raw_sample_ct)]);
  *codes_ptr = R_CAST(AlleleCode*, const_cast<uintptr_t*>(iter));
  iter = &(iter[kWordsPerVec * DivUp(entry_ct, kBytesPerVec)]);
  *vals_ptr = R_CAST(uint16_t*, const_cast<uintptr_t*>(iter));
  iter = &(iter[kWordsPerVec * DivUp(2 * entry_ct, kBytesPerVec)]);
  *cts_ptr = R_CAST(unsigned char*, const_cast<uintptr_t*>(iter));
  return &(iter[kWordsPerVec * DivUp(*sample_ct_ptr, kBytesPerVec)]);
}

// Loads every variant with PgrGetRaw() and compares the decoded raw record
// (minus hardcall phase) against the source data.
static void CheckGetRaw(Reader& r, const std::vector<VariantData>& vds, const std::vector<uint32_t>& is_multi, uint32_t sample_ct) {
  const PgenGlobalFlags all_gflags = kfPgenGlobalMultiallelicHardcallFound | kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent;
  const uintptr_t buf_word_ct = 64 + 16 * BitCtToAlignedWordCt(sample_ct) + kWordsPerVec * DivUp(16 * sample_ct, kBytesPerVec) + 2 * PgrGetRawRareAltMaxWordCt(sample_ct, 255);
  uintptr_t* buf = AllocWords(buf_word_ct * kBitsPerWord);
  const uint32_t raw_sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  for (uint32_t vidx = 0; vidx != vds.size(); ++vidx) {
    const VariantData& vd = vds[vidx];
    const uint32_t multi_dosage = is_multi[vidx] && vd.pgv.dosage_ct;
    if (multi_dosage) {
      uintptr_t* loadbuf_iter = buf;
      CHECK(PgrGetRaw(vidx, all_gflags, &r.pgr, &loadbuf_iter, nullptr) == kPglRetNotYetSupported);
    }
    uintptr_t* loadbuf_iter = buf;
    unsigned char loaded_vrtype = 0;
    PglErr reterr = PgrGetRaw(vidx, all_gflags | kfPgenGlobalMultiallelicDosagePresent, &r.pgr, &loadbuf_iter, &loaded_vrtype);
    CHECK(!reterr);
    if (reterr) {
      continue;
    }
    CHECK(S_CAST(uintptr_t, loadbuf_iter - buf) <= buf_word_ct);
    PgenVariant pgv;
    memset(&pgv, 0, sizeof(PgenVariant));
    const uintptr_t* iter = buf;
    pgv.genovec = buf;
    iter = &(iter[NypCtToAlignedWordCt(sample_ct)]);
    if (loaded_vrtype & 8) {
      pgv.patch_01_ct = iter[0];
      pgv.patch_10_ct = iter[1];
      iter = &(iter[RoundUpPow2(2, kWordsPerVec)]);
      if (pgv.patch_01_ct) {
        pgv.patch_01_set = const_cast<uintptr_t*>(iter);
        iter = &(iter[raw_sample_ctl]);
        pgv.patch_01_vals = R_CAST(AlleleCode*, const_cast<uintptr_t*>(iter));
        iter = &(iter[RoundUpPow2(raw_sample_ctl + DivUp(pgv.patch_01_ct, kBytesPerWord), kWordsPerVec) - raw_sample_ctl]);
      }
      if (pgv.patch_10_ct) {
        pgv.patch_10_set = const_cast<uintptr_t*>(iter);
        iter = &(iter[raw_sample_ctl]);
        pgv.patch_10_vals = R_CAST(AlleleCode*, const_cast<uintptr_t*>(iter));
        iter = &(iter[RoundUpPow2(raw_sample_ctl + DivUp(2 * pgv.patch_10_ct, kBytesPerWord), kWordsPerVec) - raw_sample_ctl]);
      }
    }
    if (loaded_vrtype & 0x10) {
      const uint32_t het_ct = S_CAST(uint32_t, iter[0]);
#ifdef __LP64__
      const uint32_t explicit_phasepresent_ct = iter[0] >> 32;
#else
      const uint32_t explicit_phasepresent_ct = iter[1];
#endif
      const uint32_t phaseraw_word_ct = (8 / kBytesPerWord) + 1 + (het_ct / kBitsPerWord) + DivUp(explicit_phasepresent_ct, kBitsPerWord);
      iter = &(iter[RoundUpPow2(phaseraw_word_ct, kWordsPerVec)]);
    }
    if (loaded_vrtype & 0x60) {
      pgv.dosage_present = const_cast<uintptr_t*>(iter);
      pgv.dosage_ct = PopcountWords(iter, raw_sample_ctl);
      pgv.dosage_main = R_CAST(uint16_t*, const_cast<uintptr_t*>(&(iter[raw_sample_ctaw])));
      iter = &(iter[raw_sample_ctaw + kWordsPerVec * DivUp(pgv.dosage_ct, kBytesPerVec / 2)]);
      if (loaded_vrtype & 0x80) {
        pgv.dphase_present = const_cast<uintptr_t*>(iter);
        pgv.dphase_ct = PopcountWords(iter, raw_sample_ctl);
        pgv.dphase_delta = R_CAST(int16_t*, const_cast<uintptr_t*>(&(iter[raw_sample_ctaw])));
        iter = &(iter[raw_sample_ctaw + kWordsPerVec * DivUp(pgv.dphase_ct, kBytesPerVec / 2)]);
      }
      if (vd.allele_ct > 2) {
        iter = DecodeRawRareAlt(iter, sample_ct, &pgv.multidosage_present, &pgv.multidosage_cts, &pgv.multidosage_codes, &pgv.multidosage_vals, &pgv.multidosage_sample_ct);
        if (loaded_vrtype & 0x80) {
          uint16_t* deltas;
          iter = DecodeRawRareAlt(iter, sample_ct, &pgv.multidphase_present, &pgv.multidphase_cts, &pgv.multidphase_codes, &deltas, &pgv.multidphase_sample_ct);
          pgv.multidphase_delta = R_CAST(int16_t*, deltas);
        }
      }
    }
    CHECK(iter == loadbuf_iter);
    PgenVariant src = vd.pgv;
    if (!is_multi[vidx]) {
      src.patch_01_ct = 0;
      src.patch_10_ct = 0;
      src.phasepresent_ct = 0;
      src.dphase_ct = 0;
      src.multidosage_sample_ct = 0;
      src.multidphase_sample_ct = 0;
    }
    const uint32_t expected_vrtype = ((src.patch_01_ct || src.patch_10_ct)? 8 : 0) | (src.phasepresent_ct? 0x10 : 0) | (src.dosage_ct? 0x60 : 0) | (src.dphase_ct? 0x80 : 0);
    CHECK(loaded_vrtype == expected_vrtype);
    CheckSamples(Decode(src, sample_ct, 0, 1, 1), Decode(pgv, sample_ct, 0, 1, 1), "PgrGetRaw", vidx);
  }
  free(buf);
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
  uintptr_t* sample_include_interleaved_vec = nullptr;
  if (sample_include) {
    const uint32_t raw_sample_ctv = BitCtToVecCt(sample_ct);
    sample_include_interleaved_vec = AllocWords(raw_sample_ctv * kBitsPerVec);
    FillInterleavedMaskVec(sample_include, raw_sample_ctv, sample_include_interleaved_vec);
  }
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
      CheckMDCounts(r, sample_include, sample_include_interleaved_vec, pssi, subset_ct, vidx, 2, SubsetSamples(Decode(src, sample_ct, 1, 1, 1), sample_include));
      continue;
    }
    const std::vector<SampleData> full = SubsetSamples(Decode(vd.pgv, sample_ct, 1, 1, 1), sample_include);
    CheckMDCounts(r, sample_include, sample_include_interleaved_vec, pssi, subset_ct, vidx, allele_ct, full);
    if (vd.pgv.dosage_ct) {
      CheckDCountsCollapse(r, sample_include, sample_include_interleaved_vec, pssi, subset_ct, vidx, full);
    }
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
  free(sample_include_interleaved_vec);
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
  CheckGetRaw(r, vds, is_multi, sample_ct);
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

// ***** PgrValidate() *****

// Writes a .pvar (distinct made-up alleles) and .psam for fname_prefix.pgen,
// so that run_tests.sh can run plink2 --validate on it.
static void WritePvarPsam(const char* fname_prefix, const uintptr_t* allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct) {
  std::string fname = std::string(fname_prefix) + ".pvar";
  FILE* ff = fopen(fname.c_str(), "w");
  CHECK(ff);
  fputs("#CHROM\tPOS\tID\tREF\tALT\n", ff);
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    fprintf(ff, "1\t%u\tv%u\tA\t", vidx + 1, vidx);
    const uint32_t allele_ct = allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx];
    for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
      // "C" followed by aidx in base 4 (ACGT): distinct, and never "A"
      std::string allele = "C";
      for (uint32_t uii = aidx; uii; uii /= 4) {
        allele += "ACGT"[uii % 4];
      }
      fprintf(ff, "%s%s", (aidx == 1)? "" : ",", allele.c_str());
    }
    fputc('\n', ff);
  }
  fclose(ff);
  fname = std::string(fname_prefix) + ".psam";
  ff = fopen(fname.c_str(), "w");
  CHECK(ff);
  fputs("#IID\n", ff);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    fprintf(ff, "s%u\n", sample_idx);
  }
  fclose(ff);
}

static PglErr ValidateFile(const char* fname, const uintptr_t* allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, char* errstr) {
  errstr[0] = '\0';
  Reader r;
  PglErr reterr = OpenReader(fname, allele_idx_offsets, variant_ct, sample_ct, r);
  if (!reterr) {
    uintptr_t* workspace = AllocWords(64 * PgrValidateWorkspaceWordCt(sample_ct, r.pgfi.max_allele_ct));
    reterr = PgrValidate(&r.pgr, fname, nullptr, workspace, errstr);
    free(workspace);
  }
  CloseReader(r);
  return reterr;
}

// A 4-sample, single-variant record built from scratch.  By default, sample
// 0 has a missing hardcall and a 1.0 total alt dosage, and the other samples
// have hom-ref hardcalls and no dosage.
struct TinyVariant {
  uintptr_t* genovec;
  uintptr_t* patch_01_set;
  uintptr_t* patch_10_set;
  uintptr_t* dosage_present;
  uintptr_t* dphase_present;
  uintptr_t* md_present;
  uintptr_t* mdp_present;
  AlleleCode patch_01_vals[4];
  AlleleCode patch_10_vals[8];
  uint16_t dosage_main[4];
  int16_t dphase_delta[4];
  unsigned char md_cts[4];
  AlleleCode md_codes[16];
  uint16_t md_vals[16];
  unsigned char mdp_cts[4];
  AlleleCode mdp_codes[16];
  int16_t mdp_delta[16];
  PgenVariant pgv;
};

static void InitTiny(TinyVariant& tv) {
  memset(&tv.pgv, 0, sizeof(PgenVariant));
  uintptr_t** bufs[] = {&tv.genovec, &tv.patch_01_set, &tv.patch_10_set, &tv.dosage_present, &tv.dphase_present, &tv.md_present, &tv.mdp_present};
  for (uintptr_t** bp : bufs) {
    *bp = AllocWords(8);
  }
  tv.genovec[0] = 3;  // sample 0 missing, others hom-ref
  SetBit(0, tv.dosage_present);
  tv.dosage_main[0] = 16384;
  PgenVariant& pgv = tv.pgv;
  pgv.genovec = tv.genovec;
  pgv.patch_01_set = tv.patch_01_set;
  pgv.patch_01_vals = tv.patch_01_vals;
  pgv.patch_10_set = tv.patch_10_set;
  pgv.patch_10_vals = tv.patch_10_vals;
  pgv.dosage_present = tv.dosage_present;
  pgv.dosage_main = tv.dosage_main;
  pgv.dosage_ct = 1;
  pgv.dphase_present = tv.dphase_present;
  pgv.dphase_delta = tv.dphase_delta;
  pgv.multidosage_present = tv.md_present;
  pgv.multidosage_cts = tv.md_cts;
  pgv.multidosage_codes = tv.md_codes;
  pgv.multidosage_vals = tv.md_vals;
  pgv.multidphase_present = tv.mdp_present;
  pgv.multidphase_cts = tv.mdp_cts;
  pgv.multidphase_codes = tv.mdp_codes;
  pgv.multidphase_delta = tv.mdp_delta;
}

static void FreeTiny(TinyVariant& tv) {
  uintptr_t* bufs[] = {tv.genovec, tv.patch_01_set, tv.patch_10_set, tv.dosage_present, tv.dphase_present, tv.md_present, tv.mdp_present};
  for (uintptr_t* p : bufs) {
    free(p);
  }
}

// Sample 0 gets one rare-alt dosage.
static void TinyAddRare(TinyVariant& tv, uint32_t code, uint32_t val) {
  SetBit(0, tv.md_present);
  tv.md_cts[0] = 1;
  tv.md_codes[0] = code;
  tv.md_vals[0] = val;
  tv.pgv.multidosage_sample_ct = 1;
}

// Sample 0 gets an explicit dosage-phase value, and optionally one rare-alt
// dosage-phase value.
static void TinyAddDphase(TinyVariant& tv, int32_t delta, uint32_t rare_code, int32_t rare_delta) {
  SetBit(0, tv.dphase_present);
  tv.dphase_delta[0] = delta;
  tv.pgv.dphase_ct = 1;
  if (rare_code) {
    SetBit(0, tv.mdp_present);
    tv.mdp_cts[0] = 1;
    tv.mdp_codes[0] = rare_code;
    tv.mdp_delta[0] = rare_delta;
    tv.pgv.multidphase_sample_ct = 1;
  }
}

// Writes tv as a single-variant file, then validates it.  If expected_substr
// is nullptr, validation must succeed; otherwise it must fail with
// kPglRetMalformedInput and an error message containing expected_substr.
// If insert_hphase is set, a hardcall-phase track phasing sample 0 (which
// must be the only het) is spliced into the record after writing, since the
// writer won't produce some of the records under test.
static void CheckTiny(const char* desc, TinyVariant& tv, uint32_t allele_ct, uint32_t insert_hphase, const char* expected_substr) {
  const char* fname = "tmp_tiny.pgen";
  const uintptr_t allele_idx_offsets[2] = {0, allele_ct};
  STPgenWriter spgw;
  unsigned char* alloc;
  CHECK(!OpenWriter(fname, allele_idx_offsets, 1, 4, kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent | kfPgenGlobalMultiallelicDosagePresent, &spgw, &alloc));
  PglErr reterr = SpgwAppendMultiallelicGenovecDphase16(&tv.pgv, allele_ct, &spgw);
  CHECK(!reterr);
  CHECK(!SpgwFinish(&spgw));
  CleanupSpgw(&spgw, &reterr);
  free(alloc);
  if (insert_hphase) {
    std::vector<unsigned char> file = ReadFile(fname);
    // single vblock; one 8-bit vrtype, then the record length
    const uint32_t len_byte_ct = 1 + (file[11] & 3);
    uint64_t fpos;
    memcpy(&fpos, &file[12], 8);
    unsigned char& vrtype = file[20];
    // only a raw 1-byte main track is supported
    CHECK(((vrtype & 0x1f) == 0) && (len_byte_ct == 1));
    vrtype |= 0x10;
    file[21] += 1;
    // het_ct = 1: no explicit phasepresent, phaseinfo bit 0
    file.insert(file.begin() + fpos + 1, 0);
    FILE* ff = fopen(fname, "wb");
    CHECK(ff && (fwrite(file.data(), 1, file.size(), ff) == file.size()));
    fclose(ff);
  }
  char errstr[kPglErrstrBufBlen];
  reterr = ValidateFile(fname, allele_idx_offsets, 1, 4, errstr);
  if (!expected_substr) {
    if (reterr) {
      fprintf(stderr, "%s: unexpected validation failure: %s", desc, errstr);
    }
    CHECK(!reterr);
  } else {
    const int ok = (reterr == kPglRetMalformedInput) && strstr(errstr, expected_substr);
    if (!ok) {
      fprintf(stderr, "%s: expected failure containing \"%s\", got %d: %s\n", desc, expected_substr, S_CAST(int, reterr), errstr);
    }
    CHECK(ok);
  }
}

static void ValidatorTinyTests() {
  TinyVariant tv;
  // Sanity: missing hardcall, rare-alt dosage.
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  CheckTiny("valid baseline", tv, 3, 0, nullptr);
  FreeTiny(tv);

  // Rare-alt dosages exceeding the total alt dosage.
  InitTiny(tv);
  TinyAddRare(tv, 2, 16385);
  CheckTiny("rare sum > total", tv, 3, 0, "sum to more than its total alt dosage");
  FreeTiny(tv);

  // Zero rare-alt dosage (zeros must be omitted).
  InitTiny(tv);
  TinyAddRare(tv, 2, 0);
  CheckTiny("zero rare dosage", tv, 3, 0, "rare-alt dosage is zero or greater than 2");
  FreeTiny(tv);

  // Rare-alt dosage for a sample with no dosage.
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  ClearBit(0, tv.md_present);
  SetBit(1, tv.md_present);
  CheckTiny("rare dosage without dosage", tv, 3, 0, "which has no (non-missing) dosage");
  FreeTiny(tv);

  // Same, for a sample before one that has a dosage.
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  ClearBit(0, tv.dosage_present);
  SetBit(1, tv.dosage_present);
  tv.dosage_main[0] = 0;  // sample 1, hom-ref
  CheckTiny("rare dosage without dosage (before)", tv, 3, 0, "which has no (non-missing) dosage");
  FreeTiny(tv);

  // Hardcall 0/2 (sample 0 het, patch 01 -> alt2), total alt dosage 1.0.
  // All of it on alt1 is inconsistent; all on alt2 is fine.
  for (uint32_t on_alt2 = 0; on_alt2 != 2; ++on_alt2) {
    InitTiny(tv);
    tv.genovec[0] = 1;
    SetBit(0, tv.patch_01_set);
    tv.patch_01_vals[0] = 2;
    tv.pgv.patch_01_ct = 1;
    if (on_alt2) {
      TinyAddRare(tv, 2, 16384);
    }
    CheckTiny(on_alt2? "0/2 hardcall, alt2 dosage" : "0/2 hardcall, alt1 dosage", tv, 3, 0, on_alt2? nullptr : "inconsistent with its (non-missing) hardcall");
    FreeTiny(tv);
  }
  // Hardcall 2/3 (patch 10), total alt dosage 2.0: split 1.0/1.0 between
  // alt2 and alt3 is consistent; all on alt1 isn't.
  for (uint32_t split = 0; split != 2; ++split) {
    InitTiny(tv);
    tv.genovec[0] = 2;
    SetBit(0, tv.patch_10_set);
    tv.patch_10_vals[0] = 2;
    tv.patch_10_vals[1] = 3;
    tv.pgv.patch_10_ct = 1;
    tv.dosage_main[0] = 32768;
    if (split) {
      SetBit(0, tv.md_present);
      tv.md_cts[0] = 2;
      tv.md_codes[0] = 2;
      tv.md_codes[1] = 3;
      tv.md_vals[0] = 16384;
      tv.md_vals[1] = 16384;
      tv.pgv.multidosage_sample_ct = 1;
    }
    CheckTiny(split? "2/3 hardcall, 2/3 dosage" : "2/3 hardcall, alt1 dosage", tv, 4, 0, split? nullptr : "inconsistent with its (non-missing) hardcall");
    FreeTiny(tv);
  }

  // Phased rare-alt dosages.
  // Valid: alt2 dosage 0.5, all on haplotype 1.
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  TinyAddDphase(tv, 0, 2, 8192);
  CheckTiny("valid phased rare alt", tv, 4, 0, nullptr);
  FreeTiny(tv);
  // Rare-alt dphase entry without a matching rare-alt dosage.
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  TinyAddDphase(tv, 0, 3, 4096);
  CheckTiny("dphase not a subset", tv, 4, 0, "has no corresponding rare-alt dosage");
  FreeTiny(tv);
  // Same, with the unmatched allele between two that have dosages.
  InitTiny(tv);
  SetBit(0, tv.md_present);
  tv.md_cts[0] = 2;
  tv.md_codes[0] = 2;
  tv.md_codes[1] = 4;
  tv.md_vals[0] = 4096;
  tv.md_vals[1] = 4096;
  tv.pgv.multidosage_sample_ct = 1;
  TinyAddDphase(tv, 0, 3, 2048);
  CheckTiny("dphase not a subset (interior)", tv, 5, 0, "has no corresponding rare-alt dosage");
  FreeTiny(tv);
  // Zero rare-alt delta (zeros must be omitted).
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  TinyAddDphase(tv, 0, 2, 0);
  CheckTiny("zero rare delta", tv, 4, 0, "inconsistent with its dosage");
  FreeTiny(tv);
  // |rare delta| > rare dosage.
  InitTiny(tv);
  TinyAddRare(tv, 2, 4096);
  TinyAddDphase(tv, 0, 2, 8192);
  CheckTiny("rare delta too large", tv, 4, 0, "inconsistent with its dosage");
  FreeTiny(tv);
  // Haplotype 2 has no alt allele (total alt 1.0 = 1|0), but its alt2 share
  // would be 0.5.
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  TinyAddDphase(tv, 16384, 2, -8192);
  CheckTiny("negative alt1 haplotype share", tv, 4, 0, "negative alt1 dosage on one haplotype");
  FreeTiny(tv);
  // Rare-alt dphase entry for a sample without explicit dosage phase.
  InitTiny(tv);
  TinyAddRare(tv, 2, 8192);
  SetBit(1, tv.dosage_present);
  tv.dosage_main[1] = 0;
  tv.pgv.dosage_ct = 2;
  TinyAddDphase(tv, 0, 0, 0);
  ClearBit(0, tv.dphase_present);
  SetBit(1, tv.dphase_present);
  SetBit(0, tv.mdp_present);
  tv.mdp_cts[0] = 1;
  tv.mdp_codes[0] = 2;
  tv.mdp_delta[0] = 4096;
  tv.pgv.multidphase_sample_ct = 1;
  CheckTiny("rare dphase without dphase", tv, 4, 0, "which has no (non-missing) explicit dosage-phase value");
  FreeTiny(tv);

  // Implicit dosage-phase: sample 0 is a 0/1 het with a dosage.  Phasing it
  // (spliced in after writing) requires an explicit dosage-phase value when
  // there are more than two alleles.
  // The invalid case runs last, so that run_tests.sh can check that plink2
  // --validate rejects tmp_tiny.pgen.
  for (uint32_t explicit_dphase = 2; explicit_dphase--; ) {
    InitTiny(tv);
    tv.genovec[0] = 1 | (2 << 4) | (3 << 6);  // het, hom-ref, hom-alt1, missing
    tv.dosage_main[0] = 16384;
    if (explicit_dphase) {
      TinyAddDphase(tv, 16384, 0, 0);
    }
    CheckTiny(explicit_dphase? "phased het, explicit dphase" : "phased het, no dphase", tv, 3, 1, explicit_dphase? nullptr : "no explicit dosage-phase value");
    FreeTiny(tv);
  }
  const uintptr_t tiny_allele_idx_offsets[2] = {0, 3};
  WritePvarPsam("tmp_tiny", tiny_allele_idx_offsets, 1, 4);
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
  {
    std::string prefix(fname);
    prefix.resize(prefix.size() - 5);  // remove ".pgen"
    WritePvarPsam(prefix.c_str(), allele_idx_offsets.data(), variant_ct, sample_ct);
    // Expected per-variant count of samples with a dosage or nonmissing
    // hardcall, for run_tests.sh to compare against plink2 --missing/--freq.
    FILE* ff = fopen((prefix + "_expected_nm.txt").c_str(), "w");
    CHECK(ff);
    for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
      PgenVariant src = vds[vidx].pgv;
      if (!is_multi[vidx]) {
        src.patch_01_ct = 0;
        src.patch_10_ct = 0;
      }
      uint32_t nm_ct = 0;
      for (const SampleData& sd : Decode(src, sample_ct, 0, 0, 0)) {
        nm_ct += sd.has_dosage || (sd.a1 != UINT32_MAX);
      }
      fprintf(ff, "v%u\t%u\n", vidx, nm_ct);
    }
    fclose(ff);
  }
  char errstr[kPglErrstrBufBlen];
  const PglErr validate_reterr = ValidateFile(fname, allele_idx_offsets.data(), variant_ct, sample_ct, errstr);
  if (validate_reterr) {
    fprintf(stderr, "%s: PgrValidate() failed: %s", fname, errstr);
  }
  CHECK(!validate_reterr);
  if (check_corrupt) {
    CheckCorruptValue(fname, vds, is_multi, record_ends, allele_idx_offsets.data(), sample_ct, 0);
    CheckCorruptValue(fname, vds, is_multi, record_ends, allele_idx_offsets.data(), sample_ct, 32769);
  }
  for (VariantData& vd : vds) {
    FreeVariant(vd);
  }
}

// ***** "compare" mode: semantic comparison of two filesets *****

struct PvarInfo {
  std::vector<std::string> ids;
  std::vector<std::vector<std::string>> alleles;
  std::vector<uintptr_t> allele_idx_offsets;
};

static std::vector<std::string> SplitTabs(const std::string& line) {
  std::vector<std::string> out;
  size_t start = 0;
  while (1) {
    const size_t pos = line.find('\t', start);
    if (pos == std::string::npos) {
      out.push_back(line.substr(start));
      return out;
    }
    out.push_back(line.substr(start, pos - start));
    start = pos + 1;
  }
}

static std::vector<std::string> ReadLines(const std::string& fname) {
  std::vector<std::string> out;
  FILE* ff = fopen(fname.c_str(), "r");
  if (!ff) {
    return out;
  }
  char buf[1 << 16];
  while (fgets(buf, sizeof(buf), ff)) {
    std::string line(buf);
    while ((!line.empty()) && ((line.back() == '\n') || (line.back() == '\r'))) {
      line.pop_back();
    }
    out.push_back(line);
  }
  fclose(ff);
  return out;
}

static int ReadPvar(const std::string& prefix, PvarInfo& pi) {
  int id_col = -1;
  int ref_col = -1;
  int alt_col = -1;
  pi.allele_idx_offsets.push_back(0);
  for (const std::string& line : ReadLines(prefix + ".pvar")) {
    if ((line.size() >= 2) && (line[0] == '#') && (line[1] == '#')) {
      continue;
    }
    const std::vector<std::string> toks = SplitTabs(line);
    if (line[0] == '#') {
      for (size_t i = 0; i != toks.size(); ++i) {
        if (toks[i] == "ID") {
          id_col = i;
        } else if (toks[i] == "REF") {
          ref_col = i;
        } else if (toks[i] == "ALT") {
          alt_col = i;
        }
      }
      continue;
    }
    if ((id_col < 0) || (ref_col < 0) || (alt_col < 0)) {
      return 1;
    }
    pi.ids.push_back(toks[id_col]);
    std::vector<std::string> alleles(1, toks[ref_col]);
    size_t start = 0;
    const std::string& alt = toks[alt_col];
    while (1) {
      const size_t pos = alt.find(',', start);
      alleles.push_back(alt.substr(start, (pos == std::string::npos)? std::string::npos : (pos - start)));
      if (pos == std::string::npos) {
        break;
      }
      start = pos + 1;
    }
    pi.allele_idx_offsets.push_back(pi.allele_idx_offsets.back() + alleles.size());
    pi.alleles.push_back(alleles);
  }
  return pi.ids.empty();
}

static std::vector<std::string> ReadPsamIids(const std::string& prefix) {
  std::vector<std::string> out;
  int iid_col = -1;
  for (const std::string& line : ReadLines(prefix + ".psam")) {
    const std::vector<std::string> toks = SplitTabs(line);
    if (line[0] == '#') {
      for (size_t i = 0; i != toks.size(); ++i) {
        if ((toks[i] == "IID") || (toks[i] == "#IID")) {
          iid_col = i;
        }
      }
      continue;
    }
    out.push_back(toks[iid_col]);
  }
  return out;
}

// Per-allele dosage and dosage-phase vectors (empty when absent), from a
// decoded sample.
static void DenseAlleleVals(const SampleData& sd, uint32_t allele_ct, std::vector<int32_t>& u, std::vector<int32_t>& d) {
  u.clear();
  d.clear();
  if (sd.has_dosage) {
    u.assign(allele_ct, 0);
    u[0] = 32768 - sd.dosage;
    u[1] = sd.dosage;
    for (const auto& pr : sd.md) {
      u[pr.first] = pr.second;
      u[1] -= pr.second;
    }
  }
  if (sd.has_dphase) {
    d.assign(allele_ct, 0);
    d[0] = -sd.dphase;
    d[1] = sd.dphase;
    for (const auto& pr : sd.mdp) {
      d[pr.first] = pr.second;
      d[1] -= pr.second;
    }
  }
}

// Expected effect of --hard-call-threshold (hc_thresh >= 0, in 16384ths of a
// copy), 'fill-missing-from-dosage', and --dosage-erase-threshold
// (erase_thresh >= 0) on one sample, per the plink2_help.cc definitions.
// Written independently of plink2's implementation; e.g. the nearest hardcall
// is found by brute force over all genotypes.
// flip_vars/flip_samples: expected --flip + --flip-subset effect; the
// flip_samples' alleles of flip_vars variants are complemented (A<->T,
// C<->G, N unchanged), except in variants with an allele outside {A, C, G,
// T, N}, which are skipped (cf. 'permissive').  This is applied after the
// operations above.
struct CallOps {
  int32_t hc_thresh = -1;
  uint32_t fill = 0;
  int32_t erase_thresh = -1;
  std::set<std::string> flip_vars;
  std::set<std::string> flip_samples;
};

// Returns the strand-flipped allele, or the empty string if allele isn't one
// of A, C, G, T, N.
static std::string FlipAllele(const std::string& allele) {
  if (allele.size() != 1) {
    return std::string();
  }
  switch (allele[0]) {
  case 'A':
    return "T";
  case 'C':
    return "G";
  case 'G':
    return "C";
  case 'T':
    return "A";
  case 'N':
    return "N";
  }
  return std::string();
}

static void ApplyExpectedCallOps(const CallOps& ops, uint32_t allele_ct, SampleData& sd) {
  if (!sd.has_dosage) {
    return;
  }
  std::vector<int32_t> u;
  std::vector<int32_t> d;
  DenseAlleleVals(sd, allele_ct, u, d);
  if (d.empty()) {
    d.assign(allele_ct, 0);
  }
  uint32_t a1 = sd.a1;
  uint32_t a2 = sd.a2;
  if (ops.hc_thresh >= 0) {
    int64_t dist_x2 = 0;
    std::vector<uint32_t> copies;
    for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
      const int32_t r = (u[aidx] + 8192) / 16384;
      dist_x2 += abs(u[aidx] - 16384 * r);
      for (int32_t k = 0; k != r; ++k) {
        copies.push_back(aidx);
      }
    }
    if ((dist_x2 <= 2 * ops.hc_thresh) && (copies.size() == 2)) {
      a1 = copies[0];
      a2 = copies[1];
    } else {
      a1 = UINT32_MAX;
      a2 = UINT32_MAX;
    }
  }
  if (ops.fill && (a1 == UINT32_MAX)) {
    // squared distance * 16384^2, up to a constant: sum_i (16384 h_i)^2 -
    // 2 * 16384 h_i u_i
    int64_t best_score = INT64_MAX;
    uint32_t best1 = 0;
    uint32_t best2 = 0;
    for (uint32_t c1 = 0; c1 != allele_ct; ++c1) {
      for (uint32_t c2 = c1; c2 != allele_ct; ++c2) {
        // only c1 and c2 have nonzero h
        const int64_t score = (c1 == c2)? (4LL * 16384 * 16384 - 2LL * 32768 * u[c1]) : (2LL * 16384 * 16384 - 2LL * 16384 * (u[c1] + u[c2]));
        if (score < best_score) {
          best_score = score;
          best1 = c1;
          best2 = c2;
        }
      }
    }
    int64_t l1 = 0;
    for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
      l1 += abs(u[aidx] - 16384 * S_CAST(int32_t, (aidx == best1) + (aidx == best2)));
    }
    if (l1 <= 16384) {
      a1 = best1;
      a2 = best2;
    }
  }
  if ((a1 != sd.a1) || (a2 != sd.a2)) {
    sd.a1 = a1;
    sd.a2 = a2;
    sd.phased = 0;
    sd.phaseinfo = 0;
    if (sd.has_dphase && (a1 != UINT32_MAX) && (a1 != a2)) {
      const int32_t delta_diff = d[a2] - d[a1];
      if (abs(delta_diff) > 16384) {
        sd.phased = 1;
        sd.phaseinfo = (delta_diff > 0);
      }
    }
  }
  if (ops.erase_thresh >= 0) {
    int32_t erase;
    if (!sd.has_dphase) {
      int64_t dist_x2 = 0;
      for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
        dist_x2 += abs(u[aidx] - 16384 * ((u[aidx] + 8192) / 16384));
      }
      erase = (dist_x2 <= 2 * ops.erase_thresh);
    } else {
      // haplotype dosages in 32768ths of a copy; distance in 16384ths is
      // 0.25 * sum_i |hap_i - 32768 * round|
      int64_t dist1 = 0;
      int64_t dist2 = 0;
      for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
        const int32_t hap1 = u[aidx] + d[aidx];
        const int32_t hap2 = u[aidx] - d[aidx];
        dist1 += abs(hap1 - 32768 * ((hap1 + 16384) / 32768));
        dist2 += abs(hap2 - 32768 * ((hap2 + 16384) / 32768));
      }
      erase = (dist1 <= 2 * ops.erase_thresh) && (dist2 <= 2 * ops.erase_thresh);
    }
    if (erase) {
      sd.has_dosage = 0;
      sd.dosage = 0;
      sd.md.clear();
      sd.has_dphase = 0;
      sd.dphase = 0;
      sd.mdp.clear();
    }
  }
}

// Compares orig_prefix.{pgen,pvar,psam} and new_prefix.*, matching samples by
// IID, variants by ID, and alleles by string.  Hardcalls, phase (including
// haplotype order), dosages, and dosage-phase values must agree.  If strip_phase is set, the new fileset must have no
// phase information; if strip_dosage is set, no dosages.
// If hets_missing is set, the original's heterozygous hardcalls (including
// altx/alty) must be missing in the new fileset, with their dosages erased,
// and all phase information must be gone (cf. --set-mixed-mt-missing).
// ops: expected --hard-call-threshold etc. effects, applied to the original
// (in its allele order) before comparing.
static int CompareFilesets(const char* orig_prefix, const char* new_prefix, uint32_t strip_phase, uint32_t strip_dosage, uint32_t hets_missing, const CallOps& ops) {
  strip_phase |= hets_missing;
  PvarInfo orig_pi;
  PvarInfo new_pi;
  if (ReadPvar(orig_prefix, orig_pi) || ReadPvar(new_prefix, new_pi)) {
    fprintf(stderr, "compare: failed to read .pvar files\n");
    return 1;
  }
  const std::vector<std::string> orig_iids = ReadPsamIids(orig_prefix);
  const std::vector<std::string> new_iids = ReadPsamIids(new_prefix);
  std::vector<uint32_t> new_to_orig_sample(new_iids.size());
  for (size_t new_idx = 0; new_idx != new_iids.size(); ++new_idx) {
    size_t orig_idx = 0;
    while ((orig_idx != orig_iids.size()) && (orig_iids[orig_idx] != new_iids[new_idx])) {
      ++orig_idx;
    }
    if (orig_idx == orig_iids.size()) {
      fprintf(stderr, "compare: sample %s not in original\n", new_iids[new_idx].c_str());
      return 1;
    }
    new_to_orig_sample[new_idx] = orig_idx;
  }
  const uint32_t orig_sample_ct = orig_iids.size();
  const uint32_t new_sample_ct = new_iids.size();
  const uint32_t orig_variant_ct = orig_pi.ids.size();
  const uint32_t new_variant_ct = new_pi.ids.size();
  Reader orig_r;
  Reader new_r;
  if (OpenReader((std::string(orig_prefix) + ".pgen").c_str(), orig_pi.allele_idx_offsets.data(), orig_variant_ct, orig_sample_ct, orig_r) ||
      OpenReader((std::string(new_prefix) + ".pgen").c_str(), new_pi.allele_idx_offsets.data(), new_variant_ct, new_sample_ct, new_r)) {
    fprintf(stderr, "compare: failed to open .pgen files\n");
    return 1;
  }
  ReaderBufs orig_rb;
  ReaderBufs new_rb;
  InitReaderBufs(orig_sample_ct, 255, orig_rb);
  InitReaderBufs(new_sample_ct, 255, new_rb);
  PgrSampleSubsetIndex orig_pssi;
  PgrSampleSubsetIndex new_pssi;
  PgrClearSampleSubsetIndex(&orig_r.pgr, &orig_pssi);
  PgrClearSampleSubsetIndex(&new_r.pgr, &new_pssi);
  uint64_t mismatch_ct = 0;
  uint64_t md_entry_ct = 0;
  for (uint32_t new_vidx = 0; new_vidx != new_variant_ct; ++new_vidx) {
    uint32_t orig_vidx = 0;
    while ((orig_vidx != orig_variant_ct) && (orig_pi.ids[orig_vidx] != new_pi.ids[new_vidx])) {
      ++orig_vidx;
    }
    if (orig_vidx == orig_variant_ct) {
      fprintf(stderr, "compare: variant %s not in original\n", new_pi.ids[new_vidx].c_str());
      return 1;
    }
    const std::vector<std::string>& orig_alleles = orig_pi.alleles[orig_vidx];
    const std::vector<std::string>& new_alleles = new_pi.alleles[new_vidx];
    // orig allele index -> new allele index (UINT32_MAX if dropped)
    std::vector<uint32_t> orig_to_new(orig_alleles.size(), UINT32_MAX);
    for (size_t new_aidx = 0; new_aidx != new_alleles.size(); ++new_aidx) {
      for (size_t orig_aidx = 0; orig_aidx != orig_alleles.size(); ++orig_aidx) {
        if (orig_alleles[orig_aidx] == new_alleles[new_aidx]) {
          orig_to_new[orig_aidx] = new_aidx;
        }
      }
    }
    // same, for --flip-subset samples (empty if the variant isn't flipped)
    std::vector<uint32_t> orig_to_new_flipped;
    if (ops.flip_vars.count(orig_pi.ids[orig_vidx])) {
      orig_to_new_flipped.assign(orig_alleles.size(), UINT32_MAX);
      for (size_t orig_aidx = 0; orig_aidx != orig_alleles.size(); ++orig_aidx) {
        const std::string flipped = FlipAllele(orig_alleles[orig_aidx]);
        if (flipped.empty()) {
          orig_to_new_flipped.clear();
          break;
        }
        for (size_t new_aidx = 0; new_aidx != new_alleles.size(); ++new_aidx) {
          if (new_alleles[new_aidx] == flipped) {
            orig_to_new_flipped[orig_aidx] = new_aidx;
          }
        }
      }
    }
    if (PgrGetMDp(nullptr, orig_pssi, orig_sample_ct, orig_vidx, &orig_r.pgr, &orig_rb.pgv) ||
        PgrGetMDp(nullptr, new_pssi, new_sample_ct, new_vidx, &new_r.pgr, &new_rb.pgv)) {
      fprintf(stderr, "compare: read failure at variant %s\n", new_pi.ids[new_vidx].c_str());
      return 1;
    }
    const std::vector<SampleData> orig_sds = Decode(orig_rb.pgv, orig_sample_ct, 1, 1, 1);
    const std::vector<SampleData> new_sds = Decode(new_rb.pgv, new_sample_ct, 1, 1, 1);
    md_entry_ct += new_rb.pgv.multidosage_sample_ct;
    for (uint32_t new_sidx = 0; new_sidx != new_sample_ct; ++new_sidx) {
      SampleData o = orig_sds[new_to_orig_sample[new_sidx]];
      const SampleData& n = new_sds[new_sidx];
      const std::vector<uint32_t>& allele_map = ((!orig_to_new_flipped.empty()) && ops.flip_samples.count(new_iids[new_sidx]))? orig_to_new_flipped : orig_to_new;
      ApplyExpectedCallOps(ops, orig_alleles.size(), o);
      if (hets_missing && (o.a1 != UINT32_MAX) && (o.a1 != o.a2)) {
        o.a1 = UINT32_MAX;
        o.a2 = UINT32_MAX;
        o.has_dosage = 0;
        o.has_dphase = 0;
        o.md.clear();
        o.mdp.clear();
      }
      // hardcalls, as unordered new-space allele pairs
      uint32_t expected_a1 = UINT32_MAX;
      uint32_t expected_a2 = UINT32_MAX;
      if (o.a1 != UINT32_MAX) {
        expected_a1 = allele_map[o.a1];
        expected_a2 = allele_map[o.a2];
        if (expected_a1 > expected_a2) {
          std::swap(expected_a1, expected_a2);
        }
      }
      int ok = (n.a1 == expected_a1) && (n.a2 == expected_a2);
      // phase
      if (strip_phase) {
        ok = ok && (!n.phased) && (!n.has_dphase);
      } else {
        ok = ok && (n.phased == o.phased);
        if (o.phased) {
          // phaseinfo is relative to allele order
          ok = ok && (n.phaseinfo == (o.phaseinfo ^ (allele_map[o.a1] > allele_map[o.a2])));
        }
      }
      // dosages
      std::vector<int32_t> ou;
      std::vector<int32_t> od;
      std::vector<int32_t> nu;
      std::vector<int32_t> nd;
      DenseAlleleVals(o, orig_alleles.size(), ou, od);
      DenseAlleleVals(n, new_alleles.size(), nu, nd);
      if (strip_dosage) {
        ok = ok && nu.empty() && nd.empty();
      } else {
        std::vector<int32_t> eu;
        std::vector<int32_t> ed;
        if (!ou.empty()) {
          eu.assign(new_alleles.size(), 0);
          for (size_t orig_aidx = 0; orig_aidx != ou.size(); ++orig_aidx) {
            if (allele_map[orig_aidx] != UINT32_MAX) {
              eu[allele_map[orig_aidx]] = ou[orig_aidx];
            } else if (ou[orig_aidx]) {
              ok = 0;
            }
          }
        }
        if ((!od.empty()) && (!strip_phase)) {
          ed.assign(new_alleles.size(), 0);
          for (size_t orig_aidx = 0; orig_aidx != od.size(); ++orig_aidx) {
            if (allele_map[orig_aidx] != UINT32_MAX) {
              ed[allele_map[orig_aidx]] = od[orig_aidx];
            }
          }
        }
        ok = ok && (nu == eu);
        if (!strip_phase) {
          ok = ok && (nd == ed);
        }
      }
      if (!ok) {
        if (mismatch_ct < 10) {
          fprintf(stderr, "compare: mismatch at variant %s, sample %s: hardcall %u/%u vs expected %u/%u, phased %u vs %u, dosage %s, dphase %s\n", new_pi.ids[new_vidx].c_str(), new_iids[new_sidx].c_str(), n.a1, n.a2, expected_a1, expected_a2, n.phased, o.phased, (nu.empty() == ou.empty())? "presence ok" : "presence differs", (nd.empty() == od.empty())? "presence ok" : "presence differs");
        }
        ++mismatch_ct;
      }
    }
  }
  FreeReaderBufs(orig_rb);
  FreeReaderBufs(new_rb);
  CloseReader(orig_r);
  CloseReader(new_r);
  if (mismatch_ct) {
    fprintf(stderr, "compare: %" PRIu64 " mismatch(es).\n", mismatch_ct);
    return 1;
  }
  printf("compare: %s matches %s (%u variants, %u samples, %" PRIu64 " samples with rare-alt dosages).\n", new_prefix, orig_prefix, new_variant_ct, new_sample_ct, md_entry_ct);
  return 0;
}

// Writes tmp_fixed.{pgen,pvar,psam}: fixed-width phased-dosage storage (mode
// 0x04) for 50 samples and 6 variants, which the .pvar declares to have 3
// alleles each.  Fixed-width records have no rare-alt tracks, so the readers
// must not look for any; rare-alt dosages are all zero.
static void CheckFixedWidth() {
  const uint32_t sample_ct = 50;
  const uint32_t variant_ct = 6;
  std::vector<uintptr_t> allele_idx_offsets(variant_ct + 1);
  for (uint32_t vidx = 0; vidx <= variant_ct; ++vidx) {
    allele_idx_offsets[vidx] = 3 * vidx;
  }
  std::vector<unsigned char> file = {0x6c, 0x1b, 0x04};
  for (uint32_t uii = 0; uii != 4; ++uii) {
    file.push_back((variant_ct >> (8 * uii)) & 255);
  }
  for (uint32_t uii = 0; uii != 4; ++uii) {
    file.push_back((sample_ct >> (8 * uii)) & 255);
  }
  file.push_back(0);
  std::vector<std::vector<uint32_t>> genos(variant_ct, std::vector<uint32_t>(sample_ct));
  std::vector<std::vector<uint16_t>> dosages(variant_ct, std::vector<uint16_t>(sample_ct));
  std::vector<std::vector<int16_t>> deltas(variant_ct, std::vector<int16_t>(sample_ct));
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    std::vector<unsigned char> genobytes(DivUp(sample_ct, 4), 0);
    for (uint32_t s = 0; s != sample_ct; ++s) {
      uint32_t geno = RandBelow(4);
      uint16_t dosage = 65535;
      int16_t delta = -32768;
      if (geno != 3) {
        // within half a copy of the hardcall
        const int32_t lo = MAXV(0, S_CAST(int32_t, 16384 * geno) - 8192);
        const int32_t hi = MINV(32768, S_CAST(int32_t, 16384 * geno) + 8192);
        dosage = lo + RandBelow(hi - lo + 1);
        const int32_t lim = MINV(dosage, 32768 - dosage);
        delta = S_CAST(int32_t, RandBelow(2 * lim + 1)) - lim;
      } else if (RandBelow(2)) {
        // missing hardcall, nonmissing dosage
        dosage = RandBelow(32769);
        const int32_t lim = MINV(dosage, 32768 - dosage);
        delta = S_CAST(int32_t, RandBelow(2 * lim + 1)) - lim;
      }
      genos[vidx][s] = geno;
      dosages[vidx][s] = dosage;
      deltas[vidx][s] = delta;
      genobytes[s / 4] |= geno << (2 * (s % 4));
    }
    file.insert(file.end(), genobytes.begin(), genobytes.end());
    for (uint16_t d : dosages[vidx]) {
      file.push_back(d & 255);
      file.push_back(d >> 8);
    }
    for (int16_t d : deltas[vidx]) {
      file.push_back(S_CAST(uint16_t, d) & 255);
      file.push_back(S_CAST(uint16_t, d) >> 8);
    }
  }
  FILE* ff = fopen("tmp_fixed.pgen", "wb");
  CHECK(ff && (fwrite(file.data(), 1, file.size(), ff) == file.size()));
  fclose(ff);
  WritePvarPsam("tmp_fixed", allele_idx_offsets.data(), variant_ct, sample_ct);

  char errstr[kPglErrstrBufBlen];
  const PglErr validate_reterr = ValidateFile("tmp_fixed.pgen", allele_idx_offsets.data(), variant_ct, sample_ct, errstr);
  if (validate_reterr) {
    fprintf(stderr, "tmp_fixed.pgen: PgrValidate() failed: %s", errstr);
  }
  CHECK(!validate_reterr);
  Reader r;
  CHECK(!OpenReader("tmp_fixed.pgen", allele_idx_offsets.data(), variant_ct, sample_ct, r));
  ReaderBufs rb;
  InitReaderBufs(sample_ct, 3, rb);
  PgrSampleSubsetIndex pssi;
  PgrClearSampleSubsetIndex(&r.pgr, &pssi);
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    PgenVariant& out = rb.pgv;
    CHECK(!PgrGetMDp(nullptr, pssi, sample_ct, vidx, &r.pgr, &out));
    CHECK(!out.multidosage_sample_ct);
    CHECK(!out.multidphase_sample_ct);
    const std::vector<SampleData> sds = Decode(out, sample_ct, 1, 1, 1);
    for (uint32_t s = 0; s != sample_ct; ++s) {
      const SampleData& sd = sds[s];
      const uint32_t geno = genos[vidx][s];
      CHECK((geno == 3)? (sd.a1 == UINT32_MAX) : ((sd.a1 == (geno == 2)) && (sd.a2 == (geno != 0))));
      CHECK(sd.has_dosage == (dosages[vidx][s] != 65535));
      if (sd.has_dosage) {
        CHECK(sd.dosage == dosages[vidx][s]);
        CHECK(sd.has_dphase && (sd.dphase == deltas[vidx][s]));
      }
    }
    // rare-alt dosages are zero
    CHECK(!PgrGet1D(nullptr, pssi, sample_ct, vidx, 2, &r.pgr, out.genovec, out.dosage_present, out.dosage_main, &out.dosage_ct));
    for (uint32_t dosage_idx = 0; dosage_idx != out.dosage_ct; ++dosage_idx) {
      CHECK(!out.dosage_main[dosage_idx]);
    }
    STD_ARRAY_DECL(uint32_t, 4, genocounts);
    uint64_t all_dosages[3];
    uint32_t het_ct;
    double imp_r2;
    CHECK(!PgrGetMDCounts(nullptr, nullptr, pssi, sample_ct, vidx, 1, &r.pgr, &imp_r2, &het_ct, genocounts, all_dosages));
    CHECK(!all_dosages[2]);
  }
  FreeReaderBufs(rb);
  CloseReader(r);
}

// Writes tmp_flipfail.{pgen,pvar,psam}: one 3-allele variant, where sample 0
// has a missing hardcall and dosages of 1 copy each for REF and ALT1; the
// other samples are hom-ref without dosages.  (run_tests.sh replaces the
// .pvar alleles.)
static void WriteFlipFail() {
  TinyVariant tv;
  InitTiny(tv);
  const uintptr_t allele_idx_offsets[2] = {0, 3};
  STPgenWriter spgw;
  unsigned char* alloc;
  CHECK(!OpenWriter("tmp_flipfail.pgen", allele_idx_offsets, 1, 4, kfPgenGlobalDosagePresent | kfPgenGlobalMultiallelicDosagePresent, &spgw, &alloc));
  PglErr reterr = SpgwAppendMultiallelicGenovecDphase16(&tv.pgv, 3, &spgw);
  CHECK(!reterr);
  CHECK(!SpgwFinish(&spgw));
  CleanupSpgw(&spgw, &reterr);
  free(alloc);
  FreeTiny(tv);
  WritePvarPsam("tmp_flipfail", allele_idx_offsets, 1, 4);
}

int main(int argc, char** argv) {
  // "compare <orig prefix> <new prefix> [strip-phase] [strip-dosage]
  // [hets-missing] [hct=<thresh>] [fill] [erase=<thresh>]
  // [flip-vars=<variant ID file> flip-samples=<IID file>]" (thresholds in
  // 16384ths of a copy, as plink2 truncates them; '#' lines in the ID files
  // are ignored): see CompareFilesets().
  if ((argc >= 4) && (!strcmp(argv[1], "compare"))) {
    uint32_t strip_phase = 0;
    uint32_t strip_dosage = 0;
    uint32_t hets_missing = 0;
    CallOps ops;
    for (int argi = 4; argi < argc; ++argi) {
      strip_phase |= !strcmp(argv[argi], "strip-phase");
      strip_dosage |= !strcmp(argv[argi], "strip-dosage");
      hets_missing |= !strcmp(argv[argi], "hets-missing");
      ops.fill |= !strcmp(argv[argi], "fill");
      if (!strncmp(argv[argi], "hct=", 4)) {
        ops.hc_thresh = atoi(&(argv[argi][4]));
      } else if (!strncmp(argv[argi], "erase=", 6)) {
        ops.erase_thresh = atoi(&(argv[argi][6]));
      } else if (!strncmp(argv[argi], "flip-vars=", 10) || !strncmp(argv[argi], "flip-samples=", 13)) {
        const uint32_t is_vars = (argv[argi][5] == 'v');
        for (const std::string& line : ReadLines(strchr(argv[argi], '=') + 1)) {
          if ((!line.empty()) && (line[0] != '#')) {
            (is_vars? ops.flip_vars : ops.flip_samples).insert(line);
          }
        }
      }
    }
    return CompareFilesets(argv[2], argv[3], strip_phase, strip_dosage, hets_missing, ops);
  }
  // 300 samples: 2-byte group starts for most allele counts.  70000 samples:
  // 3-byte group starts, and sample IDs past 2^16.
  RunFile("tmp_small.pgen", 300, 6 * 36, 0, 1);
  // combos 24-35: dense rare-alt dosages
  RunFile("tmp_big.pgen", 70000, 6 * 12, 24, 0);

  ValidatorTinyTests();
  CheckFixedWidth();
  WriteFlipFail();

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
    CHECK(SpgwAppendMultiallelicGenovecDphase16(&vd.pgv, 5, &spgw) == kPglRetImproperFunctionCall);
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
    CHECK(SpgwAppendMultiallelicGenovecDphase16(&vd.pgv, 5, &spgw) == kPglRetMalformedInput);
    std::swap(vd.md_codes[offset], vd.md_codes[offset + 1]);
    // Code out of range.
    const AlleleCode saved = vd.md_codes[0];
    vd.md_codes[0] = 5;
    CHECK(SpgwAppendMultiallelicGenovecDphase16(&vd.pgv, 5, &spgw) == kPglRetMalformedInput);
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
      CHECK(SpgwAppendMultiallelicGenovecDphase16(&pgv, 5, &spgw) == kPglRetMalformedInput);
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
      CHECK(SpgwAppendMultiallelicGenovecDphase16(&pgv, 5, &spgw) == kPglRetMalformedInput);
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
