// Benchmark data generator: writes a pair of .pgen/.pvar/.psam filesets with
// dosages, one partly multiallelic and one all-biallelic, from the same
// random genotype model, plus phenotype, covariate, score, and sample-subset
// files.  The all-biallelic fileset is the partly-multiallelic one with every
// non-REF allele collapsed into ALT1, so the two have identical sample and
// variant counts, missingness, phase, and dosage density.
//
// Usage:
//   gen_dosage_data <out prefix> <sample ct> <variant ct> [options]
// Options (defaults in parentheses):
//   multi=<f>    fraction of variants with 3+ alleles (0.2)
//   maxalt=<n>   max ALT alleles of a multiallelic variant, 2..15 (4)
//   dosage=<f>   fraction of (nonmissing) genotypes with fractional dosages
//                (0.3)
//   phase=<f>    fraction of genotypes that are phased (0.5)
//   missing=<f>  fraction of genotypes that are missing (0.01)
//   seed=<n>     random seed (1)
// Writes <out prefix>_multi.{pgen,pvar,psam}, <out prefix>_bi.{...}, and
// <out prefix>.{pheno,covar,score,keep}.
//
// Model, per variant: REF frequency is uniform in [0.5, 0.95]; the remaining
// mass goes mostly to ALT1, with rare ALTs (multiallelic variants only)
// splitting a uniform [0.05, 0.5] share of it.  Per genotype: two haplotype
// alleles are drawn from those frequencies.  A fractional-dosage genotype
// moves a uniform [0, 0.24] share of each haplotype's mass to another allele
// drawn from the same frequencies, which keeps the hardcall consistent with
// the dosages.  Phased genotypes get hardcall phase (if heterozygous) and
// explicit dosage-phase values (if they have dosages).
//
// Build (from this directory):
//   c++ -std=c++17 -O2 -I../../../include -o gen_dosage_data gen_dosage_data.cc
//     ../../../include/pgenlib_write.cc ../../../include/pgenlib_misc.cc
//     ../../../include/plink2_base.cc ../../../include/plink2_bits.cc
//     ../../../include/plink2_simd.cc -lz
// (run_bench.py does this.)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "pgenlib_write.h"

using namespace plink2;

static uint64_t g_rng_state = 1;

static uint64_t NextRand() {
  // splitmix64
  uint64_t z = (g_rng_state += 0x9e3779b97f4a7c15LLU);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9LLU;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebLLU;
  return z ^ (z >> 31);
}

// uniform in [0, 1)
static double RandUnif() {
  return S_CAST(double, NextRand() >> 11) * (1.0 / 9007199254740992.0);
}

static double RandNormal() {
  const double u1 = 1.0 - RandUnif();
  const double u2 = RandUnif();
  return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

static uint32_t DrawAllele(const std::vector<double>& cum_freqs) {
  const double x = RandUnif();
  const uint32_t allele_ct = cum_freqs.size();
  for (uint32_t aidx = 0; aidx + 1 != allele_ct; ++aidx) {
    if (x < cum_freqs[aidx]) {
      return aidx;
    }
  }
  return allele_ct - 1;
}

// Zeroed, cacheline-aligned word buffer with room for vector-width reads.
static uintptr_t* AllocWords(uint32_t bit_ct) {
  const uintptr_t byte_ct = RoundUpPow2(DivUp(bit_ct, 8) + 64, 64);
  void* ptr = aligned_alloc(64, byte_ct);
  memset(ptr, 0, byte_ct);
  return static_cast<uintptr_t*>(ptr);
}

// One genotype in the full (multiallelic) allele space.
struct Geno {
  uint32_t missing;
  uint32_t h1, h2;  // haplotype alleles
  uint32_t phased;
  uint32_t has_dosage;
  // each haplotype's mass: (1 - eps) on h, eps on m (16384 = 1)
  uint32_t m1, m2;
  uint32_t eps1, eps2;
};

// Per-fileset record buffers.
struct RecBufs {
  uintptr_t* genovec;
  uintptr_t* patch_01_set;
  uintptr_t* patch_10_set;
  uintptr_t* phasepresent;
  uintptr_t* phaseinfo;
  uintptr_t* dosage_present;
  uintptr_t* dphase_present;
  uintptr_t* md_present;
  uintptr_t* mdp_present;
  std::vector<AlleleCode> patch_01_vals, patch_10_vals, md_codes, mdp_codes;
  std::vector<uint16_t> dosage_main, md_vals;
  std::vector<int16_t> dphase_delta, mdp_delta;
  std::vector<unsigned char> md_cts, mdp_cts;
};

static void InitRecBufs(uint32_t sample_ct, RecBufs& rb) {
  rb.genovec = AllocWords(2 * sample_ct);
  rb.patch_01_set = AllocWords(sample_ct);
  rb.patch_10_set = AllocWords(sample_ct);
  rb.phasepresent = AllocWords(sample_ct);
  rb.phaseinfo = AllocWords(sample_ct);
  rb.dosage_present = AllocWords(sample_ct);
  rb.dphase_present = AllocWords(sample_ct);
  rb.md_present = AllocWords(sample_ct);
  rb.mdp_present = AllocWords(sample_ct);
}

// Encodes the genotypes, with allele indexes mapped through amap (identity
// for the multiallelic fileset, nonref -> 1 for the biallelic one), and
// appends the record.
static PglErr AppendRecord(const std::vector<Geno>& genos, const std::vector<uint32_t>& amap, uint32_t allele_ct, RecBufs& rb, STPgenWriter* spgwp) {
  const uint32_t sample_ct = genos.size();
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl2, rb.genovec);
  ZeroWArr(sample_ctl, rb.patch_01_set);
  ZeroWArr(sample_ctl, rb.patch_10_set);
  ZeroWArr(sample_ctl, rb.phasepresent);
  ZeroWArr(sample_ctl, rb.phaseinfo);
  ZeroWArr(sample_ctl, rb.dosage_present);
  ZeroWArr(sample_ctl, rb.dphase_present);
  ZeroWArr(sample_ctl, rb.md_present);
  ZeroWArr(sample_ctl, rb.mdp_present);
  rb.patch_01_vals.clear();
  rb.patch_10_vals.clear();
  rb.md_codes.clear();
  rb.mdp_codes.clear();
  rb.dosage_main.clear();
  rb.md_vals.clear();
  rb.dphase_delta.clear();
  rb.mdp_delta.clear();
  rb.md_cts.clear();
  rb.mdp_cts.clear();
  uint32_t phasepresent_ct = 0;
  std::vector<int32_t> u(allele_ct);
  std::vector<int32_t> d(allele_ct);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const Geno& g = genos[sample_idx];
    uintptr_t geno_code;
    if (g.missing) {
      geno_code = 3;
    } else {
      const uint32_t a1 = amap[g.h1];
      const uint32_t a2 = amap[g.h2];
      const uint32_t lo = (a1 < a2)? a1 : a2;
      const uint32_t hi = (a1 < a2)? a2 : a1;
      if (hi == 0) {
        geno_code = 0;
      } else if (lo == 0) {
        geno_code = 1;
        if (hi > 1) {
          SetBit(sample_idx, rb.patch_01_set);
          rb.patch_01_vals.push_back(hi);
        }
      } else {
        geno_code = 2;
        if (hi > 1) {
          SetBit(sample_idx, rb.patch_10_set);
          rb.patch_10_vals.push_back(lo);
          rb.patch_10_vals.push_back(hi);
        }
      }
      if (g.phased && (lo != hi)) {
        SetBit(sample_idx, rb.phasepresent);
        ++phasepresent_ct;
        if (a1 > a2) {
          SetBit(sample_idx, rb.phaseinfo);
        }
      }
      if (g.has_dosage) {
        for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
          u[aidx] = 0;
          d[aidx] = 0;
        }
        const int32_t keep1 = 16384 - g.eps1;
        const int32_t keep2 = 16384 - g.eps2;
        u[a1] += keep1;
        d[a1] += keep1;
        u[amap[g.m1]] += g.eps1;
        d[amap[g.m1]] += g.eps1;
        u[a2] += keep2;
        d[a2] -= keep2;
        u[amap[g.m2]] += g.eps2;
        d[amap[g.m2]] -= g.eps2;
        SetBit(sample_idx, rb.dosage_present);
        rb.dosage_main.push_back(32768 - u[0]);
        unsigned char md_ct = 0;
        for (uint32_t aidx = 2; aidx < allele_ct; ++aidx) {
          if (u[aidx]) {
            rb.md_codes.push_back(aidx);
            rb.md_vals.push_back(u[aidx]);
            ++md_ct;
          }
        }
        if (md_ct) {
          SetBit(sample_idx, rb.md_present);
          rb.md_cts.push_back(md_ct);
        }
        if (g.phased) {
          SetBit(sample_idx, rb.dphase_present);
          rb.dphase_delta.push_back(-d[0]);
          unsigned char mdp_ct = 0;
          for (uint32_t aidx = 2; aidx < allele_ct; ++aidx) {
            if (d[aidx]) {
              rb.mdp_codes.push_back(aidx);
              rb.mdp_delta.push_back(d[aidx]);
              ++mdp_ct;
            }
          }
          if (mdp_ct) {
            SetBit(sample_idx, rb.mdp_present);
            rb.mdp_cts.push_back(mdp_ct);
          }
        }
      }
    }
    rb.genovec[sample_idx / kBitsPerWordD2] |= geno_code << (2 * (sample_idx % kBitsPerWordD2));
  }
  const uint32_t dosage_ct = rb.dosage_main.size();
  const uint32_t dphase_ct = rb.dphase_delta.size();
  if (allele_ct == 2) {
    return SpgwAppendBiallelicGenovecDphase16(rb.genovec, rb.phasepresent, rb.phaseinfo, rb.dosage_present, rb.dphase_present, rb.dosage_main.data(), rb.dphase_delta.data(), dosage_ct, dphase_ct, spgwp);
  }
  PgenVariant pgv;
  memset(&pgv, 0, sizeof(PgenVariant));
  pgv.genovec = rb.genovec;
  pgv.patch_01_set = rb.patch_01_set;
  pgv.patch_01_vals = rb.patch_01_vals.data();
  pgv.patch_10_set = rb.patch_10_set;
  pgv.patch_10_vals = rb.patch_10_vals.data();
  pgv.phasepresent = rb.phasepresent;
  pgv.phaseinfo = rb.phaseinfo;
  pgv.dosage_present = rb.dosage_present;
  pgv.dosage_main = rb.dosage_main.data();
  pgv.dphase_present = rb.dphase_present;
  pgv.dphase_delta = rb.dphase_delta.data();
  pgv.multidosage_present = rb.md_present;
  pgv.multidosage_cts = rb.md_cts.data();
  pgv.multidosage_codes = rb.md_codes.data();
  pgv.multidosage_vals = rb.md_vals.data();
  pgv.multidphase_present = rb.mdp_present;
  pgv.multidphase_cts = rb.mdp_cts.data();
  pgv.multidphase_codes = rb.mdp_codes.data();
  pgv.multidphase_delta = rb.mdp_delta.data();
  pgv.patch_01_ct = rb.patch_01_vals.size();
  pgv.patch_10_ct = rb.patch_10_vals.size() / 2;
  pgv.phasepresent_ct = phasepresent_ct;
  pgv.dosage_ct = dosage_ct;
  pgv.dphase_ct = dphase_ct;
  pgv.multidosage_sample_ct = rb.md_cts.size();
  pgv.multidphase_sample_ct = rb.mdp_cts.size();
  return SpgwAppendMultiallelicGenovecDphase16(&pgv, allele_ct, spgwp);
}

static int OpenWriter(const std::string& fname, const std::vector<uintptr_t>& allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, uint32_t max_allele_ct, PgenGlobalFlags gflags, STPgenWriter* spgwp) {
  PreinitSpgw(spgwp);
  uintptr_t alloc_cacheline_ct;
  uint32_t max_vrec_len;
  if (SpgwInitPhase1(fname.c_str(), allele_idx_offsets.data(), nullptr, variant_ct, sample_ct, max_allele_ct, kPgenWriteBackwardSeek, gflags, 1, spgwp, &alloc_cacheline_ct, &max_vrec_len)) {
    fprintf(stderr, "gen_dosage_data: failed to open %s\n", fname.c_str());
    return 1;
  }
  unsigned char* alloc = static_cast<unsigned char*>(aligned_alloc(kCacheline, alloc_cacheline_ct * kCacheline));
  SpgwInitPhase2(max_vrec_len, spgwp, alloc);
  return 0;
}

static void WritePsam(const std::string& fname, uint32_t sample_ct) {
  FILE* f = fopen(fname.c_str(), "w");
  fputs("#IID\tSEX\n", f);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    fprintf(f, "s%u\t%u\n", sample_idx, 1 + (sample_idx % 2));
  }
  fclose(f);
}

int main(int argc, char** argv) {
  if (argc < 4) {
    fprintf(stderr, "Usage: gen_dosage_data <out prefix> <sample ct> <variant ct> [multi=<f>] [maxalt=<n>] [dosage=<f>] [phase=<f>] [missing=<f>] [seed=<n>]\n");
    return 1;
  }
  const std::string out_prefix = argv[1];
  const uint32_t sample_ct = atoi(argv[2]);
  const uint32_t variant_ct = atoi(argv[3]);
  double multi_frac = 0.2;
  uint32_t max_alt_ct = 4;
  double dosage_frac = 0.3;
  double phase_frac = 0.5;
  double missing_frac = 0.01;
  for (int argi = 4; argi < argc; ++argi) {
    const char* arg = argv[argi];
    if (!strncmp(arg, "multi=", 6)) {
      multi_frac = atof(&arg[6]);
    } else if (!strncmp(arg, "maxalt=", 7)) {
      max_alt_ct = atoi(&arg[7]);
    } else if (!strncmp(arg, "dosage=", 7)) {
      dosage_frac = atof(&arg[7]);
    } else if (!strncmp(arg, "phase=", 6)) {
      phase_frac = atof(&arg[6]);
    } else if (!strncmp(arg, "missing=", 8)) {
      missing_frac = atof(&arg[8]);
    } else if (!strncmp(arg, "seed=", 5)) {
      g_rng_state = strtoull(&arg[5], nullptr, 10);
    } else {
      fprintf(stderr, "gen_dosage_data: unrecognized option '%s'\n", arg);
      return 1;
    }
  }
  if ((!sample_ct) || (!variant_ct) || (max_alt_ct < 2) || (max_alt_ct > 15)) {
    fprintf(stderr, "gen_dosage_data: invalid argument\n");
    return 1;
  }
  // Variant structure first, so that the writers can be initialized.
  std::vector<uint32_t> allele_cts(variant_ct);
  std::vector<uintptr_t> multi_offsets(1, 0);
  std::vector<uintptr_t> bi_offsets(1, 0);
  uint32_t max_allele_ct = 2;
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    uint32_t allele_ct = 2;
    if (RandUnif() < multi_frac) {
      allele_ct = 3 + (NextRand() % (max_alt_ct - 1));
    }
    allele_cts[vidx] = allele_ct;
    if (max_allele_ct < allele_ct) {
      max_allele_ct = allele_ct;
    }
    multi_offsets.push_back(multi_offsets.back() + allele_ct);
    bi_offsets.push_back(bi_offsets.back() + 2);
  }
  const PgenGlobalFlags base_gflags = kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent;
  STPgenWriter multi_spgw;
  STPgenWriter bi_spgw;
  if (OpenWriter(out_prefix + "_multi.pgen", multi_offsets, variant_ct, sample_ct, max_allele_ct, base_gflags | ((max_allele_ct > 2)? kfPgenGlobalMultiallelicDosagePresent : kfPgenGlobal0), &multi_spgw) ||
      OpenWriter(out_prefix + "_bi.pgen", bi_offsets, variant_ct, sample_ct, 2, base_gflags, &bi_spgw)) {
    return 1;
  }
  FILE* multi_pvar = fopen((out_prefix + "_multi.pvar").c_str(), "w");
  FILE* bi_pvar = fopen((out_prefix + "_bi.pvar").c_str(), "w");
  FILE* score = fopen((out_prefix + ".score").c_str(), "w");
  fputs("#CHROM\tPOS\tID\tREF\tALT\n", multi_pvar);
  fputs("#CHROM\tPOS\tID\tREF\tALT\n", bi_pvar);
  fputs("ID\tA1\tWEIGHT\n", score);
  static const char kBases[4] = {'A', 'C', 'G', 'T'};
  RecBufs multi_rb;
  RecBufs bi_rb;
  InitRecBufs(sample_ct, multi_rb);
  InitRecBufs(sample_ct, bi_rb);
  std::vector<Geno> genos(sample_ct);
  const uint32_t variants_per_chr = DivUp(variant_ct, 22);
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    const uint32_t allele_ct = allele_cts[vidx];
    // Allele strings: distinct, REF a single base.
    const uint32_t ref_base = NextRand() % 4;
    std::vector<std::string> alleles(1, std::string(1, kBases[ref_base]));
    for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
      if (aidx < 4) {
        alleles.push_back(std::string(1, kBases[(ref_base + aidx) % 4]));
      } else {
        alleles.push_back(std::string(1, kBases[ref_base]) + kBases[(aidx - 4) % 4] + std::string((aidx - 4) / 4, 'A'));
      }
    }
    // Frequencies.
    std::vector<double> freqs(allele_ct);
    freqs[0] = 0.5 + 0.45 * RandUnif();
    double nonref = 1.0 - freqs[0];
    if (allele_ct == 2) {
      freqs[1] = nonref;
    } else {
      const double rare_share = 0.05 + 0.45 * RandUnif();
      freqs[1] = nonref * (1.0 - rare_share);
      double rare_tot = 0.0;
      std::vector<double> w(allele_ct, 0.0);
      for (uint32_t aidx = 2; aidx != allele_ct; ++aidx) {
        w[aidx] = 0.1 + RandUnif();
        rare_tot += w[aidx];
      }
      for (uint32_t aidx = 2; aidx != allele_ct; ++aidx) {
        freqs[aidx] = nonref * rare_share * w[aidx] / rare_tot;
      }
    }
    std::vector<double> cum_freqs(allele_ct);
    double cum = 0.0;
    for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
      cum += freqs[aidx];
      cum_freqs[aidx] = cum;
    }
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      Geno& g = genos[sample_idx];
      g.missing = (RandUnif() < missing_frac);
      g.h1 = DrawAllele(cum_freqs);
      g.h2 = DrawAllele(cum_freqs);
      g.phased = (RandUnif() < phase_frac);
      g.has_dosage = (!g.missing) && (RandUnif() < dosage_frac);
      g.m1 = g.h1;
      g.m2 = g.h2;
      g.eps1 = 0;
      g.eps2 = 0;
      if (g.has_dosage) {
        g.m1 = DrawAllele(cum_freqs);
        g.m2 = DrawAllele(cum_freqs);
        g.eps1 = (g.m1 == g.h1)? 0 : S_CAST(uint32_t, 0.24 * 16384 * RandUnif());
        g.eps2 = (g.m2 == g.h2)? 0 : S_CAST(uint32_t, 0.24 * 16384 * RandUnif());
      }
    }
    std::vector<uint32_t> identity(allele_ct);
    std::vector<uint32_t> collapse(allele_ct, 1);
    for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
      identity[aidx] = aidx;
    }
    collapse[0] = 0;
    PglErr reterr = AppendRecord(genos, identity, allele_ct, multi_rb, &multi_spgw);
    if (!reterr) {
      reterr = AppendRecord(genos, collapse, 2, bi_rb, &bi_spgw);
    }
    if (reterr) {
      fprintf(stderr, "gen_dosage_data: failed to write variant %u (error code %d)\n", vidx, S_CAST(int, reterr));
      return 1;
    }
    const uint32_t chr = 1 + vidx / variants_per_chr;
    const uint32_t pos = 1000 + 100 * (vidx % variants_per_chr);
    fprintf(multi_pvar, "%u\t%u\tv%u\t%s\t", chr, pos, vidx, alleles[0].c_str());
    for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
      fprintf(multi_pvar, "%s%s", (aidx == 1)? "" : ",", alleles[aidx].c_str());
    }
    fputc('\n', multi_pvar);
    // The collapsed ALT keeps ALT1's string.
    fprintf(bi_pvar, "%u\t%u\tv%u\t%s\t%s\n", chr, pos, vidx, alleles[0].c_str(), alleles[1].c_str());
    fprintf(score, "v%u\t%s\t%.4f\n", vidx, alleles[1].c_str(), RandNormal());
  }
  fclose(multi_pvar);
  fclose(bi_pvar);
  fclose(score);
  if (SpgwFinish(&multi_spgw) || SpgwFinish(&bi_spgw)) {
    fprintf(stderr, "gen_dosage_data: failed to finish .pgen files\n");
    return 1;
  }
  WritePsam(out_prefix + "_multi.psam", sample_ct);
  WritePsam(out_prefix + "_bi.psam", sample_ct);
  FILE* pheno = fopen((out_prefix + ".pheno").c_str(), "w");
  FILE* covar = fopen((out_prefix + ".covar").c_str(), "w");
  FILE* keep = fopen((out_prefix + ".keep").c_str(), "w");
  fputs("#IID\tQT\tCC\n", pheno);
  fputs("#IID\tCOV1\tCOV2\tCOV3\n", covar);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    fprintf(pheno, "s%u\t%.4f\t%u\n", sample_idx, RandNormal(), 1 + (RandUnif() < 0.3));
    fprintf(covar, "s%u\t%.4f\t%.4f\t%.4f\n", sample_idx, RandNormal(), RandNormal(), RandNormal());
    if (sample_idx % 2) {
      fprintf(keep, "s%u\n", sample_idx);
    }
  }
  fclose(pheno);
  fclose(covar);
  fclose(keep);
  return 0;
}
