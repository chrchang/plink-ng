#include "include/plink2_base.h"
#include <algorithm>

#include <Rcpp.h>
using namespace Rcpp;

typedef struct CIndexRecStruct {
  double key;  // yhat on first sort, y on second sort
  uint32_t status;  // original pos on first sort
  uint32_t yhat_int;
  bool operator<(const struct CIndexRecStruct& rhs) const {
    return (key < rhs.key);
  }
} CIndexRec;

uint32_t PopcountLookup(const uintptr_t* new_yhat_present, const uint16_t* new_yhat_64b_popcounts, const uint16_t* new_yhat_2048b_popcounts, const uint32_t* new_yhat_64kib_popcounts, uint32_t widx);

double CIndexTieheavyMain(const CIndexRec* recs, uintptr_t size, unsigned char* wkspace, uintptr_t* new_yhat_present, uint16_t* new_yhat_64b_popcounts, uint16_t* new_yhat_2048b_popcounts, uint32_t* new_yhat_64kib_popcounts);

// [[Rcpp::export]]
double CIndex(NumericVector yhat, NumericVector y, SEXP status) {
  if ((TYPEOF(status) != INTSXP) && (TYPEOF(status) != REALSXP)) {
    // could support LGLSXP too
    stop("Unsupported status type");
  }
  const uintptr_t size = yhat.size();
  if (size == 0) {
    stop("yhat cannot be empty");
  }
  if (size != static_cast<uintptr_t>(y.size())) {
    stop("yhat and y must have the same length");
  }

  // We define C-index as:
  //   sum_{all k where status[k]=1} [#(m: y[m]>y[j], yhat[m]<yhat[j]) +
  //                                  0.5 * #(m: y[m]>y[j], yhat[m]==yhat[j])]
  //   ------------------------------------------------------------------------
  //   sum_{all k where status[k]=1} #(m: y[m]>y[j])
  //
  // The computation is organized as follows:
  // 1. Replace yhat with 0-based ranks ('yhat_int'); note that this has no
  //    effect on the expression above.  Ties round down.
  // 2. Sort records in order of nondecreasing y.
  // 3. Iterate through records, starting from largest y.  Incrementally update
  //    a data structure tracking which yhat_int values have been seen so far,
  //    and in case of ties, their multiplicity.  This data structure supports
  //    fast insertion and number-of-smaller{-or-equal}-elements queries, which
  //    covers the operations we need to compute the numerator.  Formally, the
  //    total query time is still O(N^2) if status==1 isn't sparse, but the
  //    constant multiplier is so small that the std::sort calls practically
  //    always dominate for R's N<2^31 domain.  (Not literally true for current
  //    implementation, but it's trivial to branch on N and use another
  //    popcount index-levels for very large N once any application cares.)
  // More precisely, the new_yhat data structure consists of a bitarray, and
  // three index-arrays containing partial popcounts for intervals of 2^9,
  // 2^14, and 2^19 bits.

  // Intermediate data structures (all cacheline-aligned except maybe recs,
  // which is vector-aligned):
  // - recs: N CIndexRec
  // - new_yhat_present: N bits
  // - new_yhat_64b_popcounts: DivUp(N, 512) uint16s
  // - new_yhat_2048b_popcounts: DivUp(N, 512*32) uint16s
  // - new_yhat_64kib_popcounts: DivUp(N, 512*32*32) uint16s
  const uintptr_t recs_vec_ct = plink2::DivUp(sizeof(CIndexRec) * size, plink2::kBytesPerVec);
  const uintptr_t new_yhat_present_cacheline_ct = plink2::DivUp(size, plink2::kBitsPerCacheline);
  const uintptr_t new_yhat_present_vec_ct = new_yhat_present_cacheline_ct * plink2::kVecsPerCacheline;
  const uintptr_t new_yhat_64b_popcounts_cacheline_ct = plink2::DivUp(new_yhat_present_cacheline_ct, plink2::kInt16PerCacheline);
  const uintptr_t new_yhat_64b_popcounts_vec_ct = new_yhat_64b_popcounts_cacheline_ct * plink2::kVecsPerCacheline;
  const uintptr_t new_yhat_2048b_popcounts_cacheline_ct = plink2::DivUp(new_yhat_64b_popcounts_cacheline_ct, plink2::kInt16PerCacheline);
  const uintptr_t new_yhat_2048b_popcounts_vec_ct = new_yhat_2048b_popcounts_cacheline_ct * plink2::kVecsPerCacheline;
  const uintptr_t new_yhat_64kib_popcounts_vec_ct = plink2::DivUp(new_yhat_2048b_popcounts_cacheline_ct, plink2::kInt32PerVec);
  const uintptr_t vec_ct = recs_vec_ct + new_yhat_present_vec_ct + new_yhat_64b_popcounts_vec_ct + new_yhat_2048b_popcounts_vec_ct + new_yhat_64kib_popcounts_vec_ct;
  unsigned char* wkspace;
  if (plink2::cachealigned_malloc(vec_ct * plink2::kBytesPerVec, &wkspace)) {
    stop("Out of memory");
  }
  memset(wkspace, 0, (vec_ct - recs_vec_ct) * plink2::kBytesPerVec);
  unsigned char* wkspace_iter = wkspace;
  // Put new_yhat arrays first to guarantee cacheline alignment, since it
  // matters there.
  uintptr_t* new_yhat_present = reinterpret_cast<uintptr_t*>(wkspace_iter);
  wkspace_iter = &(wkspace_iter[new_yhat_present_vec_ct * plink2::kBytesPerVec]);
  uint16_t* new_yhat_64b_popcounts = reinterpret_cast<uint16_t*>(wkspace_iter);
  wkspace_iter = &(wkspace_iter[new_yhat_64b_popcounts_vec_ct * plink2::kBytesPerVec]);
  uint16_t* new_yhat_2048b_popcounts = reinterpret_cast<uint16_t*>(wkspace_iter);
  wkspace_iter = &(wkspace_iter[new_yhat_2048b_popcounts_vec_ct * plink2::kBytesPerVec]);
  uint32_t* new_yhat_64kib_popcounts = reinterpret_cast<uint32_t*>(wkspace_iter);
  wkspace_iter = &(wkspace_iter[new_yhat_64kib_popcounts_vec_ct * plink2::kBytesPerVec]);
  CIndexRec* recs = reinterpret_cast<CIndexRec*>(wkspace_iter);

  uint64_t yhat_int_sum = 0;
  {
    const double* yhatd = &(yhat[0]);
    for (uintptr_t ulii = 0; ulii != size; ++ulii) {
      recs[ulii].key = yhatd[ulii];
      recs[ulii].status = ulii;
    }
    std::sort(recs, &(recs[size]));
    recs[recs[0].status].yhat_int = 0;
    double prev_yhat = recs[0].key;
    uint32_t prev_idx = 0;
    for (uintptr_t ulii = 1; ulii != size; ++ulii) {
      const double cur_yhat = recs[ulii].key;
      if (cur_yhat != prev_yhat) {
        prev_idx = ulii;
        prev_yhat = cur_yhat;
      }
      recs[recs[ulii].status].yhat_int = prev_idx;
      yhat_int_sum += prev_idx;
    }
  }
  const double* yd = &(y[0]);
  if (TYPEOF(status) == INTSXP) {
    IntegerVector status_iv = as<IntegerVector>(status);
    if (size != static_cast<uintptr_t>(status_iv.size())) {
      plink2::aligned_free(wkspace);
      stop("y and status must have the same length");
    }
    const int32_t* status_ivi = &(status_iv[0]);
    for (uintptr_t ulii = 0; ulii != size; ++ulii) {
      recs[ulii].key = yd[ulii];
      const uint32_t uii = status_ivi[ulii];
      if (uii & (~1)) {
        plink2::aligned_free(wkspace);
        stop("all status values must be 0 or 1");
      }
      recs[ulii].status = uii;
    }
  } else {
    NumericVector status_nv = as<NumericVector>(status);
    if (size != static_cast<uintptr_t>(status_nv.size())) {
      plink2::aligned_free(wkspace);
      stop("y and status must have the same length");
    }
    const double* status_nvd = &(status_nv[0]);
    for (uintptr_t ulii = 0; ulii != size; ++ulii) {
      recs[ulii].key = yd[ulii];
      const double dxx = status_nvd[ulii];
      if (dxx == 0.0) {
        recs[ulii].status = 0;
      } else if (dxx == 1.0) {
        recs[ulii].status = 1;
      } else {
        plink2::aligned_free(wkspace);
        stop("all status values must be 0 or 1");
      }
    }
  }
  std::sort(recs, &(recs[size]));

  // With no ties, yhat_int_sum would be (size * (size - 1)) / 2.
  // k-way ties decrease the sum by (k * (k-1)) / 2, which is conveniently
  // proportional to the burden imposed by the ties (sequential scanning for
  // the next unset bit in new_yhat_present) on our usual algorithm.
  // With high enough tie-burden (threshold proportional to n log n), we switch
  // to a slightly different approach, which is usually slower but avoids
  // sequential scans in case of ties.
  const uint64_t tie_burden = (static_cast<uint64_t>(size) * (size - 1)) / 2 - yhat_int_sum;
  if (tie_burden > 32 * static_cast<uint64_t>(size) * plink2::bsrw(size)) {
    return CIndexTieheavyMain(recs, size, wkspace, new_yhat_present, new_yhat_64b_popcounts, new_yhat_2048b_popcounts, new_yhat_64kib_popcounts);
  }

  // Main loop.
  double prev_key = recs[size - 1].key;
  uint32_t recs_block_end = size;
  uint32_t new_yhat_entry_ct = 0;
  uint64_t numer_x2 = 0;
  uint64_t denom = 0;
  for (uint32_t recs_idx = size; recs_idx; ) {
    --recs_idx;
    const double cur_key = recs[recs_idx].key;
    if (cur_key != prev_key) {
      const uint32_t next_block_end = recs_idx + 1;
      for (uint32_t uii = next_block_end; uii != recs_block_end; ++uii) {
        const uint32_t yhat_int = recs[uii].yhat_int;
        // This is essentially an inlined AdvTo0Bit() call.
        uintptr_t widx = yhat_int / plink2::kBitsPerWord;
        {
          const uintptr_t remainder_bit = plink2::k1LU << (yhat_int % plink2::kBitsPerWord);
          uintptr_t cur_word = new_yhat_present[widx];
          if (!(remainder_bit & cur_word)) {
            // common case
            new_yhat_present[widx] = cur_word | remainder_bit;
          } else {
            // There's a tie, and we've already seen at least one previous
            // instance.
            // Find the rightmost clear bit in cur_word after remainder_bit, if
            // there is one.  We use the bitarray property that cur_word +
            // <bit> sets the first 0-bit in cur_word at or past the added bit,
            // may clear some other bits, and doesn't do anything else.
            uintptr_t ulii = cur_word + remainder_bit;
            if (ulii < remainder_bit) {
              // Addition 'carry' needed.
              do {
                cur_word = new_yhat_present[++widx];
              } while (cur_word == ~plink2::k0LU);
              ulii = cur_word + 1;
            }
            const uintptr_t new_word = ulii | cur_word;
            new_yhat_present[widx] = new_word;
          }
        }
        new_yhat_64b_popcounts[widx / plink2::kWordsPerCacheline] += 1;
        new_yhat_2048b_popcounts[widx / (plink2::kWordsPerCacheline * plink2::kInt16PerCacheline)] += 1;
        new_yhat_64kib_popcounts[widx / (plink2::kWordsPerCacheline * plink2::kInt16PerCacheline * plink2::kInt16PerCacheline)] += 1;
      }
      recs_block_end = next_block_end;
      new_yhat_entry_ct = size - recs_block_end;
      prev_key = cur_key;
    }
    if (!recs[recs_idx].status) {
      continue;
    }
    denom += new_yhat_entry_ct;
    const uint32_t yhat_int = recs[recs_idx].yhat_int;
    // numer_x2 += 2 * <# of set bits before yhat_int> + <# of ties>
    uint32_t widx = yhat_int / plink2::kBitsPerWord;
    uint32_t n_set_bits_before = PopcountLookup(new_yhat_present, new_yhat_64b_popcounts, new_yhat_2048b_popcounts, new_yhat_64kib_popcounts, widx);
    uint32_t remainder = yhat_int % plink2::kBitsPerWord;
    uintptr_t cur_word = new_yhat_present[widx];
    uintptr_t remainder_bit = plink2::k1LU << remainder;
    n_set_bits_before += plink2::PopcountWord(cur_word & (remainder_bit - 1));
    numer_x2 += 2 * n_set_bits_before;

    if (!(cur_word & remainder_bit)) {
      // common case
      continue;
    }
    // Similar to the insertion tie-handling code above.
    numer_x2 -= remainder;
    uintptr_t ulii = cur_word + remainder_bit;
    if (ulii < remainder_bit) {
      uint32_t widx2 = widx;
      do {
        cur_word = new_yhat_present[++widx2];
      } while (cur_word == ~plink2::k0LU);
      numer_x2 += (widx2 - widx) * plink2::kBitsPerWord;
      ulii = cur_word + 1;
    }
    numer_x2 += plink2::ctzw(ulii & (~cur_word));
  }
  plink2::aligned_free(wkspace);
  return plink2::u63tod(numer_x2) / plink2::u63tod(denom * 2);
}

uint32_t PopcountLookup(const uintptr_t* new_yhat_present, const uint16_t* new_yhat_64b_popcounts, const uint16_t* new_yhat_2048b_popcounts, const uint32_t* new_yhat_64kib_popcounts, uint32_t widx) {
  uint32_t result = 0;
  const uint32_t idx_64kib = widx / (plink2::kWordsPerCacheline * plink2::kInt16PerCacheline * plink2::kInt16PerCacheline);
  for (uint32_t uii = 0; uii != idx_64kib; ++uii) {
    result += new_yhat_64kib_popcounts[uii];
  }
  const uint32_t idx_2048b = widx / (plink2::kWordsPerCacheline * plink2::kInt16PerCacheline);
  for (uint32_t uii = idx_64kib * plink2::kInt16PerCacheline; uii != idx_2048b; ++uii) {
    result += new_yhat_2048b_popcounts[uii];
  }

  // Each entry in the bottom-level-index array is <= 512, and we're adding
  // at most 16 of them even when compiling as 32-bit, so there's no overflow
  // danger if we add entire words at a time (effectively performing 2-4
  // uint16 additions in parallel).
  const uintptr_t* walias_64b_popcounts = reinterpret_cast<const uintptr_t*>(new_yhat_64b_popcounts);
  const uint32_t widx_64b = widx / (plink2::kWordsPerCacheline * plink2::kInt16PerWord);
  uintptr_t acc = 0;
  for (uint32_t uii = idx_2048b * plink2::kWordsPerCacheline; uii != widx_64b; ++uii) {
    acc += walias_64b_popcounts[uii];
  }
  const uint32_t idx_64b = widx / plink2::kWordsPerCacheline;
  uint32_t remainder = idx_64b % plink2::kInt16PerWord;
  if (remainder) {
    acc += plink2::bzhi(walias_64b_popcounts[widx_64b], remainder * 16);
  }
  result += (acc * plink2::kMask0001) >> (plink2::kBitsPerWord - 16);

  for (uint32_t uii = idx_64b * plink2::kWordsPerCacheline; uii != widx; ++uii) {
    result += plink2::PopcountWord(new_yhat_present[uii]);
  }
  return result;
}

double CIndexTieheavyMain(const CIndexRec* recs, uintptr_t size, unsigned char* wkspace, uintptr_t* new_yhat_present, uint16_t* new_yhat_64b_popcounts, uint16_t* new_yhat_2048b_popcounts, uint32_t* new_yhat_64kib_popcounts) {
  // Alternate tie-handling strategy: maintain an additional array which tracks
  // the number of times we've seen each yhat_int value.
  uint32_t* yhat_int_multiplicities = static_cast<uint32_t*>(calloc(size, sizeof(int32_t)));
  if (!yhat_int_multiplicities) {
    plink2::aligned_free(wkspace);
    stop("Out of memory");
  }

  double prev_key = recs[size - 1].key;
  uint32_t recs_block_end = size;
  uint32_t new_yhat_entry_ct = 0;
  uint64_t numer_x2 = 0;
  uint64_t denom = 0;
  for (uint32_t recs_idx = size; recs_idx; ) {
    --recs_idx;
    const double cur_key = recs[recs_idx].key;
    if (cur_key != prev_key) {
      const uint32_t next_block_end = recs_idx + 1;
      for (uint32_t uii = next_block_end; uii != recs_block_end; ++uii) {
        const uint32_t yhat_int = recs[uii].yhat_int;
        const uint32_t multiplicity = yhat_int_multiplicities[yhat_int];
        const uint32_t adj_yhat_int = yhat_int + multiplicity;
        yhat_int_multiplicities[yhat_int] = multiplicity + 1;
        const uintptr_t widx = adj_yhat_int / plink2::kBitsPerWord;
        const uintptr_t remainder_bit = plink2::k1LU << (adj_yhat_int % plink2::kBitsPerWord);
        new_yhat_present[widx] |= remainder_bit;
        new_yhat_64b_popcounts[widx / plink2::kWordsPerCacheline] += 1;
        new_yhat_2048b_popcounts[widx / (plink2::kWordsPerCacheline * plink2::kInt16PerCacheline)] += 1;
        new_yhat_64kib_popcounts[widx / (plink2::kWordsPerCacheline * plink2::kInt16PerCacheline * plink2::kInt16PerCacheline)] += 1;
      }
      recs_block_end = next_block_end;
      new_yhat_entry_ct = size - recs_block_end;
      prev_key = cur_key;
    }
    if (!recs[recs_idx].status) {
      continue;
    }
    denom += new_yhat_entry_ct;
    const uint32_t yhat_int = recs[recs_idx].yhat_int;
    uint32_t widx = yhat_int / plink2::kBitsPerWord;
    uint32_t n_set_bits_before = PopcountLookup(new_yhat_present, new_yhat_64b_popcounts, new_yhat_2048b_popcounts, new_yhat_64kib_popcounts, widx);
    uint32_t remainder = yhat_int % plink2::kBitsPerWord;
    uintptr_t cur_word = new_yhat_present[widx];
    uintptr_t remainder_bit = plink2::k1LU << remainder;
    n_set_bits_before += plink2::PopcountWord(cur_word & (remainder_bit - 1));
    numer_x2 += 2 * n_set_bits_before + yhat_int_multiplicities[yhat_int];
  }
  free(yhat_int_multiplicities);
  plink2::aligned_free(wkspace);
  return plink2::u63tod(numer_x2) / plink2::u63tod(denom * 2);
}
