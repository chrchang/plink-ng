#ifndef __PLINK_PERM_H__

// Permutation generation and interpretation code common to many association
// tests.

void generate_cc_perm_vec(uint32_t tot_ct, uint32_t set_ct, uint32_t tot_quotient, uint64_t totq_magic, uint32_t totq_preshift, uint32_t totq_postshift, uint32_t totq_incr, uintptr_t* perm_vec, sfmt_t* sfmtp);

void generate_cc_perm1(uint32_t tot_ct, uint32_t set_ct, uint32_t tot_quotient, uint64_t totq_magic, uint32_t totq_preshift, uint32_t totq_postshift, uint32_t totq_incr, uintptr_t* perm_vec, sfmt_t* sfmtp);

void generate_cc_cluster_perm_vec(uint32_t tot_ct, uintptr_t* preimage, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts, uint32_t* tot_quotients, uint64_t* totq_magics, uint32_t* totq_preshifts, uint32_t* totq_postshifts, uint32_t* totq_incrs, uintptr_t* perm_vec, sfmt_t* sfmtp);

void generate_cc_cluster_perm1(uint32_t tot_ct, uintptr_t* preimage, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* cluster_case_cts, uint32_t* tot_quotients, uint64_t* totq_magics, uint32_t* totq_preshifts, uint32_t* totq_postshifts, uint32_t* totq_incrs, uintptr_t* perm_vec, sfmt_t* sfmtp);

// Efficient "vertical popcount" support.
void transpose_perms(uintptr_t* perm_vecs, uint32_t perm_vec_ct, uint32_t pheno_nm_ct, uint32_t* perm_vecst);

void transpose_perm1s(uintptr_t* perm_vecs, uint32_t perm_vec_ct, uint32_t pheno_nm_ct, uint32_t* perm_vecst);

#ifdef __LP64__
static inline void unroll_incr_1_4(const __m128i* acc1, __m128i* acc4, uint32_t acc1_vec_ct) {
  const __m128i m1x4 = {0x1111111111111111LLU, 0x1111111111111111LLU};
  __m128i loader;
  uint32_t vidx;
  for (vidx = 0; vidx < acc1_vec_ct; vidx++) {
    loader = *acc1++;
    *acc4 = _mm_add_epi64(*acc4, _mm_and_si128(loader, m1x4));
    acc4++;
    loader = _mm_srli_epi64(loader, 1);
    *acc4 = _mm_add_epi64(*acc4, _mm_and_si128(loader, m1x4));
    acc4++;
    loader = _mm_srli_epi64(loader, 1);
    *acc4 = _mm_add_epi64(*acc4, _mm_and_si128(loader, m1x4));
    acc4++;
    loader = _mm_srli_epi64(loader, 1);
    *acc4 = _mm_add_epi64(*acc4, _mm_and_si128(loader, m1x4));
    acc4++;
  }
}

static inline void unroll_incr_4_8(const __m128i* acc4, __m128i* acc8, uint32_t acc4_vec_ct) {
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  uint32_t vidx;
  for (vidx = 0; vidx < acc4_vec_ct; vidx++) {
    loader = *acc4++;
    *acc8 = _mm_add_epi64(*acc8, _mm_and_si128(loader, m4));
    acc8++;
    loader = _mm_srli_epi64(loader, 4);
    *acc8 = _mm_add_epi64(*acc8, _mm_and_si128(loader, m4));
    acc8++;
  }
}

static inline void unroll_zero_incr_4_8(__m128i* acc4, __m128i* acc8, uint32_t acc4_vec_ct) {
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i loader;
  uint32_t vidx;
  for (vidx = 0; vidx < acc4_vec_ct; vidx++) {
    loader = *acc4;
    *acc4++ = _mm_setzero_si128();
    *acc8 = _mm_add_epi64(*acc8, _mm_and_si128(loader, m4));
    acc8++;
    loader = _mm_srli_epi64(loader, 4);
    *acc8 = _mm_add_epi64(*acc8, _mm_and_si128(loader, m4));
    acc8++;
  }
}

static inline void unroll_incr_8_32(const __m128i* acc8, __m128i* acc32, uint32_t acc8_vec_ct) {
  const __m128i m8x32 = {0x000000ff000000ffLLU, 0x000000ff000000ffLLU};
  __m128i loader;
  uint32_t vidx;
  for (vidx = 0; vidx < acc8_vec_ct; vidx++) {
    loader = *acc8++;
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
    loader = _mm_srli_epi64(loader, 8);
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
    loader = _mm_srli_epi64(loader, 8);
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
    loader = _mm_srli_epi64(loader, 8);
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
  }
}

static inline void unroll_zero_incr_8_32(__m128i* acc8, __m128i* acc32, uint32_t acc8_vec_ct) {
  const __m128i m8x32 = {0x000000ff000000ffLLU, 0x000000ff000000ffLLU};
  __m128i loader;
  uint32_t vidx;
  for (vidx = 0; vidx < acc8_vec_ct; vidx++) {
    loader = *acc8;
    *acc8++ = _mm_setzero_si128();
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
    loader = _mm_srli_epi64(loader, 8);
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
    loader = _mm_srli_epi64(loader, 8);
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
    loader = _mm_srli_epi64(loader, 8);
    *acc32 = _mm_add_epi64(*acc32, _mm_and_si128(loader, m8x32));
    acc32++;
  }
}
#else
static inline void unroll_incr_1_4(const uintptr_t* acc1, uintptr_t* acc4, uint32_t acc1_word_ct) {
  uint32_t widx;
  uintptr_t loader;
  for (widx = 0; widx < acc1_word_ct; widx++) {
    loader = *acc1++;
    *acc4 += loader & 0x11111111U;
    acc4++;
    loader >>= 1;
    *acc4 += loader & 0x11111111U;
    acc4++;
    loader >>= 1;
    *acc4 += loader & 0x11111111U;
    acc4++;
    loader >>= 1;
    *acc4 += loader & 0x11111111U;
    acc4++;
  }
}

static inline void unroll_incr_4_8(const uintptr_t* acc4, uintptr_t* acc8, uint32_t acc4_word_ct) {
  uint32_t widx;
  uintptr_t loader;
  for (widx = 0; widx < acc4_word_ct; widx++) {
    loader = *acc4++;
    *acc8 += loader & 0x0f0f0f0fU;
    acc8++;
    loader >>= 4;
    *acc8 += loader & 0x0f0f0f0fU;
    acc8++;
  }
}

static inline void unroll_zero_incr_4_8(uintptr_t* acc4, uintptr_t* acc8, uint32_t acc4_word_ct) {
  uint32_t widx;
  uintptr_t loader;
  for (widx = 0; widx < acc4_word_ct; widx++) {
    loader = *acc4;
    *acc4++ = 0;
    *acc8 += loader & 0x0f0f0f0fU;
    acc8++;
    loader >>= 4;
    *acc8 += loader & 0x0f0f0f0fU;
    acc8++;
  }
}

static inline void unroll_incr_8_32(const uintptr_t* acc8, uintptr_t* acc32, uint32_t acc8_word_ct) {
  uint32_t widx;
  uintptr_t loader;
  for (widx = 0; widx < acc8_word_ct; widx++) {
    loader = *acc8++;
    *acc32 += (uint8_t)loader;
    acc32++;
    loader >>= 8;
    *acc32 += (uint8_t)loader;
    acc32++;
    loader >>= 8;
    *acc32 += (uint8_t)loader;
    acc32++;
    loader >>= 8;
    *acc32 += loader;
    acc32++;
  }
}

static inline void unroll_zero_incr_8_32(uintptr_t* acc8, uintptr_t* acc32, uint32_t acc8_word_ct) {
  uint32_t widx;
  uintptr_t loader;
  for (widx = 0; widx < acc8_word_ct; widx++) {
    loader = *acc8;
    *acc8++ = 0;
    *acc32 += (uint8_t)loader;
    acc32++;
    loader >>= 8;
    *acc32 += (uint8_t)loader;
    acc32++;
    loader >>= 8;
    *acc32 += (uint8_t)loader;
    acc32++;
    loader >>= 8;
    *acc32 += loader;
    acc32++;
  }
}
#endif

#endif // __PLINK_PERM_H__
