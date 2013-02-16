#include "wdist_common.h"
#include "wdist_stats.h"

// load markers in blocks to enable PERMORY-style LD exploitation
#define MODEL_BLOCKSIZE 256
#define MODEL_HALFBLOCK 128

/*
#ifdef __LP64__
void freq_cc_count_120v(__m128i* vptr, __m128i* vend, __m128i* maskvp, uint32_t* ctap, uint32_t* ctbp, uint32_t* ctcp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  const __m128i m16 = {0x0000ffff0000ffffLLU, 0x0000ffff0000ffffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i to_ct_a1;
  __m128i to_ct_b1;
  __m128i to_ct_c1;
  __m128i to_ct_a2;
  __m128i to_ct_b2;
  __m128i to_ct_c2;
  __uni16 acc_a;
  __uni16 acc_b;
  __uni16 acc_c;

  acc_a.vi = _mm_setzero_si128();
  acc_b.vi = _mm_setzero_si128();
  acc_c.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader2 = *maskvp++;
    to_ct_b1 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    to_ct_a1 = _mm_and_si128(loader2, loader);
    to_ct_c1 = _mm_and_si128(to_ct_b1, loader);
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    to_ct_a1 = _mm_add_epi64(to_ct_a1, _mm_and_si128(loader2, loader));
    to_ct_b1 = _mm_add_epi64(to_ct_b1, loader3);
    to_ct_c1 = _mm_add_epi64(to_ct_c1, _mm_and_si128(loader3, loader));
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    to_ct_a1 = _mm_add_epi64(to_ct_a1, _mm_and_si128(loader2, loader));
    to_ct_b1 = _mm_add_epi64(to_ct_b1, loader3);
    to_ct_c1 = _mm_add_epi64(to_ct_c1, _mm_and_si128(loader3, loader));

    to_ct_a1 = _mm_add_epi64(_mm_and_si128(to_ct_a1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_a1, 2), m2));
    to_ct_b1 = _mm_add_epi64(_mm_and_si128(to_ct_b1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_b1, 2), m2));
    to_ct_c1 = _mm_add_epi64(_mm_and_si128(to_ct_c1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_c1, 2), m2));

    loader = *vptr++;
    loader2 = *maskvp++;
    to_ct_b2 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    to_ct_a2 = _mm_and_si128(loader2, loader);
    to_ct_c2 = _mm_and_si128(to_ct_b2, loader);
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    to_ct_a2 = _mm_add_epi64(to_ct_a2, _mm_and_si128(loader2, loader));
    to_ct_b2 = _mm_add_epi64(to_ct_b2, loader3);
    to_ct_c2 = _mm_add_epi64(to_ct_c2, _mm_and_si128(loader3, loader));
    loader = *vptr++;
    loader2 = *maskvp++;
    loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
    to_ct_a2 = _mm_add_epi64(to_ct_a2, _mm_and_si128(loader2, loader));
    to_ct_b2 = _mm_add_epi64(to_ct_b2, loader3);
    to_ct_c2 = _mm_add_epi64(to_ct_c2, _mm_and_si128(loader3, loader));

    to_ct_a1 = _mm_add_epi64(to_ct_a1, _mm_add_epi64(_mm_and_si128(to_ct_a2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_a2, 2), m2)));
    to_ct_b1 = _mm_add_epi64(to_ct_b1, _mm_add_epi64(_mm_and_si128(to_ct_b2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_b2, 2), m2)));
    to_ct_c1 = _mm_add_epi64(to_ct_c1, _mm_add_epi64(_mm_and_si128(to_ct_c2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_c2, 2), m2)));

    acc_a.vi = _mm_add_epi64(acc_a.vi, _mm_add_epi64(_mm_and_si128(to_ct_a1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_a1, 4), m4)));
    acc_b.vi = _mm_add_epi64(acc_b.vi, _mm_add_epi64(_mm_and_si128(to_ct_b1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_b1, 4), m4)));
    acc_c.vi = _mm_add_epi64(acc_c.vi, _mm_add_epi64(_mm_and_si128(to_ct_c1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_c1, 4), m4)));
  } while (vptr < vend);
  acc_a.vi = _mm_add_epi64(_mm_and_si128(acc_a.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_a.vi, 8), m8));
  acc_b.vi = _mm_add_epi64(_mm_and_si128(acc_b.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_b.vi, 8), m8));
  acc_c.vi = _mm_add_epi64(_mm_and_si128(acc_c.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_c.vi, 8), m8));
  acc_a.vi = _mm_and_si128(_mm_add_epi64(acc_a.vi, _mm_srli_epi64(acc_a.vi, 16)), m16);
  acc_b.vi = _mm_and_si128(_mm_add_epi64(acc_b.vi, _mm_srli_epi64(acc_b.vi, 16)), m16);
  acc_c.vi = _mm_and_si128(_mm_add_epi64(acc_c.vi, _mm_srli_epi64(acc_c.vi, 16)), m16);
  acc_a.vi = _mm_add_epi64(acc_a.vi, _mm_srli_epi64(acc_a.vi, 32));
  acc_b.vi = _mm_add_epi64(acc_b.vi, _mm_srli_epi64(acc_b.vi, 32));
  acc_c.vi = _mm_add_epi64(acc_c.vi, _mm_srli_epi64(acc_c.vi, 32));
  *ctap += (uint32_t)(acc_a.u8[0] + acc_a.u8[1]);
  *ctbp += (uint32_t)(acc_b.u8[0] + acc_b.u8[1]);
  *ctcp += (uint32_t)(acc_c.u8[0] + acc_c.u8[1]);
}

void freq_hwe_haploid_count_120v(__m128i* vptr, __m128i* vend, __m128i* maskvp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
  const __m128i m16 = {0x0000ffff0000ffffLLU, 0x0000ffff0000ffffLLU};
  __m128i loader;
  __m128i loader2;
  __m128i loader3;
  __m128i to_ct_nm1;
  __m128i to_ct_hmaj1;
  __m128i to_ct_nm2;
  __m128i to_ct_hmaj2;
  __uni16 acc_nm;
  __uni16 acc_hmaj;

  acc_nm.vi = _mm_setzero_si128();
  acc_hmaj.vi = _mm_setzero_si128();
  do {
    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj1 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(_mm_and_si128(to_ct_nm1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 2), m2));
    to_ct_hmaj1 = _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 2), m2));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_andnot_si128(loader2, loader3);
    to_ct_hmaj2 = _mm_and_si128(loader, loader3);

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    loader = *vptr++;
    loader3 = _mm_srli_epi64(loader, 1);
    loader2 = _mm_xor_si128(loader, loader3); // inverted
    loader = _mm_and_si128(loader, loader3);
    loader3 = *maskvp++;
    to_ct_nm2 = _mm_add_epi64(to_ct_nm2, _mm_andnot_si128(loader2, loader3));
    to_ct_hmaj2 = _mm_add_epi64(to_ct_hmaj2, _mm_and_si128(loader, loader3));

    to_ct_nm1 = _mm_add_epi64(to_ct_nm1, _mm_add_epi64(_mm_and_si128(to_ct_nm2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_nm2, 2), m2)));
    to_ct_hmaj1 = _mm_add_epi64(to_ct_hmaj1, _mm_add_epi64(_mm_and_si128(to_ct_hmaj2, m2), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj2, 2), m2)));

    acc_nm.vi = _mm_add_epi64(acc_nm.vi, _mm_add_epi64(_mm_and_si128(to_ct_nm1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_nm1, 4), m4)));
    acc_hmaj.vi = _mm_add_epi64(acc_hmaj.vi, _mm_add_epi64(_mm_and_si128(to_ct_hmaj1, m4), _mm_and_si128(_mm_srli_epi64(to_ct_hmaj1, 4), m4)));
  } while (vptr < vend);
  acc_nm.vi = _mm_add_epi64(_mm_and_si128(acc_nm.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_nm.vi, 8), m8));
  acc_hmaj.vi = _mm_add_epi64(_mm_and_si128(acc_hmaj.vi, m8), _mm_and_si128(_mm_srli_epi64(acc_hmaj.vi, 8), m8));
  acc_nm.vi = _mm_and_si128(_mm_add_epi64(acc_nm.vi, _mm_srli_epi64(acc_nm.vi, 16)), m16);
  acc_hmaj.vi = _mm_and_si128(_mm_add_epi64(acc_hmaj.vi, _mm_srli_epi64(acc_hmaj.vi, 16)), m16);
  acc_nm.vi = _mm_add_epi64(acc_nm.vi, _mm_srli_epi64(acc_nm.vi, 32));
  acc_hmaj.vi = _mm_add_epi64(acc_hmaj.vi, _mm_srli_epi64(acc_hmaj.vi, 32));
  *ct_nmp += (uint32_t)(acc_nm.u8[0] + acc_nm.u8[1]);
  *ct_hmajp += (uint32_t)(acc_hmaj.u8[0] + acc_hmaj.u8[1]);
}
#else
void freq_hwe_count_12(uintptr_t* lptr, uintptr_t* maskp, uint32_t* ctap, uint32_t* ctbp, uint32_t* ctcp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader2 = *maskp++;
  uint32_t to_ct_a1 = loader & loader2;
  uint32_t to_ct_b1 = (loader >> 1) & loader2;
  uint32_t to_ct_c1 = loader & to_ct_b1;
  uintptr_t loader3;
  uint32_t to_ct_a2;
  uint32_t to_ct_b2;
  uint32_t to_ct_c2;
  uintptr_t partial_a;
  uintptr_t partial_b;
  uintptr_t partial_c;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;

  loader = *lptr++;
  loader2 = *maskp++;
  to_ct_a2 = loader & loader2;
  to_ct_b2 = (loader >> 1) & loader2;
  to_ct_c2 = loader & to_ct_b2;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;

  to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
  to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
  partial_a = (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
  to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
  to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
  partial_b = (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
  to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
  to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
  partial_c = (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader2 = *maskp++;
  to_ct_a1 = loader & loader2;
  to_ct_b1 = (loader >> 1) & loader2;
  to_ct_c1 = loader & to_ct_b1;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a1 += loader & loader2;
  to_ct_b1 += loader3;
  to_ct_c1 += loader & loader3;

  loader = *lptr++;
  loader2 = *maskp++;
  to_ct_a2 = loader & loader2;
  to_ct_b2 = (loader >> 1) & loader2;
  to_ct_c2 = loader & to_ct_b2;
  loader = *lptr++;
  loader2 = *maskp++;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;
  loader = *lptr;
  loader2 = *maskp;
  loader3 = (loader >> 1) & loader2;
  to_ct_a2 += loader & loader2;
  to_ct_b2 += loader3;
  to_ct_c2 += loader & loader3;

  to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
  to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
  partial_a += (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
  to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
  to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
  partial_b += (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
  to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
  to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
  partial_c += (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

  *ctap += (partial_a * 0x01010101) >> 24;
  *ctbp += (partial_b * 0x01010101) >> 24;
  *ctcp += (partial_c * 0x01010101) >> 24;
}

void freq_hwe_haploid_count_12(uintptr_t* lptr, uintptr_t* maskp, uint32_t* ct_nmp, uint32_t* ct_hmajp) {
  uintptr_t loader = *lptr++;
  uintptr_t loader3 = loader >> 1;
  uintptr_t loader2 = loader ^ (~loader3);
  uint32_t to_ct_nm1;
  uint32_t to_ct_hmaj1;
  uint32_t to_ct_nm2;
  uint32_t to_ct_hmaj2;
  uintptr_t partial_nm;
  uintptr_t partial_hmaj;
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm = (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj = (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 = loader2 & loader3;
  to_ct_hmaj1 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm1 += loader2 & loader3;
  to_ct_hmaj1 += loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 = loader2 & loader3;
  to_ct_hmaj2 = loader & loader3;

  loader = *lptr++;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp++;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  loader = *lptr;
  loader3 = loader >> 1;
  loader2 = loader & (~loader3);
  loader &= loader3;
  loader3 = *maskp;
  to_ct_nm2 += loader2 & loader3;
  to_ct_hmaj2 += loader & loader3;

  to_ct_nm1 = (to_ct_nm1 & 0x33333333) + ((to_ct_nm1 >> 2) & 0x33333333);
  to_ct_nm1 += (to_ct_nm2 & 0x33333333) + ((to_ct_nm2 >> 2) & 0x33333333);
  partial_nm += (to_ct_nm1 & 0x0f0f0f0f) + ((to_ct_nm1 >> 4) & 0x0f0f0f0f);
  to_ct_hmaj1 = (to_ct_hmaj1 & 0x33333333) + ((to_ct_hmaj1 >> 2) & 0x33333333);
  to_ct_hmaj1 += (to_ct_hmaj2 & 0x33333333) + ((to_ct_hmaj2 >> 2) & 0x33333333);
  partial_hmaj += (to_ct_hmaj1 & 0x0f0f0f0f) + ((to_ct_hmaj1 >> 4) & 0x0f0f0f0f);

  *ct_nmp += (partial_nm * 0x01010101) >> 24;
  *ct_hmajp += (partial_hmaj * 0x01010101) >> 24;
}
#endif
*/

void single_marker_cc_freqs(uintptr_t indiv_ct, uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
  // Counts the number of A2 alleles and missing calls for both cases and
  // controls, for an autosomal marker.  (The caller is expected to calculate
  // the A1 allele count.)
  // See single_marker_freqs_and_hwe() for discussion.
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & (~B)
  // missing: popcount(C)
  // A2: [popcount(A) + popcount(B)] - popcount(C)
  uint32_t tot_ctrl_ab = 0;
  uint32_t tot_ctrl_c = 0;
  uint32_t tot_case_ab = 0;
  uint32_t tot_case_c = 0;
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = *ctrl_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_ctrl_ab += popcount2_long(loader2 + loader3);
    tot_ctrl_c += popcount2_long(loader2 & (~loader3));
    loader2 = *case_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_case_ab += popcount2_long(loader2 + loader3);
    tot_case_c += popcount2_long(loader2 & (~loader3));
  }
  *ctrl_missingp = tot_ctrl_c;
  *ctrl_setp = tot_ctrl_ab - tot_ctrl_c;
  *case_missingp = tot_case_c;
  *case_setp = tot_case_ab - tot_case_c;
}

void haploid_single_marker_cc_freqs(uintptr_t indiv_ct, uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_setp, uint32_t* ctrl_missingp, uint32_t* case_setp, uint32_t* case_missingp) {
  // Counts the number of A1 and A2 alleles for both cases and controls, for a
  // haploid marker.
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := B ^ A
  // missing: popcount(C)
  // A2: popcount(A & B)
  uint32_t tot_ctrl_ab = 0;
  uint32_t tot_ctrl_c = 0;
  uint32_t tot_case_ab = 0;
  uint32_t tot_case_c = 0;
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  while (lptr < lptr_end) {
    loader = *lptr++;
    loader2 = loader >> 1;
    loader3 = loader2 ^ loader;
    loader &= loader2;
    loader2 = *ctrl_include2++;
    tot_ctrl_ab += popcount2_long(loader & loader2);
    tot_ctrl_c += popcount2_long(loader3 & loader2);
    loader2 = *case_include2++;
    tot_case_ab += popcount2_long(loader & loader2);
    tot_case_c += popcount2_long(loader3 & loader2);
  }
  *ctrl_setp = tot_ctrl_ab;
  *ctrl_missingp = tot_ctrl_c;
  *case_setp = tot_case_ab;
  *case_missingp = tot_case_c;
}

void single_marker_cc_3freqs(uintptr_t indiv_ct, uintptr_t indiv_ctl2, uintptr_t* lptr, uintptr_t* ctrl_include2, uintptr_t* case_include2, uint32_t* ctrl_hom2p, uint32_t* ctrl_hetp, uint32_t* ctrl_missingp, uint32_t* case_hom2p, uint32_t* case_hetp, uint32_t* case_missingp) {
  // Counts the number of heterozygotes, A2 homozygotes, and missing calls for
  // both cases and controls.  Assumes marker is diploid.  The caller is
  // expected to calculate the A1 allele count.
  // See single_marker_freqs_and_hwe() for discussion.
  //   A := genotype & 0x5555...
  //   B := (genotype >> 1) & 0x5555...
  //   C := A & B
  //   popcount(C) = homozyg major ct
  //   popcount(B) = het ct + homozyg major ct
  //   popcount(A) = missing_ct + homozyg major ct
  // hom2: popcount(C)
  // het: popcount(B) - popcount(C)
  // missing: popcount(A) - popcount(C)
  uint32_t tot_ctrl_a = 0;
  uint32_t tot_ctrl_b = 0;
  uint32_t tot_ctrl_c = 0;
  uint32_t tot_case_a = 0;
  uint32_t tot_case_b = 0;
  uint32_t tot_case_c = 0;
  uintptr_t* lptr_end = &(lptr[indiv_ctl2]);
  uintptr_t loader;
  uintptr_t loader2;
  uintptr_t loader3;
  /*
#ifdef __LP64__
  uintptr_t cur_decr;
  uintptr_t* lptr_12x_end;
  unfiltered_indiv_ctl2 -= unfiltered_indiv_ctl2 % 12;
  while (unfiltered_indiv_ctl2 >= 120) {
    cur_decr = 120;
  single_marker_freqs_and_hwe_loop:
    lptr_12x_end = &(lptr[cur_decr]);
    freq_hwe_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)indiv_include2, &tot_a, &tot_b, &tot_c);
    freq_hwe_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_needed) {
      freq_hwe_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[cur_decr]);
      if (hardy_needed) {
	freq_hwe_count_120v((__m128i*)lptr, (__m128i*)lptr_12x_end, (__m128i*)founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[cur_decr]);
      }
    }
    lptr = lptr_12x_end;
    indiv_include2 = &(indiv_include2[cur_decr]);
    founder_include2 = &(founder_include2[cur_decr]);
    unfiltered_indiv_ctl2 -= cur_decr;
  }
  if (unfiltered_indiv_ctl2) {
    cur_decr = unfiltered_indiv_ctl2;
    goto single_marker_freqs_and_hwe_loop;
  }
#else
  uintptr_t* lptr_twelve_end = &(lptr[unfiltered_indiv_ctl2 - unfiltered_indiv_ctl2 % 12]);
  while (lptr < lptr_twelve_end) {
    freq_hwe_count_12(lptr, indiv_include2, &tot_a, &tot_b, &tot_c);
    freq_hwe_count_12(lptr, founder_include2, &tot_a_f, &tot_b_f, &tot_c_f);
    if (hwe_needed) {
      freq_hwe_count_12(lptr, founder_ctrl_include2, &tot_a_hwe, &tot_b_hwe, &tot_c_hwe);
      founder_ctrl_include2 = &(founder_ctrl_include2[12]);
      if (hardy_needed) {
	freq_hwe_count_12(lptr, founder_case_include2, &tot_a_chwe, &tot_b_chwe, &tot_c_chwe);
	founder_case_include2 = &(founder_case_include2[12]);
      }
    }
    lptr = &(lptr[12]);
    indiv_include2 = &(indiv_include2[12]);
    founder_include2 = &(founder_include2[12]);
  }
#endif
  */
  while (lptr < lptr_end) {
    //   A := genotype & 0x5555...
    //   B := (genotype >> 1) & 0x5555...
    //   C := A & B
    //   popcount(C) = homozyg major ct
    //   popcount(B) = het ct + homozyg major ct
    //   popcount(A) = missing_ct + homozyg major ct
    loader = *lptr++;
    loader2 = *ctrl_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_ctrl_a += popcount2_long(loader2);
    tot_ctrl_b += popcount2_long(loader3);
    tot_ctrl_c += popcount2_long(loader2 & loader3);
    loader2 = *case_include2++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    tot_case_a += popcount2_long(loader2);
    tot_case_b += popcount2_long(loader3);
    tot_case_c += popcount2_long(loader2 & loader3);
  }
  *ctrl_hom2p = tot_ctrl_c;
  *ctrl_hetp = tot_ctrl_b - tot_ctrl_c;
  *ctrl_missingp = tot_ctrl_a - tot_ctrl_c;
  *case_hom2p = tot_case_c;
  *case_hetp = tot_case_b - tot_case_c;
  *case_missingp = tot_case_a - tot_case_c;
}

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* wkspace_mark = wkspace_base;
  uintptr_t unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
  uintptr_t pheno_nm_ctl2 = 2 * ((pheno_nm_ct + (BITCT - 1)) / BITCT);
  int32_t retval = 0;
  FILE* outfile = NULL;
  uint32_t model_assoc = model_modifier & MODEL_ASSOC;
  uint32_t model_adapt = model_modifier & MODEL_PERM;
  uint32_t model_maxt = model_modifier & MODEL_MPERM;
  uint32_t model_fisher = model_modifier & MODEL_FISHER;
  uint32_t model_trendonly = model_modifier & MODEL_TRENDONLY;
  uint32_t assoc_counts = model_modifier & MODEL_ASSOC_COUNTS;
  uint32_t assoc_p2 = model_modifier & MODEL_ASSOC_P2;
  uint32_t display_ci = (ci_size > 0);
  // uint32_t perms_done = 0;
  uint32_t perms_left = 0;
  uint32_t perm_pass_ct = 1;
  int32_t x_code = -1;
  int32_t y_code = -1;
  uint32_t case_ct = 0;
  uint32_t male_ct = 0;
  uint32_t nonmale_ct = 0;
  uint32_t ctrl_male_ct = 0;
  uint32_t case_male_ct = 0;
  uint32_t is_y = 0;
  uint32_t ctrl_nonmale_ct = 0;
  uint32_t case_nonmale_ct = 0;
  uint32_t load_ctrl_ct = 0;
  uint32_t load_case_ct = 0;
  uintptr_t* indiv_nonmale_ctrl_include2 = NULL;
  uintptr_t* indiv_nonmale_case_include2 = NULL;
  uintptr_t* indiv_male_ctrl_include2 = NULL;
  uintptr_t* indiv_male_case_include2 = NULL;
  uintptr_t* cur_ctrl_include2 = NULL;
  uintptr_t* cur_case_include2 = NULL;
  double dxx = 0.0;
  double dww = 0.0;
  double dvv = 0.0;
  uint32_t mu_table[MODEL_BLOCKSIZE];
  char wformat[64];
  uintptr_t marker_ctl;
  uint32_t ctrl_ct;
  uint32_t perm_pass_idx;
  uint32_t chrom_fo_idx;
  uint32_t chrom_end;
  uint32_t block_start;
  uint32_t block_size;
  uint32_t block_end;
  uint32_t marker_bidx;
  // uint32_t ld_check_idx;
  uintptr_t marker_uidx;
  uintptr_t marker_uidx2;
  uintptr_t marker_idx;
  uintptr_t indiv_uidx;
  uintptr_t indiv_idx;
  double* orig_stats;
  double* orig_odds;
  double* osptr;
  double* ooptr;
  unsigned char* loadbuf_raw;
  uintptr_t* loadbuf;
  uintptr_t* loadbuf_ptr;
  // uintptr_t* perm_exclude;
  uintptr_t* indiv_ctrl_include2;
  uintptr_t* indiv_case_include2;
  uintptr_t* indiv_nonmale_include2;
  uintptr_t* indiv_male_include2;
  uint32_t load_indiv_ct;
  uint32_t is_x;
  uint32_t is_haploid;
  uint32_t is_reverse;
  uint64_t ullii;
  uintptr_t ulii;
  uint32_t uii;
  uint32_t ujj;
  uint32_t ukk;
  uint32_t umm;
  uint32_t unn;
  uint32_t uoo;
  uint32_t upp;
  uint32_t uqq;
  uint32_t urr;
  uint32_t uss;
  double pval;
  double dyy;
  double dzz;
  double duu;
  double dtt;
  double da1;
  double da12;
  double da2;
  double du1;
  double du12;
  double du2;
  double obs_t;
  double obs_1;
  double obs_2;
  double obs_11;
  double obs_12;
  double obs_22;
  double obs_a1;
  double obs_a2;
  double obs_u1;
  double obs_u2;
  double exp_a1;
  double exp_a12;
  double exp_a2;
  double exp_u1;
  double exp_u12;
  double exp_u2;
  double mult_p;
  // double dom_p;
  // double rec_p;
  // double best_p;
  char* a1ptr;
  char* a2ptr;

  if (assoc_p2) {
    logprint("Error: --assoc p2 not yet supported.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
  }
  if (!pheno_c) {
    logprint("Error: --assoc does not yet support quantitative phenotypes.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
  }
  if (model_maxt) {
    logprint("Error: max(T) permutation is not yet implemented.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
  } else if (model_adapt) {
    logprint("Error: adaptive permutation is not yet implemented.\n");
    retval = RET_CALC_NOT_YET_SUPPORTED;
    goto model_assoc_ret_1;
    perms_left = aperm_max;
  }
  if (wkspace_alloc_uc_checked(&loadbuf_raw, unfiltered_indiv_ct4)) {
    goto model_assoc_ret_NOMEM;
  }
  if (model_assoc) {
    if (model_fisher) {
      memcpy(outname_end, ".assoc.fisher", 14);
    } else {
      memcpy(outname_end, ".assoc", 7);
    }
    if (fopen_checked(&outfile, outname, "w")) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "Writing --assoc report to %s...", outname);
    logprintb();
    fflush(stdout);
    sprintf(tbuf, " CHR %%%us         BP   A1 ", plink_maxsnp);
    if (fprintf(outfile, tbuf, "SNP") < 0) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (assoc_counts) {
      if (fputs("     C_A      C_U   A2 ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    } else {
      if (fputs("     F_A      F_U   A2 ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    }
    if (!model_fisher) {
      if (fputs("       CHISQ ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    }
    if (fputs("           P           OR ", outfile) == EOF) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (display_ci) {
      uii = (uint32_t)(ci_size * 100);
      if (uii >= 10) {
	if (fprintf(outfile, "          SE          L%u          U%u ", uii, uii) < 0) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
      } else {
	if (fprintf(outfile, "          SE           L%u           U%u ", uii, uii) < 0) {
	  goto model_assoc_ret_WRITE_FAIL;
	}
      }
    }
    if (putc('\n', outfile) == EOF) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (max_marker_allele_len == 1) {
      sprintf(wformat, "     %%%us %%10u    %%c ", plink_maxsnp);
    } else {
      sprintf(wformat, "     %%%us %%10u %%4s ", plink_maxsnp);
    }
  } else {
    uii = count_non_autosomal_markers(chrom_info_ptr, marker_exclude, 0);
    if (uii) {
      sprintf(logbuf, "Excluding %u haploid marker%s from --model analysis.\n", uii, (uii == 1)? "" : "s");
      logprintb();
    }
    marker_ct -= uii;
    if (!marker_ct) {
      logprint("Error: No markers remaining for --model analysis.\n");
      retval = RET_INVALID_CMDLINE;
      goto model_assoc_ret_1;
    }
    memcpy(outname_end, ".model", 7);
    if (fopen_checked(&outfile, outname, "w")) {
      goto model_assoc_ret_OPEN_FAIL;
    }
    sprintf(logbuf, "Writing --model report to %s...", outname);
    logprintb();
    fflush(stdout);
    sprintf(tbuf, " CHR %%%us   A1   A2     TEST            AFF          UNAFF ", plink_maxsnp);
    if (fprintf(outfile, tbuf, "SNP") < 0) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (!model_fisher) {
      if (fputs("       CHISQ   DF ", outfile) == EOF) {
	goto model_assoc_ret_WRITE_FAIL;
      }
    }
    if (fputs("           P\n", outfile) == EOF) {
      goto model_assoc_ret_WRITE_FAIL;
    }
    if (max_marker_allele_len == 1) {
      sprintf(wformat, "     %%%us    %%c    %%c  ", plink_maxsnp);
    } else {
      sprintf(wformat, "     %%%us %%4s %%4s  ", plink_maxsnp);
    }
  }
  if (wkspace_alloc_d_checked(&orig_stats, marker_ct * sizeof(double)) ||
      wkspace_alloc_d_checked(&orig_odds, marker_ct * sizeof(double))) {
    goto model_assoc_ret_NOMEM;
  }
  if (model_adapt || model_maxt) {
    // indiv_ctrl_include2, etc. point to first permutation slot instead of
    // getting their own memory
  } else {
    if (pheno_c) {
      if (wkspace_alloc_ul_checked(&indiv_ctrl_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	  wkspace_alloc_ul_checked(&indiv_case_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
	goto model_assoc_ret_NOMEM;
      }
      fill_ulong_zero(indiv_case_include2, pheno_nm_ctl2);
      indiv_uidx = 0;
      for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
	indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	if (is_set(pheno_c, indiv_uidx)) {
	  set_bit_noct(indiv_case_include2, indiv_idx * 2);
	  case_ct++;
	}
	indiv_uidx++;
      }
      vec_init_invert(pheno_nm_ct, indiv_ctrl_include2, indiv_case_include2);
      ctrl_ct = pheno_nm_ct - case_ct;
    }
    x_code = species_x_code[chrom_info_ptr->species];
    y_code = species_y_code[chrom_info_ptr->species];
    ullii = chrom_info_ptr->chrom_mask;
    if (((x_code != -1) && (ullii & (1LLU << x_code))) || (model_assoc && (((y_code != -1) && (ullii & (1LLU << y_code)))))) {
      if (wkspace_alloc_ul_checked(&indiv_nonmale_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	  wkspace_alloc_ul_checked(&indiv_male_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
	goto model_assoc_ret_NOMEM;
      }
      fill_ulong_zero(indiv_male_include2, pheno_nm_ctl2);
      if (pheno_c) {
	// argh
	if (wkspace_alloc_ul_checked(&indiv_nonmale_ctrl_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	    wkspace_alloc_ul_checked(&indiv_nonmale_case_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	    wkspace_alloc_ul_checked(&indiv_male_ctrl_include2, pheno_nm_ctl2 * sizeof(intptr_t)) ||
	    wkspace_alloc_ul_checked(&indiv_male_case_include2, pheno_nm_ctl2 * sizeof(intptr_t))) {
	  goto model_assoc_ret_NOMEM;
	}
	fill_ulong_zero(indiv_male_case_include2, pheno_nm_ctl2);
	indiv_uidx = 0;
	for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
	  indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	  if (is_set(sex_nm, indiv_uidx) && is_set(sex_male, indiv_uidx)) {
	    set_bit_noct(indiv_male_include2, indiv_idx * 2);
	    male_ct++;
	    if (is_set(pheno_c, indiv_uidx)) {
	      set_bit_noct(indiv_male_case_include2, indiv_idx * 2);
	      case_male_ct++;
	    }
	  }
	  indiv_uidx++;
	}
	vec_init_andnot(pheno_nm_ctl2, indiv_male_ctrl_include2, indiv_male_include2, indiv_male_case_include2);
	vec_init_invert(pheno_nm_ct, indiv_nonmale_include2, indiv_male_include2);
	vec_init_andnot(pheno_nm_ctl2, indiv_nonmale_case_include2, indiv_case_include2, indiv_male_case_include2);
	vec_init_andnot(pheno_nm_ctl2, indiv_nonmale_ctrl_include2, indiv_ctrl_include2, indiv_male_ctrl_include2);
	ctrl_male_ct = male_ct - case_male_ct;
	case_nonmale_ct = case_ct - case_male_ct;
	ctrl_nonmale_ct = ctrl_ct - ctrl_male_ct;
      } else {
	indiv_uidx = 0;
	for (indiv_idx = 0; indiv_idx < pheno_nm_ct; indiv_idx++) {
	  indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	  if (is_set(sex_nm, indiv_uidx) && is_set(sex_male, indiv_uidx)) {
	    set_bit_noct(indiv_male_include2, indiv_idx * 2);
	    male_ct++;
	  }
	  indiv_uidx++;
	}
	vec_init_invert(pheno_nm_ct, indiv_nonmale_include2, indiv_male_include2);
      }
      nonmale_ct = pheno_nm_ct - male_ct;
    }
  }

  if (wkspace_alloc_ul_checked(&loadbuf, MODEL_BLOCKSIZE * pheno_nm_ctl2 * sizeof(intptr_t))) {
    goto model_assoc_ret_NOMEM;
  }
  marker_ctl = (marker_ct + (BITCT - 1)) / BITCT;
  if (model_maxt || model_adapt) {
    // check if enough memory to process all permutations in one pass
    // if yes, marker-specific information arrays only need to be of size
    // MODEL_BLOCKSIZE
    // if no, they have to be of size marker_ct...
  }
  for (uii = 0; uii < MODEL_BLOCKSIZE; uii++) {
    loadbuf[(uii + 1) * pheno_nm_ctl2 - 2] = 0;
    loadbuf[(uii + 1) * pheno_nm_ctl2 - 1] = 0;
  }
  marker_uidx = next_non_set_unsafe(marker_exclude, 0);
  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
    goto model_assoc_ret_READ_FAIL;
  }
  for (perm_pass_idx = 0; perm_pass_idx < perm_pass_ct; perm_pass_idx++) {
    chrom_fo_idx = 0xffffffffU;
    marker_uidx = 0;
    marker_idx = 0;
    chrom_end = 0;
    do {
      if (marker_uidx >= chrom_end) {
	block_start = 0;
	if (model_assoc) {
	  // exploit overflow
	  chrom_fo_idx++;
	  refresh_chrom_info(chrom_info_ptr, marker_uidx, 1, 0, &chrom_end, &chrom_fo_idx, &is_x, &is_haploid);
	  uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	  is_y = (uii == (uint32_t)y_code);
	  if (is_haploid && (!is_x)) {
	    if (is_y) {
	      cur_ctrl_include2 = indiv_male_ctrl_include2;
	      cur_case_include2 = indiv_male_case_include2;
	      load_indiv_ct = male_ct;
	      load_case_ct = case_male_ct;
	    } else {
	      cur_ctrl_include2 = indiv_ctrl_include2;
	      cur_case_include2 = indiv_case_include2;
	      load_indiv_ct = pheno_nm_ct;
	      load_case_ct = case_ct;
	    }
	    load_ctrl_ct = load_indiv_ct - load_case_ct;
	  }
	} else {
	  while (1) {
	    do {
	      chrom_end = chrom_info_ptr->chrom_file_order_marker_idx[(++chrom_fo_idx) + 1U];
	    } while (marker_uidx >= chrom_end);
	    uii = chrom_info_ptr->chrom_file_order[chrom_fo_idx];
	    is_x = (uii == (uint32_t)x_code);
	    if ((!((species_haploid_mask[chrom_info_ptr->species] >> uii) & 1LLU)) || is_x) {
	      break;
	    }
	    marker_uidx = next_non_set_unsafe(marker_exclude, chrom_end);
	  }
	  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto model_assoc_ret_READ_FAIL;
	  }
	  if (!is_x) {
	    cur_ctrl_include2 = indiv_ctrl_include2;
	    cur_case_include2 = indiv_case_include2;
	    load_indiv_ct = pheno_nm_ct;
	    load_case_ct = case_ct;
	  } else {
	    cur_ctrl_include2 = indiv_nonmale_ctrl_include2;
	    cur_case_include2 = indiv_nonmale_case_include2;
	    load_indiv_ct = nonmale_ct;
	    load_case_ct = case_nonmale_ct;
	  }
	  load_ctrl_ct = load_indiv_ct - load_case_ct;
	}
	intprint2(&(wformat[2]), uii);
      } else if (model_maxt) {
	// todo: copy LD table, etc.
	// could conditionally do this for adaptive permutation, but I doubt
	// it's worth the trouble
      } else {
	block_start = 0;
      }
      block_size = block_start;
      block_end = marker_ct + block_start - marker_idx;
      if (block_end > MODEL_BLOCKSIZE) {
	block_end = MODEL_BLOCKSIZE;
      }
      do {
	if (fread(loadbuf_raw, 1, unfiltered_indiv_ct4, bedfile) < unfiltered_indiv_ct4) {
	  goto model_assoc_ret_READ_FAIL;
	}
	indiv_uidx = 0;
	ukk = pheno_nm_ct / BITCT2;
	loadbuf_ptr = &(loadbuf[block_size * pheno_nm_ctl2]);
	for (uii = 0; uii < ukk; uii++) {
	  ulii = 0;
	  for (ujj = 0; ujj < BITCT; ujj += 2) {
	    indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	    ulii |= ((uintptr_t)(((loadbuf_raw[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3))) << ujj;
	    indiv_uidx++;
	  }
	  *loadbuf_ptr++ = ulii;
	}
	ujj = 2 * (pheno_nm_ct & (BITCT2 - 1));
	if (ujj) {
	  ulii = 0;
	  for (uii = 0; uii < ujj; uii += 2) {
	    indiv_uidx = next_set_unsafe(pheno_nm, indiv_uidx);
	    ulii |= ((uintptr_t)(((loadbuf_raw[indiv_uidx / 4] >> ((indiv_uidx & 3) * 2)) & 3))) << uii;
	    indiv_uidx++;
	  }
	  *loadbuf_ptr = ulii;
	}
	mu_table[block_size] = marker_uidx;
	if (is_set(marker_exclude, ++marker_uidx)) {
	  marker_uidx = next_non_set_unsafe(marker_exclude, marker_uidx + 1);
	  if (fseeko(bedfile, bed_offset + (uint64_t)marker_uidx * unfiltered_indiv_ct4, SEEK_SET)) {
	    goto model_assoc_ret_READ_FAIL;
	  }
	}
      } while ((++block_size < block_end) && (marker_uidx < chrom_end));
      if (!perm_pass_idx) {
	// basic --assoc/--model
	osptr = &(orig_stats[marker_idx + block_start]);
	ooptr = &(orig_odds[marker_idx + block_start]);
	if (pheno_c) {
	  if (model_assoc) {
	    for (marker_bidx = block_start; marker_bidx < block_size; marker_bidx++) {
	      marker_uidx2 = mu_table[marker_bidx];
	      is_reverse = is_set(marker_reverse, marker_uidx2);
	      if (is_x) {
		if (is_reverse) {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_nonmale_ctrl_include2, indiv_nonmale_case_include2, &ujj, &uii, &umm, &ukk);
		  uii = 2 * (ctrl_nonmale_ct - uii) - ujj;
		  ukk = 2 * (case_nonmale_ct - ukk) - umm;
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_male_ctrl_include2, indiv_male_case_include2, &uoo, &unn, &uqq, &upp);
		  unn = ctrl_male_ct - unn - uoo;
		  upp = case_male_ct - upp - uqq;
		} else {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_nonmale_ctrl_include2, indiv_nonmale_case_include2, &uii, &ujj, &ukk, &umm);
		  ujj = 2 * (ctrl_nonmale_ct - ujj) - uii;
		  umm = 2 * (case_nonmale_ct - umm) - ukk;
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_male_ctrl_include2, indiv_male_case_include2, &unn, &uoo, &upp, &uqq);
		  uoo = ctrl_male_ct - uoo - unn;
		  uqq = case_male_ct - uqq - upp;
		}
		uii += unn;
		ujj += uoo;
		ukk += upp;
		umm += uqq;
	      } else if (is_haploid) {
		if (is_reverse) {
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &ujj, &uii, &umm, &ukk);
		  uii = load_ctrl_ct - uii - ujj;
		  ukk = load_case_ct - ukk - umm;
		} else {
		  haploid_single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm);
		  ujj = load_ctrl_ct - ujj - uii;
		  umm = load_case_ct - umm - ukk;
		}
	      } else {
		if (is_reverse) {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &ujj, &uii, &umm, &ukk);
		  uii = 2 * (ctrl_ct - uii) - ujj;
		  ukk = 2 * (case_ct - ukk) - umm;
		} else {
		  single_marker_cc_freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), indiv_ctrl_include2, indiv_case_include2, &uii, &ujj, &ukk, &umm);
		  ujj = 2 * (ctrl_ct - ujj) - uii;
		  umm = 2 * (case_ct - umm) - ukk;
		}
	      }
	      da1 = umm;
	      da2 = ukk;
	      du1 = ujj;
	      du2 = uii;
	      if (model_fisher) {
		pval = fisher22(uii, ujj, ukk, umm);
		*osptr = 1 - pval;
	      } else if (!assoc_p2) {
		if ((((uint64_t)umm) * (ujj + uii)) == (((uint64_t)ujj) * (umm + ukk))) {
		  *osptr = 0;
		  if ((!umm) && (!ujj)) {
		    pval = -1;
		  } else {
		    pval = chiprob_p(*osptr, 1);
		  }
		} else {
		  dxx = 1.0 / ((double)(da1 + da2 + du1 + du2));
		  dzz = (da1 + du1) * dxx;
		  dvv = (da2 + du2) * dxx;
		  dyy = (da1 + da2) * dzz;
		  dzz *= du1 + du2;
		  dww = (da1 + da2) * dvv;
		  dvv *= du1 + du2;
		  *osptr = (da1 - dyy) * (da1 - dyy) / dyy + (du1 - dzz) * (du1 - dzz) / dzz + (da2 - dww) * (da2 - dww) / dww + (du2 - dvv) * (du2 - dvv) / dvv;
		  pval = chiprob_p(*osptr, 1);
		}
	      } else {
		// todo
	      }
	      *ooptr = (da1 * du2) / (du1 * da2);
	      if ((pval <= pfilter) || (pfilter == 1.0)) {
		a1ptr = &(marker_alleles[(2 * marker_uidx2 + is_reverse) * max_marker_allele_len]);
		a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1 - is_reverse) * max_marker_allele_len]);
		if (max_marker_allele_len == 1) {
		  if (fprintf(outfile, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), marker_pos[marker_uidx2], *a1ptr) < 0) {
		    goto model_assoc_ret_WRITE_FAIL;
		  }
		  if (umm + ukk) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u ", umm);
		    } else {
		      fprintf(outfile, "%8.4g ", da1 / (da1 + da2));
		    }
		  } else {
		    fputs("      NA ", outfile);
		  }
		  if (ujj + uii) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u    ", ujj);
		    } else {
		      fprintf(outfile, "%8.4g    ", du1 / (du1 + du2));
		    }
		  } else {
		    fputs("      NA    ", outfile);
		  }
		  putc(*a2ptr, outfile);
		} else {
		  if (fprintf(outfile, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), marker_pos[marker_uidx2], a1ptr) < 0) {
		    goto model_assoc_ret_WRITE_FAIL;
		  }
		  if (umm + ukk) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u ", umm);
		    } else {
		      fprintf(outfile, "%8.4g ", da1 / (da1 + da2));
		    }
		  } else {
		    fputs("      NA ", outfile);
		  }
		  if (ujj + uii) {
		    if (assoc_counts) {
		      fprintf(outfile, "%8u", ujj);
		    } else {
		      fprintf(outfile, "%8.4g", du1 / (du1 + du2));
		    }
		  } else {
		    fputs("      NA", outfile);
		  }
		  fprintf(outfile, " %4s", a2ptr);
		}
		if (model_fisher) {
		  fprintf(outfile, " %12.4g ", pval);
		} else {
		  if (pval > -1) {
		    fprintf(outfile, " %12.4g %12.4g ", *osptr, pval);
		  } else {
		    fputs("           NA           NA ", outfile);
		  }
		}
		if (du1 * da2 == 0.0) {
		  fputs("          NA ", outfile);
		  if (display_ci) {
		    fputs("          NA           NA           NA ", outfile);
		  }
		} else {
		  fprintf(outfile, "%12.4g ", *ooptr);
		  if (display_ci) {
		    dxx = log(*ooptr);
		    dyy = sqrt(1 / da1 + 1 / da2 + 1 / du1 + 1 / du2);
		    dzz = ci_zt * dyy;
		    dww = exp(dxx - dzz);
		    dvv = exp(dxx + dzz);
		    fprintf(outfile, "%12.4g %12.4g %12.4g ", dyy, dww, dvv);
		  }
		}
		if (putc('\n', outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
	      }
	      osptr++;
	      ooptr++;
	    }
	  } else {
	    for (marker_bidx = block_start; marker_bidx < block_size; marker_bidx++) {
	      marker_uidx2 = mu_table[marker_bidx];
	      is_reverse = is_set(marker_reverse, marker_uidx2);
	      if (is_reverse) {
		single_marker_cc_3freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &ukk, &ujj, &uii, &uoo, &unn, &umm);
		uii = load_ctrl_ct - uii - ujj - ukk;
		umm = load_case_ct - umm - unn - uoo;
	      } else {
		single_marker_cc_3freqs(pheno_nm_ct, pheno_nm_ctl2, &(loadbuf[marker_bidx * pheno_nm_ctl2]), cur_ctrl_include2, cur_case_include2, &uii, &ujj, &ukk, &umm, &unn, &uoo);
		ukk = load_ctrl_ct - uii - ujj - ukk;
		uoo = load_case_ct - umm - unn - uoo;
	      }
	      da1 = uoo;
	      da12 = unn;
	      da2 = umm;
	      du1 = ukk;
	      du12 = ujj;
	      du2 = uii;
	      // invalid?
	      upp = (uoo < model_cell_ct) || (unn < model_cell_ct) || (umm < model_cell_ct) || (ukk < model_cell_ct) || (ujj < model_cell_ct) || (uii < model_cell_ct);

	      dxx = da1 + da12 + da2; // obs_A
	      dyy = du1 + du12 + du2; // obs_U
	      obs_t = dxx + dyy;
	      dzz = 1.0 / obs_t;
	      obs_11 = da1 + du1;
	      obs_12 = da12 + du12;
	      obs_22 = da2 + du2;
	      obs_1 = 2 * obs_11 + obs_12;
	      obs_2 = 2 * obs_22 + obs_12;

	      a1ptr = &(marker_alleles[(2 * marker_uidx2 + is_reverse) * max_marker_allele_len]);
	      a2ptr = &(marker_alleles[(2 * marker_uidx2 + 1 - is_reverse) * max_marker_allele_len]);
	      if (max_marker_allele_len == 1) {
		uss = sprintf(tbuf, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), *a1ptr, *a2ptr);
	      } else {
		uss = sprintf(tbuf, wformat, &(marker_ids[marker_uidx2 * max_marker_id_len]), a1ptr, a2ptr);
	      }
	      if (!model_trendonly) {
		memcpy(&(tbuf[uss]), "   GENO          ", 17);
		uqq = uss + 8;
		urr = sprintf(&(tbuf[uqq + 43]), "%u/%u/%u          ", uoo, unn, umm);
		if (urr < 24) {
		  memcpy(&(tbuf[uqq + 24 - urr]), &(tbuf[uqq + 43]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 43]), urr);
		  uqq += urr - 9;
		}
		urr = sprintf(&(tbuf[uqq + 34]), "%u/%u/%u ", ukk, ujj, uii);
		if (urr < 15) {
		  memcpy(&(tbuf[uqq + 15 - urr]), &(tbuf[uqq + 34]), urr);
		  uqq += 15;
		} else {
		  memcpy(&(tbuf[uqq]), &(tbuf[uqq + 34]), urr);
		  uqq += urr;
		}
		if (upp) {
		  dww = -9;
		} else {
		  if (model_fisher) {
		    dww = fisher23(uii, ujj, ukk, umm, unn, uoo);
		  } else {
		    dvv = dzz * dxx;
		    exp_a1 = dvv * obs_11;
		    exp_a12 = dvv * obs_12;
		    exp_a2 = dvv * obs_22;
		    dvv = dzz * dyy;
		    exp_u1 = dvv * obs_11;
		    exp_u12 = dvv * obs_12;
		    exp_u2 = dvv * obs_22;
		    dvv = ((da1 - exp_a1) * (da1 - exp_a1)) / exp_a1 +
		      ((da12 - exp_a12) * (da12 - exp_a12)) / exp_a12 +
		      ((da2 - exp_a2) * (da2 - exp_a2)) / exp_a2 +
		      ((du1 - exp_u1) * (du1 - exp_u1)) / exp_u1 +
		      ((du12 - exp_u12) * (du12 - exp_u12)) / exp_u12 +
		      ((du2 - exp_u2) * (du2 - exp_u2)) / exp_u2;
		    dww = chiprob_p(dvv, 2);
		  }
		}
		if (dww < -1) {
		  if (model_fisher) {
		    memcpy(&(tbuf[uqq]), "          NA\n", 14);
		  } else {
		    memcpy(&(tbuf[uqq]), "          NA   NA           NA\n", 32);
		  }
		} else {
		  if (model_fisher) {
		    sprintf(&(tbuf[uqq]), "%12.4g\n", dww);
		  } else {
		    sprintf(&(tbuf[uqq]), "%12.4g    2 %12.4g\n", dvv, dww);
		  }
		}
		if (fputs(tbuf, outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
	      }
	      memcpy(&(tbuf[uss]), "  TREND            ", 19);
	      uqq = uss + 8;
	      urr = sprintf(&(tbuf[uqq + 34]), "%u/%u            ", uoo * 2 + unn, umm * 2 + unn);
	      if (urr < 26) {
		memcpy(&(tbuf[uqq + 26 - urr]), &(tbuf[uqq + 34]), urr);
		uqq += 15;
	      } else {
		memcpy(&(tbuf[uqq]), &(tbuf[uqq + 34]), urr);
		uqq += urr - 11;
	      }
	      urr = sprintf(&(tbuf[uqq + 24]), "%u/%u ", ukk * 2 + ujj, uii * 2 + ujj);
	      if (urr < 15) {
		memcpy(&(tbuf[uqq + 15 - urr]), &(tbuf[uqq + 24]), urr);
		uqq += 15;
	      } else {
		memcpy(&(tbuf[uqq]), &(tbuf[uqq + 24]), urr);
		uqq += urr;
	      }
	      urr = uqq; // save this for next line
	      if ((((uint64_t)(uoo * 2 + unn)) * (uii * 2 + ujj)) == (((uint64_t)(umm * 2 + unn)) * (ukk * 2 + ujj))) {
		if (uoo * 2 + unn) {
		  duu = 0;
		  dtt = 1;
		} else {
		  dtt = -9;
		}
	      } else {
		// CA
		dww = (dyy * dzz * da12 - dxx * dzz * du12) + 2 * (dyy * dzz * da2 - dxx * dzz * du2);
		// varCA
		dvv = dxx * dyy * (obs_t * (obs_12 + 4 * obs_22) - (obs_12 + 2 * obs_22) * (obs_12 + 2 * obs_22)) * dzz * dzz * dzz;
		// CA_chisq
		duu = (dww * dww) / dvv;
		// CA_p
		dtt = chiprob_p(duu, 1);
	      }
	      if (dtt > -1) {
		if (model_fisher) {
		  sprintf(&(tbuf[uqq]), "%12.4g\n", dtt);
		} else {
		  sprintf(&(tbuf[uqq]), "%12.4g    1 %12.4g\n", duu, dtt);
		}
	      } else {
		if (!model_fisher) {
		  memcpy(&(tbuf[uqq]), "          NA   NA ", 18);
		  uqq += 18;
		}
		memcpy(&(tbuf[uqq]), "          NA\n", 14);
	      }
	      if (fputs(tbuf, outfile) == EOF) {
		goto model_assoc_ret_WRITE_FAIL;
	      }
	      if (!model_trendonly) {
		memcpy(&(tbuf[uss]), "ALLELIC", 7);
		uqq = urr;
		obs_a1 = 2 * da1 + da12;
		obs_a2 = 2 * da2 + da12;
		obs_u1 = 2 * du1 + du12;
		obs_u2 = 2 * du2 + du12;
		if (model_fisher) {
		  mult_p = fisher22(2 * uoo + unn, 2 * umm + unn, 2 * ukk + ujj, 2 * uii + ujj);
		} else {
		  exp_a1 = dxx * obs_1 * dzz;
		  exp_a2 = dxx * obs_2 * dzz;
		  exp_u1 = dyy * obs_1 * dzz;
		  exp_u2 = dyy * obs_2 * dzz;
		  dww = ((obs_a1 - exp_a1) * (obs_a1 - exp_a1)) / exp_a1 + ((obs_a2 - exp_a2) * (obs_a2 - exp_a2)) / exp_a2 + ((obs_u1 - exp_u1) * (obs_u1 - exp_u1)) / exp_u1 + ((obs_u2 - exp_u2) * (obs_u2 - exp_u2)) / exp_u2;
		  mult_p = chiprob_p(dww, 1);
		}
		if (mult_p > -1) {
		  if (model_fisher) {
		    sprintf(&(tbuf[uqq]), "%12.4g\n", mult_p);
		  } else {
		    sprintf(&(tbuf[uqq]), "%12.4g    1 %12.4g\n", dww, mult_p);
		  }
		} else {
		  if (!model_fisher) {
		    memcpy(&(tbuf[uqq]), "          NA   NA ", 18);
		    uqq += 18;
		  }
		  memcpy(&(tbuf[uqq]), "          NA\n", 14);
		}
		if (fputs(tbuf, outfile) == EOF) {
		  goto model_assoc_ret_WRITE_FAIL;
		}
	      }
	    }
	  }
	}
      }
      marker_idx += block_size - block_start;
    } while (marker_idx < marker_ct);
    if (!perm_pass_idx) {
      logprint(" done.\n");
    }
  }

  while (0) {
  model_assoc_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  model_assoc_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  model_assoc_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  model_assoc_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
 model_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
