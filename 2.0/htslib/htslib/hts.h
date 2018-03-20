/// @file htslib/hts.h
/// Format-neutral I/O, indexing, and iterator API functions.
/*
    Copyright (C) 2012-2016 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.
    Portions copyright (C) 2003-2006, 2008-2010 by Heng Li <lh3@live.co.uk>

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_HTS_H
#define HTSLIB_HTS_H

#include <stddef.h>
#include <stdint.h>

#include "hts_defs.h"
#include "hts_log.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HTS_BGZF_TYPEDEF
typedef struct BGZF BGZF;
#define HTS_BGZF_TYPEDEF
#endif
struct hFILE;
struct hts_tpool;

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/************
 * Indexing *
 ************/

#define HTS_FMT_CSI 0
#define HTS_FMT_BAI 1
#define HTS_FMT_TBI 2
#define HTS_FMT_CRAI 3

struct __hts_idx_t;
typedef struct __hts_idx_t hts_idx_t;

typedef struct {
    uint64_t u, v;
} hts_pair64_t;

typedef struct {
    uint64_t u, v;
    uint64_t max;
} hts_pair64_max_t;

    #define hts_bin_first(l) (((1<<(((l)<<1) + (l))) - 1) / 7)
    #define hts_bin_parent(l) (((l) - 1) >> 3)

    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls);
    void hts_idx_destroy(hts_idx_t *idx);
    int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped);
    void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset);

/// Save an index to a file
/** @param idx  Index to be written
    @param fn   Input BAM/BCF/etc filename, to which .bai/.csi/etc will be added
    @param fmt  One of the HTS_FMT_* index formats
    @return  0 if successful, or negative if an error occurred.
*/
int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt) HTS_RESULT_USED;

/// Save an index to a specific file
/** @param idx    Index to be written
    @param fn     Input BAM/BCF/etc filename
    @param fnidx  Output filename, or NULL to add .bai/.csi/etc to @a fn
    @param fmt    One of the HTS_FMT_* index formats
    @return  0 if successful, or negative if an error occurred.
*/
int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt) HTS_RESULT_USED;


static inline int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
{
    int l, s = min_shift, t = ((1<<((n_lvls<<1) + n_lvls)) - 1) / 7;
    for (--end, l = n_lvls; l > 0; --l, s += 3, t -= 1<<((l<<1)+l))
        if (beg>>s == end>>s) return t + (beg>>s);
    return 0;
}

static inline int hts_bin_bot(int bin, int n_lvls)
{
    int l, b;
    for (l = 0, b = bin; b; ++l, b = hts_bin_parent(b)); // compute the level of bin
    return (bin - hts_bin_first(l)) << (n_lvls - l) * 3;
}

/**************
 * Endianness *
 **************/

static inline int ed_is_big(void)
{
    long one= 1;
    return !(*((char *)(&one)));
}
static inline uint16_t ed_swap_2(uint16_t v)
{
    return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *ed_swap_2p(void *x)
{
    *(uint16_t*)x = ed_swap_2(*(uint16_t*)x);
    return x;
}
static inline uint32_t ed_swap_4(uint32_t v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *ed_swap_4p(void *x)
{
    *(uint32_t*)x = ed_swap_4(*(uint32_t*)x);
    return x;
}
static inline uint64_t ed_swap_8(uint64_t v)
{
    v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
    return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *ed_swap_8p(void *x)
{
    *(uint64_t*)x = ed_swap_8(*(uint64_t*)x);
    return x;
}

#ifdef __cplusplus
}
#endif

#endif
