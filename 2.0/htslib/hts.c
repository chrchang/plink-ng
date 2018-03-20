/*  hts.c -- format-neutral I/O, indexing, and iterator API functions.

    Copyright (C) 2008, 2009, 2012-2017 Genome Research Ltd.
    Copyright (C) 2012, 2013 Broad Institute.

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

// This just includes support for logging and index writing; the rest of the
// code has been removed.

#include <stdio.h>
#include <errno.h>
#include <stdarg.h>

#include "htslib/hts.h"
#include "htslib/bgzf.h"

#include "htslib/khash.h"
#include "htslib/ksort.h"

int hts_verbose = HTS_LOG_WARNING;

/****************
 *** Indexing ***
 ****************/

#define HTS_MIN_MARKER_DIST 0x10000

// Finds the special meta bin
//  ((1<<(3 * n_lvls + 3)) - 1) / 7 + 1
#define META_BIN(idx) ((idx)->n_bins + 1)

#define pair64_lt(a,b) ((a).u < (b).u)

KSORT_INIT(_off, hts_pair64_t, pair64_lt)
KSORT_INIT(_off_max, hts_pair64_max_t, pair64_lt)

typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    int32_t n, m;
    uint64_t *offset;
} lidx_t;

struct __hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    struct {
        uint32_t last_bin, save_bin;
        int last_coor, last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};

static char * idx_format_name(int fmt) {
    switch (fmt) {
        case HTS_FMT_CSI: return "csi";
        case HTS_FMT_BAI: return "bai";
        case HTS_FMT_TBI: return "tbi";
        case HTS_FMT_CRAI: return "crai";
        default: return "unknown";
    }
}

static inline int insert_to_b(bidx_t *b, int bin, uint64_t beg, uint64_t end)
{
    khint_t k;
    bins_t *l;
    int absent;
    k = kh_put(bin, b, bin, &absent);
    if (absent < 0) return -1; // Out of memory
    l = &kh_value(b, k);
    if (absent) {
        l->m = 1; l->n = 0;
        l->list = (hts_pair64_t*)calloc(l->m, sizeof(hts_pair64_t));
        if (!l->list) {
            kh_del(bin, b, k);
            return -1;
        }
    } else if (l->n == l->m) {
        uint32_t new_m = l->m ? l->m << 1 : 1;
        hts_pair64_t *new_list = realloc(l->list, new_m * sizeof(hts_pair64_t));
        if (!new_list) return -1;
        l->list = new_list;
        l->m = new_m;
    }
    l->list[l->n].u = beg;
    l->list[l->n++].v = end;
    return 0;
}

static inline int insert_to_l(lidx_t *l, int64_t _beg, int64_t _end, uint64_t offset, int min_shift)
{
    int i, beg, end;
    beg = _beg >> min_shift;
    end = (_end - 1) >> min_shift;
    if (l->m < end + 1) {
        size_t new_m = l->m * 2 > end + 1 ? l->m * 2 : end + 1;
        uint64_t *new_offset;

        new_offset = (uint64_t*)realloc(l->offset, new_m * sizeof(uint64_t));
        if (!new_offset) return -1;

        // fill unused memory with (uint64_t)-1
        memset(new_offset + l->m, 0xff, sizeof(uint64_t) * (new_m - l->m));
        l->m = new_m;
        l->offset = new_offset;
    }
    for (i = beg; i <= end; ++i) {
        if (l->offset[i] == (uint64_t)-1) l->offset[i] = offset;
    }
    if (l->n < end + 1) l->n = end + 1;
    return 0;
}

hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls)
{
    hts_idx_t *idx;
    idx = (hts_idx_t*)calloc(1, sizeof(hts_idx_t));
    if (idx == NULL) return NULL;
    idx->fmt = fmt;
    idx->min_shift = min_shift;
    idx->n_lvls = n_lvls;
    idx->n_bins = ((1<<(3 * n_lvls + 3)) - 1) / 7;
    idx->z.save_bin = idx->z.save_tid = idx->z.last_tid = idx->z.last_bin = 0xffffffffu;
    idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
    idx->z.last_coor = 0xffffffffu;
    if (n) {
        idx->n = idx->m = n;
        idx->bidx = (bidx_t**)calloc(n, sizeof(bidx_t*));
        if (idx->bidx == NULL) { free(idx); return NULL; }
        idx->lidx = (lidx_t*) calloc(n, sizeof(lidx_t));
        if (idx->lidx == NULL) { free(idx->bidx); free(idx); return NULL; }
    }
    return idx;
}

static void update_loff(hts_idx_t *idx, int i, int free_lidx)
{
    bidx_t *bidx = idx->bidx[i];
    lidx_t *lidx = &idx->lidx[i];
    khint_t k;
    int l;
    uint64_t offset0 = 0;
    if (bidx) {
        k = kh_get(bin, bidx, META_BIN(idx));
        if (k != kh_end(bidx))
            offset0 = kh_val(bidx, k).list[0].u;
        for (l = 0; l < lidx->n && lidx->offset[l] == (uint64_t)-1; ++l)
            lidx->offset[l] = offset0;
    } else l = 1;
    for (; l < lidx->n; ++l) // fill missing values
        if (lidx->offset[l] == (uint64_t)-1)
            lidx->offset[l] = lidx->offset[l-1];
    if (bidx == 0) return;
    for (k = kh_begin(bidx); k != kh_end(bidx); ++k) // set loff
        if (kh_exist(bidx, k))
        {
            if ( kh_key(bidx, k) < ((unsigned int)idx->n_bins) )
            {
                int bot_bin = hts_bin_bot(kh_key(bidx, k), idx->n_lvls);
                // disable linear index if bot_bin out of bounds
                kh_val(bidx, k).loff = bot_bin < lidx->n ? lidx->offset[bot_bin] : 0;
            }
            else
                kh_val(bidx, k).loff = 0;
        }
    if (free_lidx) {
        free(lidx->offset);
        lidx->m = lidx->n = 0;
        lidx->offset = 0;
    }
}

static void compress_binning(hts_idx_t *idx, int i)
{
    bidx_t *bidx = idx->bidx[i];
    khint_t k;
    int l, m;
    if (bidx == 0) return;
    // merge a bin to its parent if the bin is too small
    for (l = idx->n_lvls; l > 0; --l) {
        unsigned start = hts_bin_first(l);
        for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
            bins_t *p, *q;
            if (!kh_exist(bidx, k) || kh_key(bidx, k) >= ((unsigned int)idx->n_bins) || kh_key(bidx, k) < start) continue;
            p = &kh_value(bidx, k);
            if (l < idx->n_lvls && p->n > 1) ks_introsort(_off, p->n, p->list);
            if ((p->list[p->n - 1].v>>16) - (p->list[0].u>>16) < HTS_MIN_MARKER_DIST) {
                khint_t kp;
                kp = kh_get(bin, bidx, hts_bin_parent(kh_key(bidx, k)));
                if (kp == kh_end(bidx)) continue;
                q = &kh_val(bidx, kp);
                if (q->n + p->n > q->m) {
                    q->m = q->n + p->n;
                    kroundup32(q->m);
                    q->list = (hts_pair64_t*)realloc(q->list, q->m * sizeof(hts_pair64_t));
                }
                memcpy(q->list + q->n, p->list, p->n * sizeof(hts_pair64_t));
                q->n += p->n;
                free(p->list);
                kh_del(bin, bidx, k);
            }
        }
    }
    k = kh_get(bin, bidx, 0);
    if (k != kh_end(bidx)) ks_introsort(_off, kh_val(bidx, k).n, kh_val(bidx, k).list);
    // merge adjacent chunks that start from the same BGZF block
    for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
        bins_t *p;
        if (!kh_exist(bidx, k) || kh_key(bidx, k) >= ((unsigned int)idx->n_bins)) continue;
        p = &kh_value(bidx, k);
        for (l = 1, m = 0; l < p->n; ++l) {
            if (p->list[m].v>>16 >= p->list[l].u>>16) {
                if (p->list[m].v < p->list[l].v) p->list[m].v = p->list[l].v;
            } else p->list[++m] = p->list[l];
        }
        p->n = m + 1;
    }
}

void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)
{
    int i;
    if (idx == NULL || idx->z.finished) return; // do not run this function on an empty index or multiple times
    if (idx->z.save_tid >= 0) {
        insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, final_offset);
        insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.off_beg, final_offset);
        insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.n_mapped, idx->z.n_unmapped);
    }
    for (i = 0; i < idx->n; ++i) {
        update_loff(idx, i, (idx->fmt == HTS_FMT_CSI));
        compress_binning(idx, i);
    }
    idx->z.finished = 1;
}

int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped)
{
    int bin;
    int64_t maxpos = (int64_t) 1 << (idx->min_shift + idx->n_lvls * 3);
    if (tid<0) beg = -1, end = 0;
    if (tid >= 0 && (beg > maxpos || end > maxpos)) {
        goto pos_too_big;
    }
    if (tid >= idx->m) { // enlarge the index
        uint32_t new_m = idx->m * 2 > tid + 1 ? idx->m * 2 : tid + 1;
        bidx_t **new_bidx;
        lidx_t *new_lidx;

        new_bidx = (bidx_t**)realloc(idx->bidx, new_m * sizeof(bidx_t*));
        if (!new_bidx) return -1;
        idx->bidx = new_bidx;

        new_lidx = (lidx_t*) realloc(idx->lidx, new_m * sizeof(lidx_t));
        if (!new_lidx) return -1;
        idx->lidx = new_lidx;

        memset(&idx->bidx[idx->m], 0, (new_m - idx->m) * sizeof(bidx_t*));
        memset(&idx->lidx[idx->m], 0, (new_m - idx->m) * sizeof(lidx_t));
        idx->m = new_m;
    }
    if (idx->n < tid + 1) idx->n = tid + 1;
    if (idx->z.finished) return 0;
    if (idx->z.last_tid != tid || (idx->z.last_tid >= 0 && tid < 0)) { // change of chromosome
        if ( tid>=0 && idx->n_no_coor )
        {
            hts_log_error("NO_COOR reads not in a single block at the end %d %d", tid, idx->z.last_tid);
            return -1;
        }
        if (tid>=0 && idx->bidx[tid] != 0)
        {
            hts_log_error("Chromosome blocks not continuous");
            return -1;
        }
        idx->z.last_tid = tid;
        idx->z.last_bin = 0xffffffffu;
    } else if (tid >= 0 && idx->z.last_coor > beg) { // test if positions are out of order
        hts_log_error("Unsorted positions on sequence #%d: %d followed by %d", tid+1, idx->z.last_coor+1, beg+1);
        return -1;
    }
    if ( tid>=0 )
    {
        if (idx->bidx[tid] == 0) idx->bidx[tid] = kh_init(bin);
        if (is_mapped) {
            // shoehorn [-1,0) (VCF POS=0) into the leftmost bottom-level bin
            if (beg < 0)  beg = 0;
            if (end <= 0) end = 1;
            // idx->z.last_off points to the start of the current record
            if (insert_to_l(&idx->lidx[tid], beg, end,
                            idx->z.last_off, idx->min_shift) < 0) return -1;
        }
    }
    else idx->n_no_coor++;
    bin = hts_reg2bin(beg, end, idx->min_shift, idx->n_lvls);
    if ((int)idx->z.last_bin != bin) { // then possibly write the binning index
        if (idx->z.save_bin != 0xffffffffu) { // save_bin==0xffffffffu only happens to the first record
            if (insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin,
                            idx->z.save_off, idx->z.last_off) < 0) return -1;
        }
        if (idx->z.last_bin == 0xffffffffu && idx->z.save_bin != 0xffffffffu) { // change of chr; keep meta information
            idx->z.off_end = idx->z.last_off;
            if (insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx),
                            idx->z.off_beg, idx->z.off_end) < 0) return -1;
            if (insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx),
                            idx->z.n_mapped, idx->z.n_unmapped) < 0) return -1;
            idx->z.n_mapped = idx->z.n_unmapped = 0;
            idx->z.off_beg = idx->z.off_end;
        }
        idx->z.save_off = idx->z.last_off;
        idx->z.save_bin = idx->z.last_bin = bin;
        idx->z.save_tid = tid;
    }
    if (is_mapped) ++idx->z.n_mapped;
    else ++idx->z.n_unmapped;
    idx->z.last_off = offset;
    idx->z.last_coor = beg;
    return 0;

 pos_too_big: {
        int64_t max = end > beg ? end : beg, s = 1 << 14;
        int n_lvls = 0;
        while (max > s) {
            n_lvls++;
            s <<= 3;
        }

        if (idx->fmt == HTS_FMT_CSI) {
            hts_log_error("Region %d..%d cannot be stored in a csi index "
                "with min_shift = %d, n_lvls = %d. Try using "
                "min_shift = 14, n_lvls >= %d",
                beg, end,
                idx->min_shift, idx->n_lvls,
                n_lvls);
        } else {
            hts_log_error("Region %d..%d cannot be stored in a %s index. "
                "Try using a csi index with min_shift = 14, "
                "n_lvls >= %d",
                beg, end, idx_format_name(idx->fmt),
                n_lvls);
        }
        errno = ERANGE;
        return -1;
    }
}

void hts_idx_destroy(hts_idx_t *idx)
{
    khint_t k;
    int i;
    if (idx == 0) return;

    for (i = 0; i < idx->m; ++i) {
        bidx_t *bidx = idx->bidx[i];
        free(idx->lidx[i].offset);
        if (bidx == 0) continue;
        for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
            if (kh_exist(bidx, k))
                free(kh_value(bidx, k).list);
        kh_destroy(bin, bidx);
    }
    free(idx->bidx); free(idx->lidx); free(idx->meta);
    free(idx);
}

// The optimizer eliminates these ed_is_big() calls; still it would be good to
// TODO Determine endianness at configure- or compile-time

static inline ssize_t HTS_RESULT_USED idx_write_int32(BGZF *fp, int32_t x)
{
    if (ed_is_big()) x = ed_swap_4(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline ssize_t HTS_RESULT_USED idx_write_uint32(BGZF *fp, uint32_t x)
{
    if (ed_is_big()) x = ed_swap_4(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline ssize_t HTS_RESULT_USED idx_write_uint64(BGZF *fp, uint64_t x)
{
    if (ed_is_big()) x = ed_swap_8(x);
    return bgzf_write(fp, &x, sizeof x);
}

static int hts_idx_save_core(const hts_idx_t *idx, BGZF *fp, int fmt)
{
    int32_t i, j;

    #define check(ret) if ((ret) < 0) return -1

    check(idx_write_int32(fp, idx->n));
    if (fmt == HTS_FMT_TBI && idx->l_meta)
        check(bgzf_write(fp, idx->meta, idx->l_meta));

    for (i = 0; i < idx->n; ++i) {
        khint_t k;
        bidx_t *bidx = idx->bidx[i];
        lidx_t *lidx = &idx->lidx[i];
        // write binning index
        check(idx_write_int32(fp, bidx? kh_size(bidx) : 0));
        if (bidx)
            for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
                if (kh_exist(bidx, k)) {
                    bins_t *p = &kh_value(bidx, k);
                    check(idx_write_uint32(fp, kh_key(bidx, k)));
                    if (fmt == HTS_FMT_CSI) check(idx_write_uint64(fp, p->loff));
                    //int j;for(j=0;j<p->n;++j)fprintf(stderr,"%d,%llx,%d,%llx:%llx\n",kh_key(bidx,k),kh_val(bidx, k).loff,j,p->list[j].u,p->list[j].v);
                    check(idx_write_int32(fp, p->n));
                    for (j = 0; j < p->n; ++j) {
                        check(idx_write_uint64(fp, p->list[j].u));
                        check(idx_write_uint64(fp, p->list[j].v));
                    }
                }

        // write linear index
        if (fmt != HTS_FMT_CSI) {
            check(idx_write_int32(fp, lidx->n));
            for (j = 0; j < lidx->n; ++j)
                check(idx_write_uint64(fp, lidx->offset[j]));
        }
    }

    check(idx_write_uint64(fp, idx->n_no_coor));
    return 0;
    #undef check
}

int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)
{
    int ret, save;
    char *fnidx = (char*)calloc(1, strlen(fn) + 5);
    if (fnidx == NULL) return -1;

    strcpy(fnidx, fn);
    switch (fmt) {
    case HTS_FMT_BAI: strcat(fnidx, ".bai"); break;
    case HTS_FMT_CSI: strcat(fnidx, ".csi"); break;
    case HTS_FMT_TBI: strcat(fnidx, ".tbi"); break;
    default: abort();
    }

    ret = hts_idx_save_as(idx, fn, fnidx, fmt);
    save = errno;
    free(fnidx);
    errno = save;
    return ret;
}

int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt)
{
    BGZF *fp;

    #define check(ret) if ((ret) < 0) goto fail

    if (fnidx == NULL) return hts_idx_save(idx, fn, fmt);

    fp = bgzf_open(fnidx, (fmt == HTS_FMT_BAI)? "wu" : "w");
    if (fp == NULL) return -1;

    if (fmt == HTS_FMT_CSI) {
        check(bgzf_write(fp, "CSI\1", 4));
        check(idx_write_int32(fp, idx->min_shift));
        check(idx_write_int32(fp, idx->n_lvls));
        check(idx_write_uint32(fp, idx->l_meta));
        if (idx->l_meta) check(bgzf_write(fp, idx->meta, idx->l_meta));
    } else if (fmt == HTS_FMT_TBI) {
        check(bgzf_write(fp, "TBI\1", 4));
    } else if (fmt == HTS_FMT_BAI) {
        check(bgzf_write(fp, "BAI\1", 4));
    } else abort();

    check(hts_idx_save_core(idx, fp, fmt));

    return bgzf_close(fp);
    #undef check

fail:
    bgzf_close(fp);
    return -1;
}


void hts_set_log_level(enum htsLogLevel level)
{
    hts_verbose = level;
}

enum htsLogLevel hts_get_log_level()
{
    return hts_verbose;
}

static char get_severity_tag(enum htsLogLevel severity)
{
    switch (severity) {
    case HTS_LOG_ERROR:
        return 'E';
    case HTS_LOG_WARNING:
        return 'W';
    case HTS_LOG_INFO:
        return 'I';
    case HTS_LOG_DEBUG:
        return 'D';
    case HTS_LOG_TRACE:
        return 'T';
    default:
        break;
    }

    return '*';
}

void hts_log(enum htsLogLevel severity, const char *context, const char *format, ...)
{
    int save_errno = errno;
    if (((int)severity) <= hts_verbose) {
        va_list argptr;

        fprintf(stderr, "[%c::%s] ", get_severity_tag(severity), context);

        va_start(argptr, format);
        vfprintf(stderr, format, argptr);
        va_end(argptr);

        fprintf(stderr, "\n");
    }
    errno = save_errno;
}
