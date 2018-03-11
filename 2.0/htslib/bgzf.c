/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013-2017 Genome Research Ltd

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

// This is a heavily modified copy of commit
// 77f002106602103150db45b06ddec0f022c0b6fb : most code which isn't needed to
// generate bgzf-compressed files with libdeflate has been removed.

// #include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#ifndef _WIN32
#  include <pthread.h>
#endif
#include <sys/types.h>
#include <inttypes.h>

// coming soon...
// #include "../libdeflate/libdeflate.h"

#include "htslib/hts.h"
#include "htslib/bgzf.h"

// hFILE, hopen(), hclose(), hclose_abruptly(), hwrite(), hflush()
#include "htslib/hfile.h"

#ifndef _WIN32
// hts_tpool, hts_tpool_process, hts_tpool_result,
// hts_tpool_next_result_wait(), hts_tpool_result_data(),
// hts_tpool_delete_result(), hts_tpool_process_destroy(), hts_tpool_size(),
// hts_tpool_process_init(), hts_tpool_process_ref_incr(), hts_tpool_init(),
// hts_tpool_wake_dispatch(), hts_tpool_destroy(), hts_tpool_dispatch(),
// hts_tpool_process_flush()
#  include "htslib/thread_pool.h"
#endif

#include "htslib/hts_endian.h"  // u16_to_le()

// pool_alloc_t, pool_create(), pool_destroy(), pool_alloc(), pool_free()
#include "cram/pooled_alloc.h"

#ifndef _WIN32
#  define BGZF_MT
#endif

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


/* BGZF/GZIP header (speciallized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C

  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.

*/
static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

#ifdef BGZF_MT

typedef struct bgzf_job {
    BGZF *fp;
    unsigned char comp_data[BGZF_MAX_BLOCK_SIZE];
    size_t comp_len;
    unsigned char uncomp_data[BGZF_MAX_BLOCK_SIZE];
    size_t uncomp_len;
    int errcode;
    int64_t block_address;
    int hit_eof;
} bgzf_job;

enum mtaux_cmd {
    NONE = 0,
    SEEK,
    HAS_EOF,
    CLOSE,
};

typedef struct bgzf_mtaux_t {
    // Memory pool for bgzf_job structs, to avoid many malloc/free
    pool_alloc_t *job_pool;
    bgzf_job *curr_job;

    // Thread pool
    int n_threads;
    int own_pool;
    hts_tpool *pool;

    // Output queue holding completed bgzf_jobs
    hts_tpool_process *out_queue;

    // I/O thread.
    pthread_t io_task;
    pthread_mutex_t job_pool_m;
    int jobs_pending; // number of jobs waiting
    int flush_pending;
    void *free_block;
    int hit_eof;  // r/w entirely within main thread

    // Message passing to the reader thread; eg seek requests
    int errcode;
    uint64_t block_address;
    int eof;
    pthread_mutex_t command_m; // Set whenever fp is being updated
    pthread_cond_t command_c;
    enum mtaux_cmd command;
} mtaux_t;
#endif

typedef struct
{
    uint64_t uaddr;  // offset w.r.t. uncompressed data
    uint64_t caddr;  // offset w.r.t. compressed data
}
bgzidx1_t;

struct __bgzidx_t
{
    int noffs, moffs;       // the size of the index, n:used, m:allocated
    bgzidx1_t *offs;        // offsets
    uint64_t ublock_addr;   // offset of the current block (uncompressed data)
};

void bgzf_index_destroy(BGZF *fp);
int bgzf_index_add_block(BGZF *fp);
#ifdef BGZF_MT
static void mt_destroy(mtaux_t *mt);
#endif

static inline void packInt16(uint8_t *buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

static const char *bgzf_zerr(int errnum, z_stream *zs)
{
    static char buffer[32];

    /* Return zs->msg if available.
       zlib doesn't set this very reliably.  Looking at the source suggests
       that it may get set to a useful message for deflateInit2, inflateInit2
       and inflate when it returns Z_DATA_ERROR. For inflate with other
       return codes, deflate, deflateEnd and inflateEnd it doesn't appear
       to be useful.  For the likely non-useful cases, the caller should
       pass NULL into zs. */

    if (zs && zs->msg) return zs->msg;

    // gzerror OF((gzFile file, int *errnum)
    switch (errnum) {
    case Z_ERRNO:
        return strerror(errno);
    case Z_STREAM_ERROR:
        return "invalid parameter/compression level, or inconsistent stream state";
    case Z_DATA_ERROR:
        return "invalid or incomplete IO";
    case Z_MEM_ERROR:
        return "out of memory";
    case Z_BUF_ERROR:
        return "progress temporarily not possible, or in() / out() returned an error";
    case Z_VERSION_ERROR:
        return "zlib version mismatch";
    case Z_OK: // 0: maybe gzgets error Z_NULL
    default:
        snprintf(buffer, sizeof(buffer), "[%d] unknown", errnum);
        return buffer;  // FIXME: Not thread-safe.
    }
}

// get the compress level from the mode string: compress_level==-1 for the default level, -2 plain uncompressed
static int mode2level(const char *mode)
{
    int i, compress_level = -1;
    for (i = 0; mode[i]; ++i)
        if (mode[i] >= '0' && mode[i] <= '9') break;
    if (mode[i]) compress_level = (int)mode[i] - '0';
    if (strchr(mode, 'u')) compress_level = -2;
    return compress_level;
}
static BGZF *bgzf_write_init(const char *mode)
{
    BGZF *fp;
    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) goto mem_fail;
    fp->is_write = 1;
    int compress_level = mode2level(mode);
    if ( compress_level==-2 )
    {
        fp->is_compressed = 0;
        return fp;
    }
    fp->is_compressed = 1;

    fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
    if (fp->uncompressed_block == NULL) goto mem_fail;
    fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;

    fp->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1
    if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;
    /*
    if ( strchr(mode,'g') )
    {
        // gzip output
        fp->is_gzip = 1;
        fp->gz_stream = (z_stream*)calloc(1,sizeof(z_stream));
        if (fp->gz_stream == NULL) goto mem_fail;
        fp->gz_stream->zalloc = NULL;
        fp->gz_stream->zfree  = NULL;
        fp->gz_stream->msg    = NULL;

        int ret = deflateInit2(fp->gz_stream, fp->compress_level, Z_DEFLATED, 15|16, 8, Z_DEFAULT_STRATEGY);
        if (ret!=Z_OK) {
            hts_log_error("Call to deflateInit2 failed: %s", bgzf_zerr(ret, fp->gz_stream));
            goto fail;
        }
    }
    */
    return fp;

mem_fail:
    hts_log_error("%s", strerror(errno));

    /*
fail:
    if (fp != NULL) {
        free(fp->uncompressed_block);
        free(fp->gz_stream);
        free(fp);
    }
    */
    return NULL;
}

BGZF *bgzf_open(const char *path, const char *mode)
{
    BGZF *fp = 0;
    assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
    /*
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hopen(path, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; }
        fp->fp = fpr;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
    */
        hFILE *fpw;
        if ((fpw = hopen(path, mode)) == 0) return 0;
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
        fp->fp = fpw;
        /*
    }
    else { errno = EINVAL; return 0; }
        */

    fp->is_be = ed_is_big();
    return fp;
}

/*
int bgzf_compress(void *_dst, size_t *dlen, const void *src, size_t slen, int level)
{
    if (slen == 0) {
        // EOF block
        if (*dlen < 28) return -1;
        memcpy(_dst, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28);
        *dlen = 28;
        return 0;
    }

    uint8_t *dst = (uint8_t*)_dst;

    if (level == 0) {
        // Uncompressed data
        if (*dlen < slen+5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH) return -1;
        dst[BLOCK_HEADER_LENGTH] = 1; // BFINAL=1, BTYPE=00; see RFC1951
        u16_to_le(slen,  &dst[BLOCK_HEADER_LENGTH+1]); // length
        u16_to_le(~slen, &dst[BLOCK_HEADER_LENGTH+3]); // ones-complement length
        memcpy(dst + BLOCK_HEADER_LENGTH+5, src, slen);
        *dlen = slen+5 + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

    } else {
        level = level > 0 ? level : 6; // libdeflate doesn't honour -1 as default
        // NB levels go up to 12 here.
        struct libdeflate_compressor *z = libdeflate_alloc_compressor(level);
        if (!z) return -1;

        // Raw deflate
        size_t clen =
            libdeflate_deflate_compress(z, src, slen,
                                        dst + BLOCK_HEADER_LENGTH,
                                        *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);

        if (clen <= 0) {
            hts_log_error("Call to libdeflate_deflate_compress failed");
            libdeflate_free_compressor(z);
            return -1;
        }

        *dlen = clen + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

        libdeflate_free_compressor(z);
    }

    // write the header
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
    packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes

    // write the footer
    uint32_t crc = libdeflate_crc32(0, src, slen);
    packInt32((uint8_t*)&dst[*dlen - 8], crc);
    packInt32((uint8_t*)&dst[*dlen - 4], slen);
    return 0;
}
*/
int bgzf_compress(void *_dst, size_t *dlen, const void *src, size_t slen, int level)
{
    uint32_t crc;
    z_stream zs;
    uint8_t *dst = (uint8_t*)_dst;

    // compress the body
    zs.zalloc = NULL; zs.zfree = NULL;
    zs.msg = NULL;
    zs.next_in  = (Bytef*)src;
    zs.avail_in = slen;
    zs.next_out = dst + BLOCK_HEADER_LENGTH;
    zs.avail_out = *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
    int ret = deflateInit2(&zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY); // -15 to disable zlib header/footer
    if (ret!=Z_OK) {
        hts_log_error("Call to deflateInit2 failed: %s", bgzf_zerr(ret, &zs));
        return -1;
    }
    if ((ret = deflate(&zs, Z_FINISH)) != Z_STREAM_END) {
        hts_log_error("Deflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? &zs : NULL));
        return -1;
    }
    if ((ret = deflateEnd(&zs)) != Z_OK) {
        hts_log_error("Call to deflateEnd failed: %s", bgzf_zerr(ret, NULL));
        return -1;
    }
    *dlen = zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    // write the header
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
    packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes
    // write the footer
    crc = crc32(crc32(0L, NULL, 0L), (Bytef*)src, slen);
    packInt32((uint8_t*)&dst[*dlen - 8], crc);
    packInt32((uint8_t*)&dst[*dlen - 4], slen);
    return 0;
}


// Deflate the block in fp->uncompressed_block into fp->compressed_block. Also adds an extra field that stores the compressed block length.
static int deflate_block(BGZF *fp, int block_length)
{
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    ret = bgzf_compress(fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level);

    if ( ret != 0 )
    {
        hts_log_debug("Compression error %d", ret);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    fp->block_offset = 0;
    return comp_size;
}

#ifdef BGZF_MT

void *bgzf_encode_func(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;

    j->comp_len = BGZF_MAX_BLOCK_SIZE;
    int ret = bgzf_compress(j->comp_data, &j->comp_len,
                            j->uncomp_data, j->uncomp_len,
                            j->fp->compress_level);
    if (ret != 0)
        j->errcode |= BGZF_ERR_ZLIB;

    return arg;
}

/*
 * Takes compressed blocks off the results queue and calls hwrite to
 * punt them to the output stream.
 *
 * Returns NULL when no more are left, or -1 on error
 */
static void *bgzf_mt_writer(void *vp) {
    BGZF *fp = (BGZF *)vp;
    mtaux_t *mt = fp->mt;
    hts_tpool_result *r;

    // Iterates until result queue is shutdown, where it returns NULL.
    while ((r = hts_tpool_next_result_wait(mt->out_queue))) {
        bgzf_job *j = (bgzf_job *)hts_tpool_result_data(r);
        assert(j);

        if (hwrite(fp->fp, j->comp_data, j->comp_len) != (ssize_t)j->comp_len) {
            fp->errcode |= BGZF_ERR_IO;
            goto err;
        }

        /*
         * Periodically call hflush (which calls fsync when on a file).
         * This avoids the fsync being done at the bgzf_close stage,
         * which can sometimes cause signficant delays.  As this is in
         * a separate thread, spreading the sync delays throughout the
         * program execution seems better.
         * Frequency of 1/512 has been chosen by experimentation
         * across local XFS, NFS and Lustre tests.
         */
        if (++mt->flush_pending % 512 == 0)
            if (hflush(fp->fp) != 0)
                goto err;


        hts_tpool_delete_result(r, 0);

        // Also updated by main thread
        pthread_mutex_lock(&mt->job_pool_m);
        pool_free(mt->job_pool, j);
        mt->jobs_pending--;
        pthread_mutex_unlock(&mt->job_pool_m);
    }

    if (hflush(fp->fp) != 0)
        goto err;

    hts_tpool_process_destroy(mt->out_queue);

    return NULL;

 err:
    hts_tpool_process_destroy(mt->out_queue);
    return (void *)-1;
}

int bgzf_thread_pool(BGZF *fp, hts_tpool *pool, int qsize) {
    // No gain from multi-threading when not compressed
    if (!fp->is_compressed)
        return 0;

    mtaux_t *mt;
    mt = (mtaux_t*)calloc(1, sizeof(mtaux_t));
    if (!mt) return -1;
    fp->mt = mt;

    mt->pool = pool;
    mt->n_threads = hts_tpool_size(pool);
    if (!qsize)
        qsize = mt->n_threads*2;
    if (!(mt->out_queue = hts_tpool_process_init(mt->pool, qsize, 0))) {
        free(mt);
        return -1;
    }
    hts_tpool_process_ref_incr(mt->out_queue);

    mt->job_pool = pool_create(sizeof(bgzf_job));

    pthread_mutex_init(&mt->job_pool_m, NULL);
    pthread_mutex_init(&mt->command_m, NULL);
    pthread_cond_init(&mt->command_c, NULL);
    mt->flush_pending = 0;
    mt->jobs_pending = 0;
    mt->free_block = fp->uncompressed_block; // currently in-use block
    pthread_create(&mt->io_task, NULL,
                   // fp->is_write ? bgzf_mt_writer : bgzf_mt_reader, fp);
                   bgzf_mt_writer, fp);
    return 0;
}

int bgzf_mt(BGZF *fp, int n_threads, __attribute__((unused)) int n_sub_blks)
{
    // No gain from multi-threading when not compressed
    if (!fp->is_compressed || fp->is_gzip)
        return 0;

    if (n_threads < 1) return -1;
    hts_tpool *p = hts_tpool_init(n_threads);
    if (!p)
        return -1;

    if (bgzf_thread_pool(fp, p, 0) != 0) {
        hts_tpool_destroy(p);
        return -1;
    }

    fp->mt->own_pool = 1;

    return 0;
}

static void mt_destroy(mtaux_t *mt)
{
    pthread_mutex_lock(&mt->command_m);
    mt->command = CLOSE;
    pthread_cond_signal(&mt->command_c);
    hts_tpool_wake_dispatch(mt->out_queue); // unstick the reader
    pthread_mutex_unlock(&mt->command_m);

    // Destroying the queue first forces the writer to exit.
    hts_tpool_process_destroy(mt->out_queue);
    pthread_join(mt->io_task, NULL);

    pthread_mutex_destroy(&mt->job_pool_m);
    pthread_mutex_destroy(&mt->command_m);
    pthread_cond_destroy(&mt->command_c);
    if (mt->curr_job)
        pool_free(mt->job_pool, mt->curr_job);

    if (mt->own_pool)
        hts_tpool_destroy(mt->pool);

    pool_destroy(mt->job_pool);

    free(mt);
    fflush(stderr);
}

static int mt_queue(BGZF *fp)
{
    mtaux_t *mt = fp->mt;

    // Also updated by writer thread
    pthread_mutex_lock(&mt->job_pool_m);
    bgzf_job *j = pool_alloc(mt->job_pool);
    mt->jobs_pending++;
    pthread_mutex_unlock(&mt->job_pool_m);

    j->fp = fp;
    j->errcode = 0;
    j->uncomp_len  = fp->block_offset;
    memcpy(j->uncomp_data, fp->uncompressed_block, j->uncomp_len);

    // Need non-block vers & job_pending?
    hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_encode_func, j);

    fp->block_offset = 0;
    return 0;
}

static int mt_flush_queue(BGZF *fp)
{
    mtaux_t *mt = fp->mt;

    // Drain the encoder jobs.
    // We cannot use hts_tpool_flush here as it can cause deadlock if
    // the queue is full up of decoder tasks.  The best solution would
    // be to have one input queue per type of job, but we don't right now.
    //hts_tpool_flush(mt->pool);
    pthread_mutex_lock(&mt->job_pool_m);
    while (mt->jobs_pending != 0) {
        pthread_mutex_unlock(&mt->job_pool_m);
        usleep(10000); // FIXME: replace by condition variable
        pthread_mutex_lock(&mt->job_pool_m);
    }
    pthread_mutex_unlock(&mt->job_pool_m);

    // Wait on bgzf_mt_writer to drain the queue
    if (hts_tpool_process_flush(mt->out_queue) != 0)
        return -1;

    return (fp->errcode == 0)? 0 : -1;
}

static int lazy_flush(BGZF *fp)
{
    if (fp->mt)
        return fp->block_offset ? mt_queue(fp) : 0;
    else
        return bgzf_flush(fp);
}

#else  // ~ #ifdef BGZF_MT

int bgzf_mt(__attribute__((unused)) BGZF *fp, __attribute__((unused)) int n_threads, __attribute__((unused))  int n_sub_blks)
{
    return 0;
}

static inline int lazy_flush(BGZF *fp)
{
    return bgzf_flush(fp);
}

#endif // ~ #ifdef BGZF_MT

int bgzf_flush(BGZF *fp)
{
    if (!fp->is_write) return 0;
#ifdef BGZF_MT
    if (fp->mt) {
        int ret = 0;
        if (fp->block_offset) ret = mt_queue(fp);
        return ret ? ret : mt_flush_queue(fp);
    }
#endif
    while (fp->block_offset > 0) {
        int block_length;
        if ( fp->idx_build_otf )
        {
            bgzf_index_add_block(fp);
            fp->idx->ublock_addr += fp->block_offset;
        }
        block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, fp->compressed_block, block_length) != block_length) {
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }
        fp->block_address += block_length;
    }
    return 0;
}

int bgzf_flush_try(BGZF *fp, ssize_t size)
{
    if (fp->block_offset + size > BGZF_BLOCK_SIZE) return lazy_flush(fp);
    return 0;
}

ssize_t bgzf_write(BGZF *fp, const void *data, size_t length)
{
    if ( !fp->is_compressed )
        return hwrite(fp->fp, data, length);

    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
    assert(fp->is_write);
    while (remaining > 0) {
        uint8_t* buffer = (uint8_t*)fp->uncompressed_block;
        int copy_length = BGZF_BLOCK_SIZE - fp->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (fp->block_offset == BGZF_BLOCK_SIZE) {
            if (lazy_flush(fp) != 0) return -1;
        }
    }
    return length - remaining;
}

ssize_t bgzf_block_write(BGZF *fp, const void *data, size_t length)
{
    if ( !fp->is_compressed )
        return hwrite(fp->fp, data, length);
    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
    assert(fp->is_write);
    uint64_t current_block; //keep track of current block
    uint64_t ublock_size; // amount of uncompressed data to be fed into next block
    while (remaining > 0) {
        current_block = fp->idx->moffs - fp->idx->noffs;
        ublock_size = current_block + 1 < (uint64_t)fp->idx->moffs ? fp->idx->offs[current_block+1].uaddr-fp->idx->offs[current_block].uaddr : BGZF_MAX_BLOCK_SIZE;
        uint8_t* buffer = (uint8_t*)fp->uncompressed_block;
        int copy_length = ublock_size - fp->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if ((uint64_t)fp->block_offset == ublock_size) {
            if (lazy_flush(fp) != 0) return -1;
            if (fp->idx->noffs > 0)
                fp->idx->noffs--;  // decrement noffs to track the blocks
        }
    }
    return length - remaining;
}


ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length)
{
    ssize_t ret = hwrite(fp->fp, data, length);
    if (ret < 0) fp->errcode |= BGZF_ERR_IO;
    return ret;
}

int bgzf_close(BGZF* fp)
{
    int ret, block_length;
    if (fp == 0) return -1;
    if (fp->is_write && fp->is_compressed) {
        if (bgzf_flush(fp) != 0) return -1;
        fp->compress_level = -1;
        block_length = deflate_block(fp, 0); // write an empty block
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, fp->compressed_block, block_length) < 0
            || hflush(fp->fp) != 0) {
            hts_log_error("File write failed");
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
    }
#ifdef BGZF_MT
    if (fp->mt) {
        if (!fp->mt->free_block)
            fp->uncompressed_block = NULL;
        mt_destroy(fp->mt);
    }
#endif
    if ( fp->is_gzip )
    {
        if (fp->gz_stream == NULL) ret = Z_OK;
        else if (!fp->is_write) ret = inflateEnd(fp->gz_stream);
        else ret = deflateEnd(fp->gz_stream);
        if (ret != Z_OK) {
            hts_log_error("Call to inflateEnd/deflateEnd failed: %s", bgzf_zerr(ret, NULL));
        }
        free(fp->gz_stream);
    }
    ret = hclose(fp->fp);
    if (ret != 0) return -1;
    bgzf_index_destroy(fp);
    free(fp->uncompressed_block);
    // free_cache(fp);
    free(fp);
    return 0;
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

void bgzf_index_destroy(BGZF *fp)
{
    if ( !fp->idx ) return;
    free(fp->idx->offs);
    free(fp->idx);
    fp->idx = NULL;
    fp->idx_build_otf = 0;
}

int bgzf_index_build_init(BGZF *fp)
{
    bgzf_index_destroy(fp);
    fp->idx = (bgzidx_t*) calloc(1,sizeof(bgzidx_t));
    if ( !fp->idx ) return -1;
    fp->idx_build_otf = 1;  // build index on the fly
    return 0;
}

int bgzf_index_add_block(BGZF *fp)
{
    fp->idx->noffs++;
    if ( fp->idx->noffs > fp->idx->moffs )
    {
        fp->idx->moffs = fp->idx->noffs;
        kroundup32(fp->idx->moffs);
        fp->idx->offs = (bgzidx1_t*) realloc(fp->idx->offs, fp->idx->moffs*sizeof(bgzidx1_t));
        if ( !fp->idx->offs ) return -1;
    }
    fp->idx->offs[ fp->idx->noffs-1 ].uaddr = fp->idx->ublock_addr;
    fp->idx->offs[ fp->idx->noffs-1 ].caddr = fp->block_address;
    return 0;
}

static inline int hwrite_uint64(uint64_t x, hFILE *f)
{
    if (ed_is_big()) x = ed_swap_8(x);
    if (hwrite(f, &x, sizeof(x)) != sizeof(x)) return -1;
    return 0;
}

static char * get_name_suffix(const char *bname, const char *suffix)
{
    size_t len = strlen(bname) + strlen(suffix) + 1;
    char *buff = malloc(len);
    if (!buff) return NULL;
    snprintf(buff, len, "%s%s", bname, suffix);
    return buff;
}

int bgzf_index_dump_hfile(BGZF *fp, struct hFILE *idx, const char *name)
{
    // Note that the index contains one extra record when indexing files opened
    // for reading. The terminating record is not present when opened for writing.
    // This is not a bug.

    int i;

    if (!fp->idx) {
        hts_log_error("Called for BGZF handle with no index");
        errno = EINVAL;
        return -1;
    }

    if (bgzf_flush(fp) != 0) return -1;

    if (hwrite_uint64(fp->idx->noffs - 1, idx) < 0) goto fail;
    for (i=1; i<fp->idx->noffs; i++)
    {
        if (hwrite_uint64(fp->idx->offs[i].caddr, idx) < 0) goto fail;
        if (hwrite_uint64(fp->idx->offs[i].uaddr, idx) < 0) goto fail;
    }
    return 0;

 fail:
    hts_log_error("Error writing to %s : %s", name ? name : "index", strerror(errno));
    return -1;
}

int bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix)
{
    const char *name = bname, *msg = NULL;
    char *tmp = NULL;
    hFILE *idx = NULL;

    if (!fp->idx) {
        hts_log_error("Called for BGZF handle with no index");
        errno = EINVAL;
        return -1;
    }

    if ( suffix )
    {
        tmp = get_name_suffix(bname, suffix);
        if ( !tmp ) return -1;
        name = tmp;
    }

    idx = hopen(name, "wb");
    if ( !idx ) {
        msg = "Error opening";
        goto fail;
    }

    if (bgzf_index_dump_hfile(fp, idx, name) != 0) goto fail;

    if (hclose(idx) < 0)
    {
        idx = NULL;
        msg = "Error on closing";
        goto fail;
    }

    free(tmp);
    return 0;

 fail:
    if (msg != NULL) {
        hts_log_error("%s %s : %s", msg, name, strerror(errno));
    }
    if (idx) hclose_abruptly(idx);
    free(tmp);
    return -1;
}
