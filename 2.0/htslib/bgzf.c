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

// This is a modified copy of commit
// 53241915fa8b6a3e807f2ac85d8701f27a7fe528 ; some code which isn't needed to
// read or generate bgzf-compressed files with libdeflate has been removed.

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

#include "../libdeflate/libdeflate.h"

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#ifndef _WIN32
#  include "htslib/thread_pool.h"
#endif
#include "htslib/hts_endian.h"
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

#ifdef BGZF_MT
static void mt_destroy(mtaux_t *mt);
#endif

static inline void packInt16(uint8_t *buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

static BGZF *bgzf_read_init(hFILE *hfpr)
{
    BGZF *fp;
    uint8_t magic[18];
    ssize_t n = hpeek(hfpr, magic, 18);
    if (n < 0) return NULL;

    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) return NULL;

    fp->is_write = 0;
    fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
    if (fp->uncompressed_block == NULL) { free(fp); return NULL; }
    fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;
    fp->is_compressed = (n==18 && magic[0]==0x1f && magic[1]==0x8b);
    // fp->is_gzip = ( !fp->is_compressed || ((magic[3]&4) && memcmp(&magic[12], "BC\2\0",4)==0) ) ? 0 : 1;
    if ( (!fp->is_compressed) || (!((magic[3]&4) && memcmp(&magic[12], "BC\2\0",4)==0))) {
      // Use zlibWrapper to handle non-BGZF cases for now.
      return NULL;
    }
    return fp;
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

    fp->compress_level = compress_level < 0? -1 : compress_level; // Z_DEFAULT_COMPRESSION==-1
    if (fp->compress_level > 9) fp->compress_level = -1;
    return fp;

mem_fail:
    hts_log_error("%s", strerror(errno));

    if (fp != NULL) {
        free(fp->uncompressed_block);
        free(fp);
    }
    return NULL;
}

BGZF *bgzf_open(const char *path, const char *mode)
{
    BGZF *fp = 0;
    // assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hopen(path, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; }
        fp->fp = fpr;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        hFILE *fpw;
        if ((fpw = hopen(path, mode)) == 0) return 0;
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
        fp->fp = fpw;
    }
    else { errno = EINVAL; return 0; }

    fp->is_be = ed_is_big();
    return fp;
}

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

static int bgzf_uncompress(uint8_t *dst, size_t *dlen, const uint8_t *src, size_t slen) {
    struct libdeflate_decompressor *z = libdeflate_alloc_decompressor();
    if (!z) {
        hts_log_error("Call to libdeflate_alloc_decompressor failed");
        return -1;
    }

    int ret = libdeflate_deflate_decompress(z, src, slen, dst, *dlen, dlen);
    libdeflate_free_decompressor(z);

    if (ret != LIBDEFLATE_SUCCESS) {
        hts_log_error("Inflate operation failed: %d", ret);
        return -1;
    }

    return 0;
}

// Inflate the block in fp->compressed_block into fp->uncompressed_block
static int inflate_block(BGZF* fp, int block_length)
{
    size_t dlen = BGZF_MAX_BLOCK_SIZE;
    int ret = bgzf_uncompress(fp->uncompressed_block, &dlen,
                              (const uint8_t*)fp->compressed_block + 18, block_length - 18);
    if (ret < 0) {
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }

    // Check CRC of uncompressed block matches the gzip header.
    // NB: we may wish to switch out the zlib crc32 for something more performant.
    // See PR#361 and issue#467
    uint32_t c1 = libdeflate_crc32(0L, (unsigned char *)fp->uncompressed_block, dlen);
    uint32_t c2 = le_to_u32((uint8_t *)fp->compressed_block + block_length-8);
    if (c1 != c2) {
        fp->errcode |= BGZF_ERR_CRC;
        return -1;
    }

    return dlen;
}

// Returns: 0 on success (BGZF header); -1 on non-BGZF GZIP header; -2 on error
static int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}

/*
 * Absolute htell in this compressed file.
 *
 * Do not confuse with the external bgzf_tell macro which returns the virtual
 * offset.
 */
static off_t bgzf_htell(BGZF *fp) {
#ifdef BGZF_MT
    if (fp->mt) {
        pthread_mutex_lock(&fp->mt->job_pool_m);
        off_t pos = fp->block_address + fp->block_clength;
        pthread_mutex_unlock(&fp->mt->job_pool_m);
        return pos;
    } else {
#endif
        return htell(fp->fp);
#ifdef BGZF_MT
    }
#endif
}

int bgzf_read_block(BGZF *fp)
{
#ifdef BGZF_MT
    hts_tpool_result *r;

    if (fp->mt) {
    again:
        if (fp->mt->hit_eof) {
            // Further reading at EOF will always return 0
            fp->block_length = 0;
            return 0;
        }
        r = hts_tpool_next_result_wait(fp->mt->out_queue);
        bgzf_job *j = r ? (bgzf_job *)hts_tpool_result_data(r) : NULL;

        if (!j || j->errcode == BGZF_ERR_MT) {
            if (!fp->mt->free_block) {
                fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
                if (fp->uncompressed_block == NULL) return -1;
                fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;
            } // else it's already allocated with malloc, maybe even in-use.
            mt_destroy(fp->mt);
            fp->mt = NULL;
            hts_tpool_delete_result(r, 0);

            goto single_threaded;
        }

        if (j->errcode) {
            fp->errcode = j->errcode;
            return -1;
        }

        if (j->hit_eof) {
            if (!fp->last_block_eof && !fp->no_eof_block) {
                fp->no_eof_block = 1;
                hts_log_warning("EOF marker is absent. The input is probably truncated");
            }
            fp->mt->hit_eof = 1;
        }

        // Zero length blocks in the middle of a file are (wrongly)
        // considered as EOF by many callers.  We work around this by
        // trying again to see if we hit a genuine EOF.
        if (!j->hit_eof && j->uncomp_len == 0) {
            fp->last_block_eof = 1;
            hts_tpool_delete_result(r, 0);
            goto again;
        }

        // block_length=0 and block_offset set by bgzf_seek.
        if (fp->block_length != 0) fp->block_offset = 0;
        if (!j->hit_eof) fp->block_address = j->block_address;
        fp->block_clength = j->comp_len;
        fp->block_length = j->uncomp_len;
        // bgzf_read() can change fp->block_length
        fp->last_block_eof = (fp->block_length == 0);

        // Steal the data block as it's quicker than a memcpy.
        // We just need to make sure we delay the pool free.
        if (fp->mt->curr_job) {
            pthread_mutex_lock(&fp->mt->job_pool_m);
            pool_free(fp->mt->job_pool, fp->mt->curr_job);
            pthread_mutex_unlock(&fp->mt->job_pool_m);
        }
        fp->uncompressed_block = j->uncomp_data;
        fp->mt->curr_job = j;
        if (fp->mt->free_block) {
            free(fp->mt->free_block); // clear up last non-mt block
            fp->mt->free_block = NULL;
        }

        hts_tpool_delete_result(r, 0);
        return 0;
    }
#endif  // BGZF_MT

    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size, block_length, remaining;

#ifdef BGZF_MT
 single_threaded:
#endif
    size = 0;

    // Reading compressed file
    int64_t block_address;
    block_address = bgzf_htell(fp);
    // if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;

    // loop to skip empty bgzf blocks
    while (1)
    {
        count = hread(fp->fp, header, sizeof(header));
        if (count == 0) { // no data read
            if (!fp->last_block_eof && !fp->no_eof_block) {
                fp->no_eof_block = 1;
                hts_log_warning("EOF marker is absent. The input is probably truncated");
            }
            fp->block_length = 0;
            return 0;
        }
        int ret;
        if ( count != sizeof(header) || (ret=check_header(header))==-2 )
        {
            fp->errcode |= BGZF_ERR_HEADER;
            return -1;
        }
        size = count;
        block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
        if (block_length < BLOCK_HEADER_LENGTH)
        {
            fp->errcode |= BGZF_ERR_HEADER;
            return -1;
        }
        compressed_block = (uint8_t*)fp->compressed_block;
        memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
        remaining = block_length - BLOCK_HEADER_LENGTH;
        count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
        if (count != remaining) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        size += count;
        if ((count = inflate_block(fp, block_length)) < 0) {
            hts_log_debug("Inflate block operation failed");
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->last_block_eof = (count == 0);
        if ( count ) break;     // otherwise an empty bgzf block
    }
    if (fp->block_length != 0) fp->block_offset = 0; // Do not reset offset if this read follows a seek.
    fp->block_address = block_address;
    fp->block_length = count;
    // cache_block(fp, size);
    return 0;
}

ssize_t bgzf_read(BGZF *fp, void *data, size_t length)
{
    ssize_t bytes_read = 0;
    uint8_t *output = (uint8_t*)data;
    if (length <= 0) return 0;
    assert(fp->is_write == 0);
    while ((size_t)bytes_read < length) {
        int copy_length, available = fp->block_length - fp->block_offset;
        uint8_t *buffer;
        if (available <= 0) {
            int ret = bgzf_read_block(fp);
            if (ret != 0) {
                hts_log_error("Read block operation failed with error %d", ret);
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            available = fp->block_length - fp->block_offset;
            if (available <= 0) break;
        }
        copy_length = length - bytes_read < (unsigned int)available? length - bytes_read : (unsigned int)available;
        buffer = (uint8_t*)fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;

        // For raw gzip streams this avoids short reads.
        if (fp->block_offset == fp->block_length) {
            fp->block_address = bgzf_htell(fp);
            fp->block_offset = fp->block_length = 0;
        }
    }

    fp->uncompressed_address += bytes_read;

    return bytes_read;
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

// Our input block has already been decoded by bgzf_mt_read_block().
// We need to split that into a fetch block (compressed) and make this
// do the actual decompression step.
void *bgzf_decode_func(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;

    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
    int ret = bgzf_uncompress(j->uncomp_data, &j->uncomp_len,
                              j->comp_data+18, j->comp_len-18);
    if (ret != 0)
        j->errcode |= BGZF_ERR_ZLIB;

    return arg;
}

/*
 * Nul function so we can dispatch a job with the correct serial
 * to mark failure or to indicate an empty read (EOF).
 */
void *bgzf_nul_func(void *arg) { return arg; }

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


/*
 * Reads a compressed block of data using hread and dispatches it to
 * the thread pool for decompression.  This is the analogue of the old
 * non-threaded bgzf_read_block() function, but without modifying fp
 * in any way (except for the read offset).  All output goes via the
 * supplied bgzf_job struct.
 *
 * Returns NULL when no more are left, or -1 on error
 */
int bgzf_mt_read_block(BGZF *fp, bgzf_job *j)
{
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size = 0, block_length, remaining;

    // NOTE: Guaranteed to be compressed as we block multi-threading in
    // uncompressed mode.  However it may be gzip compression instead
    // of bgzf.

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);

    // if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hpeek(fp->fp, header, sizeof(header));
    if (count == 0) // no data read
        return -1;
    int ret;
    if ( count != sizeof(header) || (ret=check_header(header))==-2 )
    {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if (ret == -1) {
        j->errcode |= BGZF_ERR_MT;
        return -1;
    }

    count = hread(fp->fp, header, sizeof(header));
    if (count != sizeof(header)) // no data read
        return -1;

    size = count;
    block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
    if (block_length < BLOCK_HEADER_LENGTH) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    compressed_block = (uint8_t*)j->comp_data;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        j->errcode |= BGZF_ERR_IO;
        return -1;
    }
    size += count;
    j->comp_len = block_length;
    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
    j->block_address = block_address;
    j->fp = fp;
    j->errcode = 0;

    return 0;
}


static int bgzf_check_EOF_common(BGZF *fp)
{
    uint8_t buf[28];
    off_t offset = htell(fp->fp);
    if (hseek(fp->fp, -28, SEEK_END) < 0) {
        if (errno == ESPIPE) { hclearerr(fp->fp); return 2; }
#ifdef _WIN32
        if (errno == EINVAL) { hclearerr(fp->fp); return 2; }
#endif
        else return -1;
    }
    if ( hread(fp->fp, buf, 28) != 28 ) return -1;
    if ( hseek(fp->fp, offset, SEEK_SET) < 0 ) return -1;
    return (memcmp("\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", buf, 28) == 0)? 1 : 0;
}

/*
 * Checks EOF from the reader thread.
 */
static void bgzf_mt_eof(BGZF *fp) {
    mtaux_t *mt = fp->mt;

    pthread_mutex_lock(&mt->job_pool_m);
    mt->eof = bgzf_check_EOF_common(fp);
    pthread_mutex_unlock(&mt->job_pool_m);
    pthread_cond_signal(&mt->command_c);
}


/*
 * Performs the seek (called by reader thread).
 *
 * This simply drains the entire queue, throwing away blocks, seeks,
 * and starts it up again.  Brute force, but maybe sufficient.
 */
static void bgzf_mt_seek(BGZF *fp) {
    mtaux_t *mt = fp->mt;

    hts_tpool_process_reset(mt->out_queue, 0);
    pthread_mutex_lock(&mt->job_pool_m);
    mt->command = NONE;
    mt->errcode = 0;

    if (hseek(fp->fp, mt->block_address, SEEK_SET) < 0)
        mt->errcode = BGZF_ERR_IO;

    pthread_mutex_unlock(&mt->job_pool_m);
    pthread_cond_signal(&mt->command_c);
}

static void *bgzf_mt_reader(void *vp) {
    BGZF *fp = (BGZF *)vp;
    mtaux_t *mt = fp->mt;

restart:
    pthread_mutex_lock(&mt->job_pool_m);
    bgzf_job *j = pool_alloc(mt->job_pool);
    pthread_mutex_unlock(&mt->job_pool_m);
    j->errcode = 0;
    j->comp_len = 0;
    j->uncomp_len = 0;
    j->hit_eof = 0;

    while (bgzf_mt_read_block(fp, j) == 0) {
        // Dispatch
        hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_decode_func, j);

        // Check for command
        pthread_mutex_lock(&mt->command_m);
        switch (mt->command) {
        case SEEK:
            bgzf_mt_seek(fp);
            pthread_mutex_unlock(&mt->command_m);
            goto restart;

        case HAS_EOF:
            bgzf_mt_eof(fp);
            break;

        case CLOSE:
            pthread_cond_signal(&mt->command_c);
            pthread_mutex_unlock(&mt->command_m);
            hts_tpool_process_destroy(mt->out_queue);
            return NULL;

        default:
            break;
        }
        pthread_mutex_unlock(&mt->command_m);

        // Allocate buffer for next block
        pthread_mutex_lock(&mt->job_pool_m);
        j = pool_alloc(mt->job_pool);
        pthread_mutex_unlock(&mt->job_pool_m);
        j->errcode = 0;
        j->comp_len = 0;
        j->uncomp_len = 0;
        j->hit_eof = 0;
    }

    if (j->errcode == BGZF_ERR_MT) {
        // Attempt to multi-thread decode a raw gzip stream cannot be done.
        // We tear down the multi-threaded decoder and revert to the old code.
        hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_nul_func, j);
        hts_tpool_process_ref_decr(mt->out_queue);
        return &j->errcode;
    }

    // Dispatch an empty block so EOF is spotted.
    // We also use this mechanism for returning errors, in which case
    // j->errcode is set already.

    j->hit_eof = 1;
    hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_nul_func, j);
    if (j->errcode != 0) {
        hts_tpool_process_destroy(mt->out_queue);
        return &j->errcode;
    }

    // We hit EOF so can stop reading, but we may get a subsequent
    // seek request.  In this case we need to restart the reader.
    //
    // To handle this we wait on a condition variable and then
    // monitor the command. (This could be either seek or close.)
    for (;;) {
        pthread_mutex_lock(&mt->command_m);
        if (mt->command == NONE)
            pthread_cond_wait(&mt->command_c, &mt->command_m);
        switch(mt->command) {
        default:
            pthread_mutex_unlock(&mt->command_m);
            break;

        case SEEK:
            bgzf_mt_seek(fp);
            pthread_mutex_unlock(&mt->command_m);
            goto restart;

        case HAS_EOF:
            bgzf_mt_eof(fp);
            pthread_mutex_unlock(&mt->command_m);
            continue;

        case CLOSE:
            pthread_cond_signal(&mt->command_c);
            pthread_mutex_unlock(&mt->command_m);
            hts_tpool_process_destroy(mt->out_queue);
            return NULL;
        }
    }
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
                   fp->is_write ? bgzf_mt_writer : bgzf_mt_reader, fp);

    return 0;
}

int bgzf_mt(BGZF *fp, int n_threads, __attribute__((unused)) int n_sub_blks)
{
    // No gain from multi-threading when not compressed
    if (!fp->is_compressed)
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

int bgzf_mt(__attribute__((unused)) BGZF *fp, __attribute__((unused)) int n_threads, __attribute__((unused)) int n_sub_blks)
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
        block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed");
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

int bgzf_close(BGZF* fp)
{
    int ret, block_length;
    if (fp == 0) return -1;
    if (fp->is_write && fp->is_compressed) {
        if (bgzf_flush(fp) != 0) return -1;
        fp->compress_level = -1;
        block_length = deflate_block(fp, 0); // write an empty block
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed");
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
    ret = hclose(fp->fp);
    if (ret != 0) return -1;
    // bgzf_index_destroy(fp);
    free(fp->uncompressed_block);
    // free_cache(fp);
    free(fp);
    return 0;
}

int64_t bgzf_seek(BGZF* fp, int64_t pos, int where)
{
    int block_offset;
    int64_t block_address;

    if (fp->is_write || where != SEEK_SET || fp->is_gzip) {
        fp->errcode |= BGZF_ERR_MISUSE;
        return -1;
    }
    block_offset = pos & 0xFFFF;
    block_address = pos >> 16;

#ifdef BGZF_MT
    if (fp->mt) {
        // The reader runs asynchronous and does loops of:
        //    Read block
        //    Check & process command
        //    Dispatch decode job
        //
        // Once at EOF it then switches to loops of
        //    Wait for command
        //    Process command (possibly switching back to above loop).
        //
        // To seek we therefore send the reader thread a SEEK command,
        // waking it up if blocked in dispatch and signalling if
        // waiting for a command.  We then wait for the response so we
        // know the seek succeeded.
        pthread_mutex_lock(&fp->mt->command_m);
        fp->mt->hit_eof = 0;
        fp->mt->command = SEEK;
        fp->mt->block_address = block_address;
        pthread_cond_signal(&fp->mt->command_c);
        hts_tpool_wake_dispatch(fp->mt->out_queue);
        pthread_cond_wait(&fp->mt->command_c, &fp->mt->command_m);

        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = block_address;
        fp->block_offset = block_offset;

        pthread_mutex_unlock(&fp->mt->command_m);
    } else {
#endif
        if (hseek(fp->fp, block_address, SEEK_SET) < 0) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = block_address;
        fp->block_offset = block_offset;
#ifdef BGZF_MT
    }
#endif

    return 0;
}
