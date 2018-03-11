/// @file htslib/bgzf.h
/// Low-level routines for direct BGZF operations.
/*
   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013, 2014,2017 Genome Research Ltd

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

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

// This is a heavily modified copy of commit
// 77f002106602103150db45b06ddec0f022c0b6fb : most code which isn't needed to
// generate bgzf-compressed files with libdeflate has been removed.

#ifndef HTSLIB_BGZF_H
#define HTSLIB_BGZF_H

#include <stdint.h>
#include <stdio.h>

// Do not use zstd wrapper here.
// (This dependency will soon be replaced with a copy of libdeflate.)
#ifdef STATIC_ZLIB
#  include "../../../zlib-1.2.11/zlib.h"
#else
#  include <zlib.h>
#endif

#include <sys/types.h>

#include "hts_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BGZF_BLOCK_SIZE     0xff00 // make sure compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE 0x10000

#define BGZF_ERR_ZLIB   1
#define BGZF_ERR_HEADER 2
#define BGZF_ERR_IO     4
#define BGZF_ERR_MISUSE 8
#define BGZF_ERR_MT     16 // stream cannot be multi-threaded
#define BGZF_ERR_CRC    32

struct hFILE;
struct hts_tpool;
struct bgzf_mtaux_t;
typedef struct __bgzidx_t bgzidx_t;
typedef struct bgzf_cache_t bgzf_cache_t;

struct BGZF {
    // Reserved bits should be written as 0; read as "don't care"
    unsigned errcode:16, reserved:1, is_write:1, no_eof_block:1, is_be:1;
    signed compress_level:9;
    unsigned last_block_eof:1, is_compressed:1, is_gzip:1;
    int cache_size;
    int block_length, block_clength, block_offset;
    int64_t block_address, uncompressed_address;
    void *uncompressed_block, *compressed_block;
    bgzf_cache_t *cache;
    struct hFILE *fp; // actual file handle
    struct bgzf_mtaux_t *mt; // only used for multi-threading
    bgzidx_t *idx;      // BGZF index
    int idx_build_otf;  // build index on the fly, set by bgzf_index_build_init()
    z_stream *gz_stream;// for gzip-compressed files
};
#ifndef HTS_BGZF_TYPEDEF
typedef struct BGZF BGZF;
#define HTS_BGZF_TYPEDEF
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

    /******************
     * Basic routines *
     ******************/

    /**
     * Open the specified file for reading or writing.
     */
    BGZF* bgzf_open(const char* path, const char *mode);

    /**
     * Close the BGZF and free all associated resources.
     *
     * @param fp    BGZF file handler
     * @return      0 on success and -1 on error
     */
    int bgzf_close(BGZF *fp);

    /**
     * Write _length_ bytes from _data_ to the file.  If no I/O errors occur,
     * the complete _length_ bytes will be written (or queued for writing).
     *
     * @param fp     BGZF file handler
     * @param data   data array to write
     * @param length size of data to write
     * @return       number of bytes written (i.e., _length_); negative on error
     */
    ssize_t bgzf_write(BGZF *fp, const void *data, size_t length) HTS_RESULT_USED;

    /**
     * Write _length_ bytes from _data_ to the file, the index will be used to
     * decide the amount of uncompressed data to be writen to each bgzip block.
     * If no I/O errors occur, the complete _length_ bytes will be written (or
     * queued for writing).
     * @param fp     BGZF file handler
     * @param data   data array to write
     * @param length size of data to write
     * @return       number of bytes written (i.e., _length_); negative on error
     */
    ssize_t bgzf_block_write(BGZF *fp, const void *data, size_t length);

    /**
     * Write _length_ bytes directly to the underlying stream without
     * compressing.  Bypasses BGZF blocking, so must be used with care
     * in specialised circumstances only.
     *
     * @param fp     BGZF file handler
     * @param data   data array to write
     * @param length number of raw bytes to write
     * @return       number of bytes actually written; -1 on error
     */
    ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length) HTS_RESULT_USED;

    /**
     * Write the data in the buffer to the file.
     *
     * @param fp     BGZF file handle
     * @return       0 on success and -1 on error
     */
    int bgzf_flush(BGZF *fp) HTS_RESULT_USED;

    /**
     * Return a virtual file pointer to the current location in the file.
     * No interpretation of the value should be made, other than a subsequent
     * call to bgzf_seek can be used to position the file at the same point.
     * Return value is non-negative on success.
     */
    #define bgzf_tell(fp) (((fp)->block_address << 16) | ((fp)->block_offset & 0xFFFF))

    /*********************
     * Advanced routines *
     *********************/

    /**
     * Flush the file if the remaining buffer size is smaller than _size_
     * @return      0 if flushing succeeded or was not needed; negative on error
     */
    int bgzf_flush_try(BGZF *fp, ssize_t size) HTS_RESULT_USED;

    /**
     * Enable multi-threading (when compiled with -DBGZF_MT) via a shared
     * thread pool.  This means both encoder and decoder can balance
     * usage across a single pool of worker jobs.
     *
     * @param fp          BGZF file handler; must be opened for writing
     * @param pool        The thread pool (see hts_create_threads)
     */
    int bgzf_thread_pool(BGZF *fp, struct hts_tpool *pool, int qsize);

    /**
     * Enable multi-threading (only effective when the library was compiled
     * with -DBGZF_MT)
     *
     * @param fp          BGZF file handler; must be opened for writing
     * @param n_threads   #threads used for writing
     * @param n_sub_blks  #blocks processed by each thread; a value 64-256 is recommended
     */
    int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks);

    /**
     * Compress a single BGZF block.
     *
     * @param dst    output buffer (must have size >= BGZF_MAX_BLOCK_SIZE)
     * @param dlen   size of output buffer; updated on return to the number
     *               of bytes actually written to dst
     * @param src    buffer to be compressed
     * @param slen   size of data to compress (must be <= BGZF_BLOCK_SIZE)
     * @param level  compression level
     * @return       0 on success and negative on error
     */
    int bgzf_compress(void *dst, size_t *dlen, const void *src, size_t slen, int level);

    /*******************
     * bgzidx routines *
     *******************/

    /**
     * Tell BGZF to build index while compressing.
     *
     * @param fp          BGZF file handler; can be opened for reading or writing.
     *
     * Returns 0 on success and -1 on error.
     */
    int bgzf_index_build_init(BGZF *fp);

    /// Save BGZF index
    /**
     * @param fp          BGZF file handler
     * @param bname       base name
     * @param suffix      suffix to add to bname (can be NULL)
     * @return 0 on success and -1 on error.
     */
    int bgzf_index_dump(BGZF *fp,
                        const char *bname, const char *suffix) HTS_RESULT_USED;

#ifdef __cplusplus
}
#endif

#endif
