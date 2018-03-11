/// @file htslib/hfile.h
/// Buffered low-level input/output streams.
/*
    Copyright (C) 2013-2016 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

// This is a heavily modified copy of commit
// 77f002106602103150db45b06ddec0f022c0b6fb : most code which isn't needed to
// generate bgzf-compressed files with libdeflate has been removed.

#ifndef HTSLIB_HFILE_H
#define HTSLIB_HFILE_H

#include <string.h>

#include <sys/types.h>

#include "hts_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct hFILE_backend;
/// Low-level input/output stream handle
/** The fields of this structure are declared here solely for the benefit
of the hFILE-related inline functions.  They may change in future releases.
User code should not use them directly; you should imagine that hFILE is an
opaque incomplete type.
*/
typedef struct hFILE {
    // @cond internal
    char *buffer, *begin, *end, *limit;
    const struct hFILE_backend *backend;
    off_t offset;
    unsigned at_eof:1, mobile:1, readonly:1;
    int has_errno;
    // @endcond
} hFILE;

/// Open the named file or URL as a stream
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.

The usual `fopen(3)` _mode_ letters are supported: one of
`r` (read), `w` (write), `a` (append), optionally followed by any of
`+` (update), `e` (close on `exec(2)`), `x` (create exclusively),
`:` (indicates scheme-specific variable arguments follow).
*/
hFILE *hopen(const char *filename, const char *mode, ...) HTS_RESULT_USED;

/// Flush (for output streams) and close the stream
/** @return  0 if successful, or `EOF` (with _errno_ set) if an error occurred.
*/
int hclose(hFILE *fp) HTS_RESULT_USED;

/// Close the stream, without flushing or propagating errors
/** For use while cleaning up after an error only.  Preserves _errno_.
*/
void hclose_abruptly(hFILE *fp);

/// Write a block of characters to the file
/** @return  Either _nbytes_, or negative if an error occurred.

In the absence of I/O errors, the full _nbytes_ will be written.
*/
static inline ssize_t HTS_RESULT_USED
hwrite(hFILE *fp, const void *buffer, size_t nbytes)
{
    extern ssize_t hwrite2(hFILE *, const void *, size_t, size_t);
    extern int hfile_set_blksize(hFILE *fp, size_t bufsiz);

    if(!fp->mobile){
      if (fp->limit - fp->begin < (ssize_t)nbytes){
            hfile_set_blksize(fp, fp->limit - fp->buffer + nbytes);
            fp->end = fp->limit;
        }
    }

    size_t n = fp->limit - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(fp->begin, buffer, n);
    fp->begin += n;
    return (n==nbytes)? (ssize_t) n : hwrite2(fp, buffer, nbytes, n);
}

/// For writing streams, flush buffered output to the underlying stream
/** @return  0 if successful, or `EOF` if an error occurred.

This includes low-level flushing such as via `fdatasync(2)`.
*/
int hflush(hFILE *fp) HTS_RESULT_USED;

#ifdef __cplusplus
}
#endif

#endif
