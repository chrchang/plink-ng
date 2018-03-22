/*  hfile.c -- buffered low-level input/output streams.

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

// This is a modified copy of commit
// 53241915fa8b6a3e807f2ac85d8701f27a7fe528 ; some code which isn't needed to
// read or generate bgzf-compressed files with libdeflate has been removed.

// #include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <limits.h>

#include <pthread.h>

#include "htslib/hfile.h"
#include "hfile_internal.h"

#ifndef ENOTSUP
#define ENOTSUP EINVAL
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif
#ifndef EPROTONOSUPPORT
#define EPROTONOSUPPORT ENOSYS
#endif

#ifndef SSIZE_MAX /* SSIZE_MAX is POSIX 1 */
#define SSIZE_MAX LONG_MAX
#endif

/* hFILE fields are used as follows:

   char *buffer;     // Pointer to the start of the I/O buffer
   char *begin;      // First not-yet-read character / unused position
   char *end;        // First unfilled/unfillable position
   char *limit;      // Pointer to the first position past the buffer

   const hFILE_backend *backend;  // Methods to refill/flush I/O buffer

   off_t offset;     // Offset within the stream of buffer position 0
   unsigned at_eof:1;// For reading, whether EOF has been seen
   unsigned mobile:1;// Buffer is a mobile window or fixed full contents
   unsigned readonly:1;// Whether opened as "r" rather than "r+"/"w"/"a"
   int has_errno;    // Error number from the last failure on this stream

For reading, begin is the first unread character in the buffer and end is the
first unfilled position:

   -----------ABCDEFGHIJKLMNO---------------
   ^buffer    ^begin         ^end           ^limit

For writing, begin is the first unused position and end is unused so remains
equal to buffer:

   ABCDEFGHIJKLMNOPQRSTUVWXYZ---------------
   ^buffer                   ^begin         ^limit
   ^end

Thus if begin > end then there is a non-empty write buffer, if begin < end
then there is a non-empty read buffer, and if begin == end then both buffers
are empty.  In all cases, the stream's file position indicator corresponds
to the position pointed to by begin.

The above is the normal scenario of a mobile window.  For in-memory
streams (eg via hfile_init_fixed) the buffer can be used as the full
contents without any separate backend behind it.  These always have at_eof
set, offset set to 0, need no read() method, and should just return EINVAL
for seek():

   abcdefghijkLMNOPQRSTUVWXYZ------
   ^buffer    ^begin         ^end  ^limit
*/

hFILE *hfile_init(size_t struct_size, const char *mode, size_t capacity)
{
    hFILE *fp = (hFILE *) malloc(struct_size);
    if (fp == NULL) goto error;

    if (capacity == 0) capacity = 32768;
    // FIXME For now, clamp input buffer sizes so mpileup doesn't eat memory
    if (strchr(mode, 'r') && capacity > 32768) capacity = 32768;

    fp->buffer = (char *) malloc(capacity);
    if (fp->buffer == NULL) goto error;

    fp->begin = fp->end = fp->buffer;
    fp->limit = &fp->buffer[capacity];

    fp->offset = 0;
    fp->at_eof = 0;
    fp->mobile = 1;
    fp->readonly = (strchr(mode, 'r') && ! strchr(mode, '+'));
    fp->has_errno = 0;
    return fp;

error:
    hfile_destroy(fp);
    return NULL;
}

void hfile_destroy(hFILE *fp)
{
    int save = errno;
    if (fp) free(fp->buffer);
    free(fp);
    errno = save;
}

static inline int writebuffer_is_nonempty(hFILE *fp)
{
    return fp->begin > fp->end;
}

/* Refills the read buffer from the backend (once, so may only partially
   fill the buffer), returning the number of additional characters read
   (which might be 0), or negative when an error occurred.  */
static ssize_t refill_buffer(hFILE *fp)
{
    ssize_t n;

    // Move any unread characters to the start of the buffer
    if (fp->mobile && fp->begin > fp->buffer) {
        fp->offset += fp->begin - fp->buffer;
        memmove(fp->buffer, fp->begin, fp->end - fp->begin);
        fp->end = &fp->buffer[fp->end - fp->begin];
        fp->begin = fp->buffer;
    }

    // Read into the available buffer space at fp->[end,limit)
    if (fp->at_eof || fp->end == fp->limit) n = 0;
    else {
        n = fp->backend->read(fp, fp->end, fp->limit - fp->end);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
    }

    fp->end += n;
    return n;
}

/*
 * Changes the buffer size for an hFILE.  Ideally this is done
 * immediately after opening.  If performed later, this function may
 * fail if we are reducing the buffer size and the current offset into
 * the buffer is beyond the new capacity.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int hfile_set_blksize(hFILE *fp, size_t bufsiz) {
    char *buffer;
    ptrdiff_t curr_used;
    if (!fp) return -1;
    curr_used = (fp->begin > fp->end ? fp->begin : fp->end) - fp->buffer;
    if (bufsiz == 0) bufsiz = 32768;

    // Ensure buffer resize will not erase live data
    if ((ptrdiff_t)bufsiz < curr_used)
        return -1;

    if (!(buffer = (char *) realloc(fp->buffer, bufsiz))) return -1;

    fp->begin  = buffer + (fp->begin - fp->buffer);
    fp->end    = buffer + (fp->end   - fp->buffer);
    fp->buffer = buffer;
    fp->limit  = &fp->buffer[bufsiz];

    return 0;
}

ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes)
{
    size_t n = fp->end - fp->begin;
    while (n < nbytes) {
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;
        else if (ret == 0) break;
        else n += ret;
    }

    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp->begin, n);
    return n;
}

/* Called only from hread(); when called, our buffer is empty and nread bytes
   have already been placed in the destination buffer.  */
ssize_t hread2(hFILE *fp, void *destv, size_t nbytes, size_t nread)
{
    const size_t capacity = fp->limit - fp->buffer;
    int buffer_invalidated = 0;
    char *dest = (char *) destv;
    dest += nread, nbytes -= nread;

    // Read large requests directly into the destination buffer
    while (nbytes * 2 >= capacity && !fp->at_eof) {
        ssize_t n = fp->backend->read(fp, dest, nbytes);
        if (n < 0) { fp->has_errno = errno; return n; }
        else if (n == 0) fp->at_eof = 1;
        else buffer_invalidated = 1;
        fp->offset += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    if (buffer_invalidated) {
        // Our unread buffer is empty, so begin == end, but our already-read
        // buffer [buffer,begin) is likely non-empty and is no longer valid as
        // its contents are no longer adjacent to the file position indicator.
        // Discard it so that hseek() can't try to take advantage of it.
        fp->offset += fp->begin - fp->buffer;
        fp->begin = fp->end = fp->buffer;
    }

    while (nbytes > 0 && !fp->at_eof) {
        size_t n;
        ssize_t ret = refill_buffer(fp);
        if (ret < 0) return ret;

        n = fp->end - fp->begin;
        if (n > nbytes) n = nbytes;
        memcpy(dest, fp->begin, n);
        fp->begin += n;
        dest += n, nbytes -= n;
        nread += n;
    }

    return nread;
}

/* Flushes the write buffer, fp->[buffer,begin), out through the backend
   returning 0 on success or negative if an error occurred.  */
static ssize_t flush_buffer(hFILE *fp)
{
    const char *buffer = fp->buffer;
    while (buffer < fp->begin) {
        ssize_t n = fp->backend->write(fp, buffer, fp->begin - buffer);
        if (n < 0) { fp->has_errno = errno; return n; }
        buffer += n;
        fp->offset += n;
    }

    fp->begin = fp->buffer;  // Leave the buffer empty
    return 0;
}

int hflush(hFILE *fp)
{
    if (flush_buffer(fp) < 0) return EOF;
    if (fp->backend->flush) {
        if (fp->backend->flush(fp) < 0) { fp->has_errno = errno; return EOF; }
    }
    return 0;
}

/* Called only from hwrite() and hputs2(); when called, our buffer is full and
   ncopied bytes from the source have already been copied to our buffer.  */
ssize_t hwrite2(hFILE *fp, const void *srcv, size_t totalbytes, size_t ncopied)
{
    const char *src = (const char *) srcv;
    ssize_t ret;
    const size_t capacity = fp->limit - fp->buffer;
    size_t remaining = totalbytes - ncopied;
    src += ncopied;

    ret = flush_buffer(fp);
    if (ret < 0) return ret;

    // Write large blocks out directly from the source buffer
    while (remaining * 2 >= capacity) {
        ssize_t n = fp->backend->write(fp, src, remaining);
        if (n < 0) { fp->has_errno = errno; return n; }
        fp->offset += n;
        src += n, remaining -= n;
    }

    // Just buffer any remaining characters
    memcpy(fp->begin, src, remaining);
    fp->begin += remaining;

    return totalbytes;
}

off_t hseek(hFILE *fp, off_t offset, int whence)
{
    off_t curpos, pos;

    if (writebuffer_is_nonempty(fp) && fp->mobile) {
        int ret = flush_buffer(fp);
        if (ret < 0) return ret;
    }

    curpos = htell(fp);

    // Relative offsets are given relative to the hFILE's stream position,
    // which may differ from the backend's physical position due to buffering
    // read-ahead.  Correct for this by converting to an absolute position.
    if (whence == SEEK_CUR) {
        if (curpos + offset < 0) {
            // Either a negative offset resulted in a position before the
            // start of the file, or we overflowed when given a positive offset
            fp->has_errno = errno = (offset < 0)? EINVAL : EOVERFLOW;
            return -1;
        }

        whence = SEEK_SET;
        offset = curpos + offset;
    }
    // For fixed immobile buffers, convert everything else to SEEK_SET too
    // so that seeking can be avoided for all (within range) requests.
    else if (! fp->mobile && whence == SEEK_END) {
        size_t length = fp->end - fp->buffer;
        if (offset > 0 || -offset > (off_t)length) {
            fp->has_errno = errno = EINVAL;
            return -1;
        }

        whence = SEEK_SET;
        offset = length + offset;
    }

    // Avoid seeking if the desired position is within our read buffer.
    // (But not when the next operation may be a write on a mobile buffer.)
    if (whence == SEEK_SET && (! fp->mobile || fp->readonly) &&
        offset >= fp->offset && offset - fp->offset <= fp->end - fp->buffer) {
        fp->begin = &fp->buffer[offset - fp->offset];
        return offset;
    }

    pos = fp->backend->seek(fp, offset, whence);
    if (pos < 0) { fp->has_errno = errno; return pos; }

    // Seeking succeeded, so discard any non-empty read buffer
    fp->begin = fp->end = fp->buffer;
    fp->at_eof = 0;

    fp->offset = pos;
    return pos;
}

int hclose(hFILE *fp)
{
    int err = fp->has_errno;

    if (writebuffer_is_nonempty(fp) && hflush(fp) < 0) err = fp->has_errno;
    if (fp->backend->close(fp) < 0) err = errno;
    hfile_destroy(fp);

    if (err) {
        errno = err;
        return EOF;
    }
    else return 0;
}

void hclose_abruptly(hFILE *fp)
{
    int save = errno;
    if (fp->backend->close(fp) < 0) { /* Ignore subsequent errors */ }
    hfile_destroy(fp);
    errno = save;
}


/***************************
 * File descriptor backend *
 ***************************/

#ifndef _WIN32
#include <sys/socket.h>
#include <sys/stat.h>
#define HAVE_STRUCT_STAT_ST_BLKSIZE
#else
/*
#include <winsock2.h>
#define HAVE_CLOSESOCKET
#define HAVE_SETMODE
*/
#endif
#include <fcntl.h>
#include <unistd.h>

/* For Unix, it doesn't matter whether a file descriptor is a socket.
   However Windows insists on send()/recv() and its own closesocket()
   being used when fd happens to be a socket.  */

typedef struct {
    hFILE base;
    int fd;
    unsigned is_socket:1;
} hFILE_fd;

static ssize_t fd_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket? recv(fp->fd, buffer, nbytes, 0)
                         : read(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
    return n;
}

static ssize_t fd_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    ssize_t n;
    do {
        n = fp->is_socket?  send(fp->fd, buffer, nbytes, 0)
                         : write(fp->fd, buffer, nbytes);
    } while (n < 0 && errno == EINTR);
#ifdef _WIN32
        // On windows we have no SIGPIPE.  Instead write returns
        // EINVAL.  We check for this and our fd being a pipe.
        // If so, we raise SIGTERM instead of SIGPIPE.  It's not
        // ideal, but I think the only alternative is extra checking
        // in every single piece of code.
        if (n < 0 && errno == EINVAL &&
            GetLastError() == ERROR_NO_DATA &&
            GetFileType((HANDLE)_get_osfhandle(fp->fd)) == FILE_TYPE_PIPE) {
            raise(SIGTERM);
        }
#endif
    return n;
}

static off_t fd_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    return lseek(fp->fd, offset, whence);
}

static int fd_flush(__attribute__((unused)) hFILE *fpv)
{
    int ret = 0;
    do {
#ifdef HAVE_FDATASYNC
        hFILE_fd *fp = (hFILE_fd *) fpv;
        ret = fdatasync(fp->fd);
#elif defined(HAVE_FSYNC)
        hFILE_fd *fp = (hFILE_fd *) fpv;
        ret = fsync(fp->fd);
#endif
        // Ignore invalid-for-fsync(2) errors due to being, e.g., a pipe,
        // and operation-not-supported errors (Mac OS X)
        if (ret < 0 && (errno == EINVAL || errno == ENOTSUP)) ret = 0;
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static int fd_close(hFILE *fpv)
{
    hFILE_fd *fp = (hFILE_fd *) fpv;
    int ret;
    do {
#ifdef HAVE_CLOSESOCKET
        ret = fp->is_socket? closesocket(fp->fd) : close(fp->fd);
#else
        ret = close(fp->fd);
#endif
    } while (ret < 0 && errno == EINTR);
    return ret;
}

static const struct hFILE_backend fd_backend =
{
    fd_read, fd_write, fd_seek, fd_flush, fd_close
};

static size_t blksize(int fd)
{
#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
    struct stat sbuf;
    if (fstat(fd, &sbuf) != 0) return 0;
    return sbuf.st_blksize;
#else
    return 0;
#endif
}

static hFILE *hopen_fd(const char *filename, const char *mode)
{
    hFILE_fd *fp = NULL;
    int fd = open(filename, hfile_oflags(mode), 0666);
    if (fd < 0) goto error;

    fp = (hFILE_fd *) hfile_init(sizeof (hFILE_fd), mode, blksize(fd));
    if (fp == NULL) goto error;

    fp->fd = fd;
    fp->is_socket = 0;
    fp->base.backend = &fd_backend;
    return &fp->base;

error:
    if (fd >= 0) { int save = errno; (void) close(fd); errno = save; }
    hfile_destroy((hFILE *) fp);
    return NULL;
}

int hfile_oflags(const char *mode)
{
    int rdwr = 0, flags = 0;
    const char *s;
    for (s = mode; *s; s++)
        switch (*s) {
        case 'r': rdwr = O_RDONLY;  break;
        case 'w': rdwr = O_WRONLY; flags |= O_CREAT | O_TRUNC;  break;
        case 'a': rdwr = O_WRONLY; flags |= O_CREAT | O_APPEND;  break;
        case '+': rdwr = O_RDWR;  break;
#ifdef O_CLOEXEC
        case 'e': flags |= O_CLOEXEC;  break;
#endif
#ifdef O_EXCL
        case 'x': flags |= O_EXCL;  break;
#endif
        default:  break;
        }

#ifdef O_BINARY
    flags |= O_BINARY;
#endif

    return rdwr | flags;
}

hFILE *hopen(const char *fname, const char *mode, ...)
{
    return hopen_fd(fname, mode);
}
