/* pigz.c -- parallel implementation of gzip
 * Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2013 Mark Adler
 * Version 2.3  3 Mar 2013  Mark Adler
 */

// Modified by Christopher Chang (chrchang@alumni.caltech.edu) to export
// parallel_compress() as a library function (except on Windows, where the
// yarn.c threadling library doesn't yet work).

/*
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the author be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.  (THIS IS AN ALTERED
     VERSION.)
  3. This notice may not be removed or altered from any source distribution.

  Mark Adler
  madler@alumni.caltech.edu

  Mark accepts donations for providing this software.  Donations are not
  required or expected.  Any amount that you feel is appropriate would be
  appreciated.  You can use this link:

  https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=536055

 */

/* Version history:
   1.0    17 Jan 2007  First version, pipe only
   1.1    28 Jan 2007  Avoid void * arithmetic (some compilers don't get that)
                       Add note about requiring zlib 1.2.3
                       Allow compression level 0 (no compression)
                       Completely rewrite parallelism -- add a write thread
                       Use deflateSetDictionary() to make use of history
                       Tune argument defaults to best performance on four cores
   1.2.1   1 Feb 2007  Add long command line options, add all gzip options
                       Add debugging options
   1.2.2  19 Feb 2007  Add list (--list) function
                       Process file names on command line, write .gz output
                       Write name and time in gzip header, set output file time
                       Implement all command line options except --recursive
                       Add --keep option to prevent deleting input files
                       Add thread tracing information with -vv used
                       Copy crc32_combine() from zlib (shared libraries issue)
   1.3    25 Feb 2007  Implement --recursive
                       Expand help to show all options
                       Show help if no arguments or output piping are provided
                       Process options in GZIP environment variable
                       Add progress indicator to write thread if --verbose
   1.4     4 Mar 2007  Add --independent to facilitate damaged file recovery
                       Reallocate jobs for new --blocksize or --processes
                       Do not delete original if writing to stdout
                       Allow --processes 1, which does no threading
                       Add NOTHREAD define to compile without threads
                       Incorporate license text from zlib in source code
   1.5    25 Mar 2007  Reinitialize jobs for new compression level
                       Copy attributes and owner from input file to output file
                       Add decompression and testing
                       Add -lt (or -ltv) to show all entries and proper lengths
                       Add decompression, testing, listing of LZW (.Z) files
                       Only generate and show trace log if DEBUG defined
                       Take "-" argument to mean read file from stdin
   1.6    30 Mar 2007  Add zlib stream compression (--zlib), and decompression
   1.7    29 Apr 2007  Decompress first entry of a zip file (if deflated)
                       Avoid empty deflate blocks at end of deflate stream
                       Show zlib check value (Adler-32) when listing
                       Don't complain when decompressing empty file
                       Warn about trailing junk for gzip and zlib streams
                       Make listings consistent, ignore gzip extra flags
                       Add zip stream compression (--zip)
   1.8    13 May 2007  Document --zip option in help output
   2.0    19 Oct 2008  Complete rewrite of thread usage and synchronization
                       Use polling threads and a pool of memory buffers
                       Remove direct pthread library use, hide in yarn.c
   2.0.1  20 Oct 2008  Check version of zlib at compile time, need >= 1.2.3
   2.1    24 Oct 2008  Decompress with read, write, inflate, and check threads
                       Remove spurious use of ctime_r(), ctime() more portable
                       Change application of job->calc lock to be a semaphore
                       Detect size of off_t at run time to select %lu vs. %llu
                       #define large file support macro even if not __linux__
                       Remove _LARGEFILE64_SOURCE, _FILE_OFFSET_BITS is enough
                       Detect file-too-large error and report, blame build
                       Replace check combination routines with those from zlib
   2.1.1  28 Oct 2008  Fix a leak for files with an integer number of blocks
                       Update for yarn 1.1 (yarn_prefix and yarn_abort)
   2.1.2  30 Oct 2008  Work around use of beta zlib in production systems
   2.1.3   8 Nov 2008  Don't use zlib combination routines, put back in pigz
   2.1.4   9 Nov 2008  Fix bug when decompressing very short files
   2.1.5  20 Jul 2009  Added 2008, 2009 to --license statement
                       Allow numeric parameter immediately after -p or -b
                       Enforce parameter after -p, -b, -s, before other options
                       Enforce numeric parameters to have only numeric digits
                       Try to determine the number of processors for -p default
                       Fix --suffix short option to be -S to match gzip [Bloch]
                       Decompress if executable named "unpigz" [Amundsen]
                       Add a little bit of testing to Makefile
   2.1.6  17 Jan 2010  Added pigz.spec to distribution for RPM systems [Brown]
                       Avoid some compiler warnings
                       Process symbolic links if piping to stdout [Hoffstätte]
                       Decompress if executable named "gunzip" [Hoffstätte]
                       Allow ".tgz" suffix [Chernookiy]
                       Fix adler32 comparison on .zz files
   2.1.7  17 Dec 2011  Avoid unused parameter warning in reenter()
                       Don't assume 2's complement ints in compress_thread()
                       Replicate gzip -cdf cat-like behavior
                       Replicate gzip -- option to suppress option decoding
                       Test output from make test instead of showing it
                       Updated pigz.spec to install unpigz, pigz.1 [Obermaier]
                       Add PIGZ environment variable [Mueller]
                       Replicate gzip suffix search when decoding or listing
                       Fix bug in load() to set in_left to zero on end of file
                       Do not check suffix when input file won't be modified
                       Decompress to stdout if name is "*cat" [Hayasaka]
                       Write data descriptor signature to be like Info-ZIP
                       Update and sort options list in help
                       Use CC variable for compiler in Makefile
                       Exit with code 2 if a warning has been issued
                       Fix thread synchronization problem when tracing
                       Change macro name MAX to MAX2 to avoid library conflicts
                       Determine number of processors on HP-UX [Lloyd]
   2.2    31 Dec 2011  Check for expansion bound busting (e.g. modified zlib)
                       Make the "threads" list head global variable volatile
                       Fix construction and printing of 32-bit check values
                       Add --rsyncable functionality
   2.2.1   1 Jan 2012  Fix bug in --rsyncable buffer management
   2.2.2   1 Jan 2012  Fix another bug in --rsyncable buffer management
   2.2.3  15 Jan 2012  Remove volatile in yarn.c
                       Reduce the number of input buffers
                       Change initial rsyncable hash to comparison value
                       Improve the efficiency of arriving at a byte boundary
                       Add thread portability #defines from yarn.c
                       Have rsyncable compression be independent of threading
                       Fix bug where constructed dictionaries not being used
   2.2.4  11 Mar 2012  Avoid some return value warnings
                       Improve the portability of printing the off_t type
                       Check for existence of compress binary before using
                       Update zlib version checking to 1.2.6 for new functions
                       Fix bug in zip (-K) output
                       Fix license in pigz.spec
                       Remove thread portability #defines in pigz.c
   2.2.5  28 Jul 2012  Avoid race condition in free_pool()
                       Change suffix to .tar when decompressing or listing .tgz
                       Print name of executable in error messages
                       Show help properly when the name is unpigz or gunzip
                       Fix permissions security problem before output is closed
   2.3     3 Mar 2013  Don't complain about missing suffix when not writing output file
                       Put all global variables in one global structure for readability
                       Do not decompress concatenated zlib streams -- only gzip streams
                       Add option for compression level 11 to use zopfli
                       Fix handling of junk after compressed data
 */

/* To-do:
    - make source portable for Windows, VMS, etc. (see gzip source code)
    - make build portable (currently good for Unixish)
 */

/*
   pigz compresses using threads to make use of multiple processors and cores.
   The input is broken up into 128 KB chunks with each compressed in parallel.
   The individual check value for each chunk is also calculated in parallel.
   The compressed data is written in order to the output, and a combined check
   value is calculated from the individual check values.

   The compressed data format generated is in the gzip, zlib, or single-entry
   zip format using the deflate compression method.  The compression produces
   partial raw deflate streams which are concatenated by a single write thread
   and wrapped with the appropriate header and trailer, where the trailer
   contains the combined check value.

   Each partial raw deflate stream is terminated by an empty stored block
   (using the Z_SYNC_FLUSH option of zlib), in order to end that partial bit
   stream at a byte boundary, unless that partial stream happens to already end
   at a byte boundary (the latter requires zlib 1.2.6 or later).  Ending on a
   byte boundary allows the partial streams to be concatenated simply as
   sequences of bytes.  This adds a very small four to five byte overhead
   (average 3.75 bytes) to the output for each input chunk.

   The default input block size is 128K, but can be changed with the -b option.
   The number of compress threads is set by default to 8, which can be changed
   using the -p option.  Specifying -p 1 avoids the use of threads entirely.
   pigz will try to determine the number of processors in the machine, in which
   case if that number is two or greater, pigz will use that as the default for
   -p instead of 8.

   The input blocks, while compressed independently, have the last 32K of the
   previous block loaded as a preset dictionary to preserve the compression
   effectiveness of deflating in a single thread.  This can be turned off using
   the --independent or -i option, so that the blocks can be decompressed
   independently for partial error recovery or for random access.

   Decompression can't be parallelized, at least not without specially prepared
   deflate streams for that purpose.  As a result, pigz uses a single thread
   (the main thread) for decompression, but will create three other threads for
   reading, writing, and check calculation, which can speed up decompression
   under some circumstances.  Parallel decompression can be turned off by
   specifying one process (-dp 1 or -tp 1).

   pigz requires zlib 1.2.1 or later to allow setting the dictionary when doing
   raw deflate.  Since zlib 1.2.3 corrects security vulnerabilities in zlib
   version 1.2.1 and 1.2.2, conditionals check for zlib 1.2.3 or later during
   the compilation of pigz.c.  zlib 1.2.4 includes some improvements to
   Z_FULL_FLUSH and deflateSetDictionary() that permit identical output for
   pigz with and without threads, which is not possible with zlib 1.2.3.  This
   may be important for uses of pigz -R where small changes in the contents
   should result in small changes in the archive for rsync.  Note that due to
   the details of how the lower levels of compression result in greater speed,
   compression level 3 and below does not permit identical pigz output with
   and without threads.

   pigz uses the POSIX pthread library for thread control and communication,
   through the yarn.h interface to yarn.c.  yarn.c can be replaced with
   equivalent implementations using other thread libraries.  pigz can be
   compiled with NOTHREAD #defined to not use threads at all (in which case
   pigz will not be able to live up to the "parallel" in its name).
 */

/*
   Details of parallel compression implementation:

   When doing parallel compression, pigz uses the main thread to read the input
   in 'size' sized chunks (see -b), and puts those in a compression job list,
   each with a sequence number to keep track of the ordering.  If it is not the
   first chunk, then that job also points to the previous input buffer, from
   which the last 32K will be used as a dictionary (unless -i is specified).
   This sets a lower limit of 32K on 'size'.

   pigz launches up to 'procs' compression threads (see -p).  Each compression
   thread continues to look for jobs in the compression list and perform those
   jobs until instructed to return.  When a job is pulled, the dictionary, if
   provided, will be loaded into the deflate engine and then that input buffer
   is dropped for reuse.  Then the input data is compressed into an output
   buffer that grows in size if necessary to hold the compressed data. The job
   is then put into the write job list, sorted by the sequence number. The
   compress thread however continues to calculate the check value on the input
   data, either a CRC-32 or Adler-32, possibly in parallel with the write
   thread writing the output data.  Once that's done, the compress thread drops
   the input buffer and also releases the lock on the check value so that the
   write thread can combine it with the previous check values.  The compress
   thread has then completed that job, and goes to look for another.

   All of the compress threads are left running and waiting even after the last
   chunk is processed, so that they can support the next input to be compressed
   (more than one input file on the command line).  Once pigz is done, it will
   call all the compress threads home (that'll do pig, that'll do).

   Before starting to read the input, the main thread launches the write thread
   so that it is ready pick up jobs immediately.  The compress thread puts the
   write jobs in the list in sequence sorted order, so that the first job in
   the list is always has the lowest sequence number.  The write thread waits
   for the next write job in sequence, and then gets that job.  The job still
   holds its input buffer, from which the write thread gets the input buffer
   length for use in check value combination.  Then the write thread drops that
   input buffer to allow its reuse.  Holding on to the input buffer until the
   write thread starts also has the benefit that the read and compress threads
   can't get way ahead of the write thread and build up a large backlog of
   unwritten compressed data.  The write thread will write the compressed data,
   drop the output buffer, and then wait for the check value to be unlocked
   by the compress thread.  Then the write thread combines the check value for
   this chunk with the total check value for eventual use in the trailer.  If
   this is not the last chunk, the write thread then goes back to look for the
   next output chunk in sequence.  After the last chunk, the write thread
   returns and joins the main thread.  Unlike the compress threads, a new write
   thread is launched for each input stream.  The write thread writes the
   appropriate header and trailer around the compressed data.

   The input and output buffers are reused through their collection in pools.
   Each buffer has a use count, which when decremented to zero returns the
   buffer to the respective pool.  Each input buffer has up to three parallel
   uses: as the input for compression, as the data for the check value
   calculation, and as a dictionary for compression.  Each output buffer has
   only one use, which is as the output of compression followed serially as
   data to be written.  The input pool is limited in the number of buffers, so
   that reading does not get way ahead of compression and eat up memory with
   more input than can be used.  The limit is approximately two times the
   number of compression threads.  In the case that reading is fast as compared
   to compression, that number allows a second set of buffers to be read while
   the first set of compressions are being performed.  The number of output
   buffers is not directly limited, but is indirectly limited by the release of
   input buffers to about the same number.
 */

#ifdef _WIN32
// stopgap non-parallel code for Windows

#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#ifdef DYNAMIC_ZLIB
  #include <zlib.h>
#else
  #include "../zlib-1.2.11/zlib.h"
#endif

#include "pigz.h"

#ifdef _WIN64
  #define putc_unlocked _fputc_nolock
#else
  #define putc_unlocked putc
#endif

void pigz_init(uint32_t setprocs) {
  return;
}

void parallel_compress(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, uint32_t(* emitn)(uint32_t, unsigned char*)) {
  // minor issue: this currently writes \n instead of \r\n linebreaks.
  uint32_t overflow_ct = 0;
  gzFile gz_outfile = gzopen(out_fname, do_append? "ab": "wb");
  unsigned char* write_ptr;
  uint32_t last_size;
  if (!gz_outfile) {
    putc_unlocked('\n', stdout);
    fflush(stdout);
    fprintf(stderr, "Error: Failed to open %s.\n", out_fname);
    exit(2);
  }
  do {
    last_size = emitn(overflow_ct, overflow_buf);
    if (last_size > PIGZ_BLOCK_SIZE) {
      overflow_ct = last_size - PIGZ_BLOCK_SIZE;
      last_size = PIGZ_BLOCK_SIZE;
    } else {
      overflow_ct = 0;
    }
    if (last_size) {
      if (!gzwrite(gz_outfile, overflow_buf, last_size)) {
	putc_unlocked('\n', stdout);
	fflush(stdout);
	fputs("Error: File write failure.\n", stderr);
	gzclose(gz_outfile);
	exit(6);
      }
    }
    if (overflow_ct) {
      write_ptr = &(overflow_buf[PIGZ_BLOCK_SIZE]);
      while (overflow_ct > PIGZ_BLOCK_SIZE) {
	if (!gzwrite(gz_outfile, write_ptr, PIGZ_BLOCK_SIZE)) {
	  putc_unlocked('\n', stdout);
	  fflush(stdout);
	  fputs("Error: File write failure.\n", stderr);
	  gzclose(gz_outfile);
	  exit(6);
	}
	write_ptr = &(write_ptr[PIGZ_BLOCK_SIZE]);
	overflow_ct -= PIGZ_BLOCK_SIZE;
      }
      memcpy(overflow_buf, write_ptr, overflow_ct);
    }
  } while (last_size);
  if (gzclose(gz_outfile) != Z_OK) {
    putc_unlocked('\n', stdout);
    fflush(stdout);
    fputs("Error: File write failure.\n", stderr);
    exit(6);
  }
}

int32_t pzwrite_init(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr) {
    ps_ptr->outfile = fopen(out_fname, do_append? "ab" : "wb");
    ps_ptr->gz_outfile = NULL;
    if (!ps_ptr->outfile) {
        putc_unlocked('\n', stdout);
	fflush(stdout);
        fprintf(stderr, "Error: Failed to open %s.\n", out_fname);
        return 2; // RET_OPEN_FAIL
    }
    ps_ptr->overflow_buf = overflow_buf;
    return 0;
}

void compressed_pzwrite_init(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr) {
    ps_ptr->outfile = NULL;
    ps_ptr->gz_outfile = gzopen(out_fname, do_append? "ab" : "wb");
    if (!ps_ptr->gz_outfile) {
        putc_unlocked('\n', stdout);
        fflush(stdout);
        fprintf(stderr, "Error: Failed to open %s.\n", out_fname);
        exit(2);
    }
    ps_ptr->overflow_buf = overflow_buf;
}

int32_t flex_pzwrite_init(uint32_t output_gz, char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr) {
    if (!output_gz) {
        return pzwrite_init(out_fname, overflow_buf, do_append, ps_ptr);
    } else {
        compressed_pzwrite_init(out_fname, overflow_buf, do_append, ps_ptr);
        return 0;
    }
}

int32_t force_pzwrite(Pigz_state* ps_ptr, char** writep_ptr, uint32_t write_min) {
    unsigned char* writep = (unsigned char*)(*writep_ptr);
    if (ps_ptr->overflow_buf != writep) {
        if (!fwrite(ps_ptr->overflow_buf, writep - ps_ptr->overflow_buf, 1, ps_ptr->outfile)) {
	    return 6; // RET_WRITE_FAIL
	}
        *writep_ptr = (char*)(ps_ptr->overflow_buf);
    }
    return 0;
}

void force_compressed_pzwrite(Pigz_state* ps_ptr, char** writep_ptr, uint32_t write_min) {
    unsigned char* writep = (unsigned char*)(*writep_ptr);
    if (ps_ptr->overflow_buf != writep) {
        if (!gzwrite(ps_ptr->gz_outfile, ps_ptr->overflow_buf, writep - ps_ptr->overflow_buf)) {
	    putc_unlocked('\n', stdout);
	    fflush(stdout);
	    fputs("Error: File write failure.\n", stderr);
            gzclose(ps_ptr->gz_outfile);
            exit(6);
        }
        *writep_ptr = (char*)(ps_ptr->overflow_buf);
    }
}

int32_t flex_pzputs_std(Pigz_state* ps_ptr, char** writep_ptr, char* ss, uint32_t sslen) {
    unsigned char* writep = (unsigned char*)(*writep_ptr);
    unsigned char* readp = (unsigned char*)ss;
    uint32_t cur_write_space = 2 * PIGZ_BLOCK_SIZE - ((uintptr_t)(writep - ps_ptr->overflow_buf));
    while (sslen > cur_write_space) {
        memcpy(writep, readp, cur_write_space);
	if (is_uncompressed_pzwrite(ps_ptr)) {
	    if (!fwrite(ps_ptr->overflow_buf, 2 * PIGZ_BLOCK_SIZE, 1, ps_ptr->outfile)) {
                return 6;
	    }
	} else {
	    if (!gzwrite(ps_ptr->gz_outfile, ps_ptr->overflow_buf, 2 * PIGZ_BLOCK_SIZE)) {
	        putc_unlocked('\n', stdout);
		fflush(stdout);
	        fputs("Error: File write failure.\n", stderr);
	        gzclose(ps_ptr->gz_outfile);
	        exit(6);
	    }
	}
        writep = ps_ptr->overflow_buf;
        readp = &(readp[cur_write_space]);
	sslen -= cur_write_space;
	cur_write_space = 2 * PIGZ_BLOCK_SIZE;
    }
    memcpy(writep, readp, sslen);
    *writep_ptr = (char*)(&(writep[sslen]));
    return flex_pzwrite(ps_ptr, writep_ptr);
}

int32_t pzwrite_close_null(Pigz_state* ps_ptr, char* writep) {
    force_pzwrite(ps_ptr, &writep, 0);
    int32_t ii = ferror(ps_ptr->outfile);
    int32_t jj = fclose(ps_ptr->outfile);
    ps_ptr->overflow_buf = NULL;
    return ii || jj;
}

void compressed_pzwrite_close_null(Pigz_state* ps_ptr, char* writep) {
    force_compressed_pzwrite(ps_ptr, &writep, 0);
    ps_ptr->overflow_buf = NULL;
    if (gzclose(ps_ptr->gz_outfile) != Z_OK) {
        putc_unlocked('\n', stdout);
	fflush(stdout);
        fputs("Error: File write failure.\n", stderr);
        exit(6);
    }
}

int32_t flex_pzwrite_close_null(Pigz_state* ps_ptr, char* writep) {
    if (is_uncompressed_pzwrite(ps_ptr)) {
        return pzwrite_close_null(ps_ptr, writep);
    } else {
        compressed_pzwrite_close_null(ps_ptr, writep);
        return 0;
    }
}
#else

#define VERSION "pigz 2.3\n"

/* use large file functions if available */
#define _FILE_OFFSET_BITS 64

/* included headers and what is expected from each */
#include <stdio.h>      /* fflush(), fprintf(), fputs(), getchar(), putc(), */
                        /* puts(), printf(), vasprintf(), stderr, EOF, NULL,
                           SEEK_END, size_t, off_t */
#include <stdlib.h>     /* exit(), malloc(), free(), realloc(), atol(), */
                        /* atoi(), getenv() */
#include <stdarg.h>     /* va_start(), va_end(), va_list */
#include <string.h>     /* memset(), memchr(), memcpy(), strcmp(), strcpy() */
                        /* strncpy(), strlen(), strcat(), strrchr() */
#include <errno.h>      /* errno, EEXIST */
#include <assert.h>     /* assert() */
#include <time.h>       /* ctime(), time(), time_t, mktime() */
#include <signal.h>     /* signal(), SIGINT */
#include <sys/types.h>  /* ssize_t */
#include <sys/stat.h>   /* chmod(), stat(), fstat(), lstat(), struct stat, */
                        /* S_IFDIR, S_IFLNK, S_IFMT, S_IFREG */
#include <sys/time.h>   /* utimes(), gettimeofday(), struct timeval */
#include <unistd.h>     /* unlink(), _exit(), read(), write(), close(), */
                        /* lseek(), isatty(), chown() */
#include <fcntl.h>      /* open(), O_CREAT, O_EXCL, O_RDONLY, O_TRUNC, */
                        /* O_WRONLY */
#include <dirent.h>     /* opendir(), readdir(), closedir(), DIR, */
                        /* struct dirent */
#include <limits.h>     /* PATH_MAX, UINT_MAX */
#if __STDC_VERSION__-0 >= 199901L || __GNUC__-0 >= 3
#  include <inttypes.h> /* intmax_t */
#endif

#ifdef __hpux
#  include <sys/param.h>
#  include <sys/pstat.h>
#endif

#ifdef DYNAMIC_ZLIB
  #include <zlib.h>
#else
  #include "../zlib-1.2.11/zlib.h" /* deflateInit2(), deflateReset(), deflate(), */
                        /* deflateEnd(), deflateSetDictionary(), crc32(),
                           inflateBackInit(), inflateBack(), inflateBackEnd(),
                           Z_DEFAULT_COMPRESSION, Z_DEFAULT_STRATEGY,
                           Z_DEFLATED, Z_NO_FLUSH, Z_NULL, Z_OK,
                           Z_SYNC_FLUSH, z_stream */
#endif
#if !defined(ZLIB_VERNUM) || ZLIB_VERNUM < 0x1230
#  error Need zlib version 1.2.3 or later
#endif

#ifndef NOTHREAD
#  include "yarn.h"     /* thread, launch(), join(), join_all(), */
                        /* lock, new_lock(), possess(), twist(), wait_for(),
                           release(), peek_lock(), free_lock(), yarn_name */
#endif

#include "pigz.h"


/* for local functions and globals */
#define local static

/* prevent end-of-line conversions on MSDOSish operating systems */
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <io.h>       /* setmode(), O_BINARY */
#  define SET_BINARY_MODE(fd) setmode(fd, O_BINARY)
#else
#  define SET_BINARY_MODE(fd)
#endif

/* release an allocated pointer, if allocated, and mark as unallocated */
#define RELEASE(ptr) \
    do { \
        if ((ptr) != NULL) { \
            free(ptr); \
            ptr = NULL; \
        } \
    } while (0)

/* sliding dictionary size for deflate */
#define DICT 32768U

/* largest power of 2 that fits in an unsigned int -- used to limit requests
   to zlib functions that use unsigned int lengths */
#define MAXP2 (UINT_MAX - (UINT_MAX >> 1))

/* rsyncable constants -- RSYNCBITS is the number of bits in the mask for
   comparison.  For random input data, there will be a hit on average every
   1<<RSYNCBITS bytes.  So for an RSYNCBITS of 12, there will be an average of
   one hit every 4096 bytes, resulting in a mean block size of 4096.  RSYNCMASK
   is the resulting bit mask.  RSYNCHIT is what the hash value is compared to
   after applying the mask.

   The choice of 12 for RSYNCBITS is consistent with the original rsyncable
   patch for gzip which also uses a 12-bit mask.  This results in a relatively
   small hit to compression, on the order of 1.5% to 3%.  A mask of 13 bits can
   be used instead if a hit of less than 1% to the compression is desired, at
   the expense of more blocks transmitted for rsync updates.  (Your mileage may
   vary.)

   This implementation of rsyncable uses a different hash algorithm than what
   the gzip rsyncable patch uses in order to provide better performance in
   several regards.  The algorithm is simply to shift the hash value left one
   bit and exclusive-or that with the next byte.  This is masked to the number
   of hash bits (RSYNCMASK) and compared to all ones except for a zero in the
   top bit (RSYNCHIT). This rolling hash has a very small window of 19 bytes
   (RSYNCBITS+7).  The small window provides the benefit of much more rapid
   resynchronization after a change, than does the 4096-byte window of the gzip
   rsyncable patch.

   The comparison value is chosen to avoid matching any repeated bytes or short
   sequences.  The gzip rsyncable patch on the other hand uses a sum and zero
   for comparison, which results in certain bad behaviors, such as always
   matching everywhere in a long sequence of zeros.  Such sequences occur
   frequently in tar files.

   This hash efficiently discards history older than 19 bytes simply by
   shifting that data past the top of the mask -- no history needs to be
   retained to undo its impact on the hash value, as is needed for a sum.

   The choice of the comparison value (RSYNCHIT) has the virtue of avoiding
   extremely short blocks.  The shortest block is five bytes (RSYNCBITS-7) from
   hit to hit, and is unlikely.  Whereas with the gzip rsyncable algorithm,
   blocks of one byte are not only possible, but in fact are the most likely
   block size.

   Thanks and acknowledgement to Kevin Day for his experimentation and insights
   on rsyncable hash characteristics that led to some of the choices here.
 */
#define RSYNCBITS 12
#define RSYNCMASK ((1U << RSYNCBITS) - 1)
#define RSYNCHIT (RSYNCMASK >> 1)

/* initial pool counts and sizes -- INBUFS is the limit on the number of input
   spaces as a function of the number of processors (used to throttle the
   creation of compression jobs), OUTPOOL is the initial size of the output
   data buffer, chosen to make resizing of the buffer very unlikely and to
   allow prepending with a dictionary for use as an input buffer for zopfli */
#define INBUFS(p) (((p)<<1)+3)
#define OUTPOOL(s) ((s)+((s)>>4)+DICT)

/* input buffer size */
#define BUF 32768U

/* globals (modified by main thread only when it's the only thread) */
local struct {
    char *prog;             /* name by which pigz was invoked */
    int outd;               /* output file descriptor */
    char *outf;             /* output file name (allocated if not NULL) */
    int verbosity;          /* 0 = quiet, 1 = normal, 2 = verbose, 3 = trace */
    time_t mtime;           /* time stamp from input file for gzip header */
    int level;              /* compression level */
    int procs;              /* maximum number of compression threads (>= 1) */
    size_t block;           /* uncompressed input size per thread (>= 32K) */
    int warned;             /* true if a warning has been given */
} g;

/* display a complaint with the program name on stderr */
local int complain(const char *fmt, ...)
{
    va_list ap;

    if (g.verbosity > 0) {
        fprintf(stderr, "%s: ", g.prog);
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
        putc_unlocked('\n', stderr);
        fflush(stderr);
        g.warned = 1;
    }
    return 0;
}

/* exit with error, delete output file if in the middle of writing it */
local int bail(const char *why, const char *what)
{
    if (g.outd != -1 && g.outf != NULL)
        unlink(g.outf);
    complain("abort: %s%s", why, what);
    exit(1);
    return 0;
}

/* write len bytes, repeating write() calls as needed */
local void writen(int desc, unsigned char *buf, size_t len)
{
    ssize_t ret;

    while (len) {
        ret = write(desc, buf, len);
        if (ret < 1) {
            complain("write error code %d", errno);
            bail("write error on ", g.outf);
        }
        buf += ret;
        len -= ret;
    }
}

/* put a 4-byte integer into a byte array in LSB order or MSB order */
#define PUT2L(a,b) (*(a)=(b)&0xff,(a)[1]=(b)>>8)
#define PUT4L(a,b) (PUT2L(a,(b)&0xffff),PUT2L((a)+2,(b)>>16))
#define PUT4M(a,b) (*(a)=(b)>>24,(a)[1]=(b)>>16,(a)[2]=(b)>>8,(a)[3]=(b))

/* write a gzip, zlib, or zip header using the information in the globals */
local unsigned long put_header(void)
{
    unsigned long len;
    unsigned char head[30];

    head[0] = 31;
    head[1] = 139;
    head[2] = 8;                /* deflate */
    head[3] = 0;
    PUT4L(head + 4, g.mtime);
    head[8] = g.level >= 9 ? 2 : (g.level == 1 ? 4 : 0);
    head[9] = 3;                /* unix */
    writen(g.outd, head, 10);
    len = 10;
    return len;
}

/* write a gzip, zlib, or zip trailer */
local void put_trailer(unsigned long ulen, unsigned long clen,
                       unsigned long check, unsigned long head)
{
    unsigned char tail[46];

    PUT4L(tail, check);
    PUT4L(tail + 4, ulen);
    writen(g.outd, tail, 8);
}

/* compute check value depending on format */
#define CHECK(a,b,c) crc32(a,b,c)

#ifndef NOTHREAD
/* -- threaded portions of pigz -- */

/* -- check value combination routines for parallel calculation -- */

#define COMB(a,b,c) crc32_comb(a,b,c)
/* combine two crc-32's or two adler-32's (copied from zlib 1.2.3 so that pigz
   can be compatible with older versions of zlib) */

/* we copy the combination routines from zlib here, in order to avoid
   linkage issues with the zlib 1.2.3 builds on Sun, Ubuntu, and others */

local unsigned long gf2_matrix_times(unsigned long *mat, unsigned long vec)
{
    unsigned long sum;

    sum = 0;
    while (vec) {
        if (vec & 1)
            sum ^= *mat;
        vec >>= 1;
        mat++;
    }
    return sum;
}

local void gf2_matrix_square(unsigned long *square, unsigned long *mat)
{
    int n;

    for (n = 0; n < 32; n++)
        square[n] = gf2_matrix_times(mat, mat[n]);
}

local unsigned long crc32_comb(unsigned long crc1, unsigned long crc2,
                               size_t len2)
{
    int n;
    unsigned long row;
    unsigned long even[32];     /* even-power-of-two zeros operator */
    unsigned long odd[32];      /* odd-power-of-two zeros operator */

    /* degenerate case */
    if (len2 == 0)
        return crc1;

    /* put operator for one zero bit in odd */
    odd[0] = 0xedb88320UL;          /* CRC-32 polynomial */
    row = 1;
    for (n = 1; n < 32; n++) {
        odd[n] = row;
        row <<= 1;
    }

    /* put operator for two zero bits in even */
    gf2_matrix_square(even, odd);

    /* put operator for four zero bits in odd */
    gf2_matrix_square(odd, even);

    /* apply len2 zeros to crc1 (first square will put the operator for one
       zero byte, eight zero bits, in even) */
    do {
        /* apply zeros operator for this bit of len2 */
        gf2_matrix_square(even, odd);
        if (len2 & 1)
            crc1 = gf2_matrix_times(even, crc1);
        len2 >>= 1;

        /* if no more bits set, then done */
        if (len2 == 0)
            break;

        /* another iteration of the loop with odd and even swapped */
        gf2_matrix_square(odd, even);
        if (len2 & 1)
            crc1 = gf2_matrix_times(odd, crc1);
        len2 >>= 1;

        /* if no more bits set, then done */
    } while (len2 != 0);

    /* return combined crc */
    crc1 ^= crc2;
    return crc1;
}

#define BASE 65521U     /* largest prime smaller than 65536 */
#define LOW16 0xffff    /* mask lower 16 bits */

/* initialize a pool (pool structure itself provided, not allocated) -- the
   limit is the maximum number of spaces in the pool, or -1 to indicate no
   limit, i.e., to never wait for a buffer to return to the pool */
local void new_pool(struct pool *pool, size_t size, int limit)
{
    pool->have = new_lock(0);
    pool->head = NULL;
    pool->size = size;
    pool->limit = limit;
    pool->made = 0;
}

/* get a space from a pool -- the use count is initially set to one, so there
   is no need to call use_space() for the first use */
local struct space *get_space(struct pool *pool)
{
    struct space *space;

    /* if can't create any more, wait for a space to show up */
    possess(pool->have);
    if (pool->limit == 0)
        wait_for(pool->have, NOT_TO_BE, 0);

    /* if a space is available, pull it from the list and return it */
    if (pool->head != NULL) {
        space = pool->head;
        possess(space->use);
        pool->head = space->next;
        twist(pool->have, BY, -1);      /* one less in pool */
        twist(space->use, TO, 1);       /* initially one user */
        space->len = 0;
        return space;
    }

    /* nothing available, don't want to wait, make a new space */
    assert(pool->limit != 0);
    if (pool->limit > 0)
        pool->limit--;
    pool->made++;
    release(pool->have);
    space = (struct space*)malloc(sizeof(struct space));
    if (space == NULL)
        bail("not enough memory", "");
    space->use = new_lock(1);           /* initially one user */
    space->buf = (unsigned char*)malloc(pool->size);
    if (space->buf == NULL)
        bail("not enough memory", "");
    space->size = pool->size;
    space->len = 0;
    space->pool = pool;                 /* remember the pool this belongs to */
    return space;
}

/* compute next size up by multiplying by about 2**(1/3) and round to the next
   power of 2 if we're close (so three applications results in doubling) -- if
   small, go up to at least 16, if overflow, go to max size_t value */
local size_t grow(size_t size)
{
    size_t was, top;
    int shift;

    was = size;
    size += size >> 2;
    top = size;
    for (shift = 0; top > 7; shift++)
        top >>= 1;
    if (top == 7)
        size = (size_t)1 << (shift + 3);
    if (size < 16)
        size = 16;
    if (size <= was)
        size = (size_t)0 - 1;
    return size;
}

/* increase the size of the buffer in space */
local void grow_space(struct space *space)
{
    size_t more;

    /* compute next size up */
    more = grow(space->size);
    if (more == space->size)
        bail("not enough memory", "");

    /* reallocate the buffer */
    space->buf = (unsigned char*)realloc(space->buf, more);
    if (space->buf == NULL)
        bail("not enough memory", "");
    space->size = more;
}

/* increment the use count to require one more drop before returning this space
   to the pool */
local void use_space(struct space *space)
{
    possess(space->use);
    twist(space->use, BY, +1);
}

/* drop a space, returning it to the pool if the use count is zero */
local void drop_space(struct space *space)
{
    int use;
    struct pool *pool;

    possess(space->use);
    use = peek_lock(space->use);
    assert(use != 0);
    if (use == 1) {
        pool = space->pool;
        possess(pool->have);
        space->next = pool->head;
        pool->head = space;
        twist(pool->have, BY, +1);
    }
    twist(space->use, BY, -1);
}

/* free the memory and lock resources of a pool -- return number of spaces for
   debugging and resource usage measurement */
local int free_pool(struct pool *pool)
{
    int count;
    struct space *space;

    possess(pool->have);
    count = 0;
    while ((space = pool->head) != NULL) {
        pool->head = space->next;
        free(space->buf);
        free_lock(space->use);
        free(space);
        count++;
    }
    assert(count == pool->made);
    release(pool->have);
    free_lock(pool->have);
    return count;
}

/* input and output buffer pools */
local struct pool in_pool;
local struct pool out_pool;
local struct pool dict_pool;
local struct pool lens_pool;

/* -- parallel compression -- */

/* compress or write job (passed from compress list to write list) -- if seq is
   equal to -1, compress_thread is instructed to return; if more is false then
   this is the last chunk, which after writing tells write_thread to return */
struct job {
    long seq;                   /* sequence number */
    int more;                   /* true if this is not the last chunk */
    struct space *in;           /* input data to compress */
    struct space *out;          /* dictionary or resulting compressed data */
    struct space *lens;         /* coded list of flush block lengths */
    unsigned long check;        /* check value for input data */
    lock *calc;                 /* released when check calculation complete */
    struct job *next;           /* next job in the list (either list) */
};

/* list of compress jobs (with tail for appending to list) */
local lock *compress_have = NULL;   /* number of compress jobs waiting */
local struct job *compress_head, **compress_tail;

/* list of write jobs */
local lock *write_first;            /* lowest sequence number in list */
local struct job *write_head;

/* number of compression threads running */
local int cthreads = 0;

/* write thread if running */
local thread *writeth = NULL;

/* setup job lists (call from main thread) */
local void setup_jobs(void)
{
    /* set up only if not already set up*/
    if (compress_have != NULL)
        return;

    /* allocate locks and initialize lists */
    compress_have = new_lock(0);
    compress_head = NULL;
    compress_tail = &compress_head;
    write_first = new_lock(-1);
    write_head = NULL;

    /* initialize buffer pools (initial size for out_pool not critical, since
       buffers will be grown in size if needed -- initial size chosen to make
       this unlikely -- same for lens_pool) */
    new_pool(&in_pool, g.block, INBUFS(g.procs));
    new_pool(&out_pool, OUTPOOL(g.block), -1);
    new_pool(&dict_pool, DICT, -1);
    new_pool(&lens_pool, g.block >> (RSYNCBITS - 1), -1);
}

/* command the compress threads to all return, then join them all (call from
   main thread), free all the thread-related resources */
local void finish_jobs(void)
{
    struct job job;
    int caught;

    /* only do this once */
    if (compress_have == NULL)
        return;

    /* command all of the extant compress threads to return */
    possess(compress_have);
    job.seq = -1;
    job.next = NULL;
    compress_head = &job;
    compress_tail = &(job.next);
    twist(compress_have, BY, +1);       /* will wake them all up */

    /* join all of the compress threads, verify they all came back */
    caught = join_all();
    assert(caught == cthreads);
    cthreads = 0;

    /* free the resources */
    caught = free_pool(&lens_pool);
    caught = free_pool(&dict_pool);
    caught = free_pool(&out_pool);
    caught = free_pool(&in_pool);
    free_lock(write_first);
    free_lock(compress_have);
    compress_have = NULL;
    close(g.outd);
    g.outd = -1;
}

/* compress all strm->avail_in bytes at strm->next_in to out->buf, updating
   out->len, grow the size of the buffer (out->size) if necessary -- respect
   the size limitations of the zlib stream data types (size_t may be larger
   than unsigned) */
local void deflate_engine(z_stream *strm, struct space *out, int flush)
{
    size_t room;

    do {
        room = out->size - out->len;
        if (room == 0) {
            grow_space(out);
            room = out->size - out->len;
        }
        strm->next_out = out->buf + out->len;
        strm->avail_out = room < UINT_MAX ? (unsigned)room : UINT_MAX;
        (void)deflate(strm, flush);
        out->len = strm->next_out - out->buf;
    } while (strm->avail_out == 0);
    assert(strm->avail_in == 0);
}

/* get the next compression job from the head of the list, compress and compute
   the check value on the input, and put a job in the write list with the
   results -- keep looking for more jobs, returning when a job is found with a
   sequence number of -1 (leave that job in the list for other incarnations to
   find) */
local void compress_thread(void *dummy)
{
    struct job *job;                /* job pulled and working on */
    struct job *here, **prior;      /* pointers for inserting in write list */
    unsigned long check;            /* check value of input */
    unsigned char *next;            /* pointer for blocks, check value data */
    size_t left;                    /* input left to process */
    size_t len;                     /* remaining bytes to compress/check */
#if ZLIB_VERNUM >= 0x1260
    int bits;                       /* deflate pending bits */
#endif
    z_stream strm;                  /* deflate stream */

    (void)dummy;

    /* initialize the deflate stream for this thread */
    strm.zfree = Z_NULL;
    strm.zalloc = Z_NULL;
    strm.opaque = Z_NULL;
    if (deflateInit2(&strm, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        bail("not enough memory", "");

    /* keep looking for work */
    for (;;) {
        /* get a job (like I tell my son) */
        possess(compress_have);
        wait_for(compress_have, NOT_TO_BE, 0);
        job = compress_head;
        assert(job != NULL);
        if (job->seq == -1)
            break;
        compress_head = job->next;
        if (job->next == NULL)
            compress_tail = &compress_head;
        twist(compress_have, BY, -1);

        /* got a job -- initialize and set the compression level (note that if
           deflateParams() is called immediately after deflateReset(), there is
           no need to initialize the input/output for the stream) */
	(void)deflateReset(&strm);
	(void)deflateParams(&strm, g.level, Z_DEFAULT_STRATEGY);

        /* set dictionary if provided, release that input or dictionary buffer
           (not NULL if dict is true and if this is not the first work unit) */
        if (job->out != NULL) {
            len = job->out->len;
            left = len < DICT ? len : DICT;
	    deflateSetDictionary(&strm, job->out->buf + (len - left),
				 left);
            drop_space(job->out);
        }

        /* set up input and output */
        job->out = get_space(&out_pool);
	strm.next_in = job->in->buf;
	strm.next_out = job->out->buf;

        /* compress each block, either flushing or finishing */
        next = job->lens == NULL ? NULL : job->lens->buf;
        left = job->in->len;
        job->out->len = 0;
        do {
            /* decode next block length from blocks list */
            len = next == NULL ? 128 : *next++;
            if (len < 128)                          /* 64..32831 */
                len = (len << 8) + (*next++) + 64;
            else if (len == 128)                    /* end of list */
                len = left;
            else if (len < 192)                     /* 1..63 */
                len &= 0x3f;
            else {                                  /* 32832..4227135 */
                len = ((len & 0x3f) << 16) + (*next++ << 8) + 32832U;
                len += *next++;
            }
            left -= len;

	    /* run MAXP2-sized amounts of input through deflate -- this
	       loop is needed for those cases where the unsigned type is
	       smaller than the size_t type, or when len is close to the
	       limit of the size_t type */
	    while (len > MAXP2) {
		strm.avail_in = MAXP2;
		deflate_engine(&strm, job->out, Z_NO_FLUSH);
		len -= MAXP2;
	    }

	    /* run the last piece through deflate -- end on a byte
	       boundary, using a sync marker if necessary, or finish the
	       deflate stream if this is the last block */
	    strm.avail_in = (unsigned)len;
	    if (left || job->more) {
#if ZLIB_VERNUM >= 0x1260
		deflate_engine(&strm, job->out, Z_BLOCK);

		/* add enough empty blocks to get to a byte boundary */
		(void)deflatePending(&strm, Z_NULL, &bits);
		if (bits & 1)
		    deflate_engine(&strm, job->out, Z_SYNC_FLUSH);
		else if (bits & 7) {
		    do {        /* add static empty blocks */
			bits = deflatePrime(&strm, 10, 2);
			assert(bits == Z_OK);
			(void)deflatePending(&strm, Z_NULL, &bits);
		    } while (bits & 7);
		    deflate_engine(&strm, job->out, Z_BLOCK);
		}
#else
		deflate_engine(&strm, job->out, Z_SYNC_FLUSH);
#endif
	    }
	    else
		deflate_engine(&strm, job->out, Z_FINISH);
        } while (left);
        if (job->lens != NULL) {
            drop_space(job->lens);
            job->lens = NULL;
        }

        /* reserve input buffer until check value has been calculated */
        use_space(job->in);

        /* insert write job in list in sorted order, alert write thread */
        possess(write_first);
        prior = &write_head;
        while ((here = *prior) != NULL) {
            if (here->seq > job->seq)
                break;
            prior = &(here->next);
        }
        job->next = here;
        *prior = job;
        twist(write_first, TO, write_head->seq);

        /* calculate the check value in parallel with writing, alert the write
           thread that the calculation is complete, and drop this usage of the
           input buffer */
        len = job->in->len;
        next = job->in->buf;
        check = CHECK(0L, Z_NULL, 0);
        while (len > MAXP2) {
            check = CHECK(check, next, MAXP2);
            len -= MAXP2;
            next += MAXP2;
        }
        check = CHECK(check, next, (unsigned)len);
        drop_space(job->in);
        job->check = check;
        possess(job->calc);
        twist(job->calc, TO, 1);

        /* done with that one -- go find another job */
    }

    /* found job with seq == -1 -- free deflate memory and return to join */
    release(compress_have);
    (void)deflateEnd(&strm);
}

/* collect the write jobs off of the list in sequence order and write out the
   compressed data until the last chunk is written -- also write the header and
   trailer and combine the individual check values of the input buffers */
local void write_thread(void *dummy)
{
    long seq;                       /* next sequence number looking for */
    struct job *job;                /* job pulled and working on */
    size_t len;                     /* input length */
    int more;                       /* true if more chunks to write */
    unsigned long head;             /* header length */
    unsigned long ulen;             /* total uncompressed size (overflow ok) */
    unsigned long clen;             /* total compressed size (overflow ok) */
    unsigned long check;            /* check value of uncompressed data */

    (void)dummy;

    /* build and write header */
    head = put_header();

    /* process output of compress threads until end of input */
    ulen = clen = 0;
    check = CHECK(0L, Z_NULL, 0);
    seq = 0;
    do {
        /* get next write job in order */
        possess(write_first);
        wait_for(write_first, TO_BE, seq);
        job = write_head;
        write_head = job->next;
        twist(write_first, TO, write_head == NULL ? -1 : write_head->seq);

        /* update lengths, save uncompressed length for COMB */
        more = job->more;
        len = job->in->len;
        drop_space(job->in);
        ulen += (unsigned long)len;
        clen += (unsigned long)(job->out->len);

        /* write the compressed data and drop the output buffer */
        writen(g.outd, job->out->buf, job->out->len);
        drop_space(job->out);

        /* wait for check calculation to complete, then combine, once
           the compress thread is done with the input, release it */
        possess(job->calc);
        wait_for(job->calc, TO_BE, 1);
        release(job->calc);
        check = COMB(check, job->check, len);

        /* free the job */
        free_lock(job->calc);
        free(job);

        /* get the next buffer in sequence */
        seq++;
    } while (more);

    /* write trailer */
    put_trailer(ulen, clen, check, head);

    /* verify no more jobs, prepare for next use */
    possess(compress_have);
    assert(compress_head == NULL && peek_lock(compress_have) == 0);
    release(compress_have);
    possess(write_first);
    assert(write_head == NULL);
    twist(write_first, TO, -1);
}

/* compress ind to outd, using multiple threads for the compression and check
   value calculations and one other thread for writing the output -- compress
   threads will be launched and left running (waiting actually) to support
   subsequent calls of parallel_compress() */
void parallel_compress(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, uint32_t(* emitn)(uint32_t, unsigned char*))
{
    // overflow_buf must have size >= PIGZ_BLOCK_SIZE + maximum emission 

    // if overflow_ct is nonzero, this points to the first uncompressed
    // character in overflow_buf
    unsigned char* read_ptr = NULL;

    uint32_t overflow_ct;
    long seq;                       /* sequence number */
    struct space *curr;             /* input data to compress */
    struct space *next;             /* input data that follows curr */
    struct space *dict;             /* dictionary for next compression */
    struct job *job;                /* job for compress, then write */

    int more;                       /* true if more input to read */
    size_t len;                     /* for various length computations */
    uint32_t cur_len;

    g.outf = out_fname;
    g.outd = open(g.outf, O_WRONLY | (do_append? O_APPEND : (O_CREAT | O_TRUNC)), 0644);

    /* if first time or after an option change, setup the job lists */
    setup_jobs();

    /* start write thread */
    writeth = launch(write_thread, NULL);

    /* read from input and start compress threads (write thread will pick up
     the output of the compress threads) */
    seq = 0;
    next = get_space(&in_pool);
    cur_len = emitn(0, overflow_buf);
    if (cur_len > PIGZ_BLOCK_SIZE) {
        memcpy(next->buf, overflow_buf, PIGZ_BLOCK_SIZE);
	next->len = PIGZ_BLOCK_SIZE;
	read_ptr = &(overflow_buf[PIGZ_BLOCK_SIZE]);
    } else {
	memcpy(next->buf, overflow_buf, cur_len);
        next->len = cur_len;
    }
    overflow_ct = cur_len - next->len;

    dict = NULL;
    do {
        /* create a new job */
        job = (struct job*)malloc(sizeof(struct job));
        if (job == NULL)
            bail("not enough memory", "");
        job->calc = new_lock(0);

        /* update input spaces */
        curr = next;

        /* get more input if we don't already have some */
	next = get_space(&in_pool);
	if (overflow_ct >= PIGZ_BLOCK_SIZE) {
	    // no need to call emitn(), since we still have >= 128K of text
	    // from the previous call to compress
	    memcpy(next->buf, read_ptr, PIGZ_BLOCK_SIZE);
	    next->len = PIGZ_BLOCK_SIZE;
	    read_ptr = &(read_ptr[PIGZ_BLOCK_SIZE]);
	    overflow_ct -= PIGZ_BLOCK_SIZE;
	} else {
	    if (overflow_ct) {
	        memcpy(overflow_buf, read_ptr, overflow_ct);
	    }
	    cur_len = emitn(overflow_ct, overflow_buf);
	    if (cur_len > PIGZ_BLOCK_SIZE) {
	        memcpy(next->buf, overflow_buf, PIGZ_BLOCK_SIZE);
		next->len = PIGZ_BLOCK_SIZE;
		read_ptr = &(overflow_buf[PIGZ_BLOCK_SIZE]);
	    } else {
	        memcpy(next->buf, overflow_buf, cur_len);
		next->len = cur_len;
	    }
	    overflow_ct = cur_len - next->len;
	}

        /* if rsyncable, generate block lengths and prepare curr for job to
           likely have less than size bytes (up to the last hash hit) */
        job->lens = NULL;

        /* compress curr->buf to curr->len -- compress thread will drop curr */
        job->in = curr;

        /* set job->more if there is more to compress after curr */
        more = next->len != 0;
        job->more = more;

        /* provide dictionary for this job, prepare dictionary for next job */
        job->out = dict;
        if (more) {
            if (curr->len >= DICT || job->out == NULL) {
                dict = curr;
                use_space(dict);
            }
            else {
                dict = get_space(&dict_pool);
                len = DICT - curr->len;
                memcpy(dict->buf, job->out->buf + (job->out->len - len), len);
                memcpy(dict->buf + len, curr->buf, curr->len);
                dict->len = DICT;
            }
        }

        /* preparation of job is complete */
        job->seq = seq;
        if (++seq < 1)
            bail("input too long: ", "");

        /* start another compress thread if needed */
        if (cthreads < seq && cthreads < g.procs) {
            (void)launch(compress_thread, NULL);
            cthreads++;
        }

        /* put job at end of compress list, let all the compressors know */
        possess(compress_have);
        job->next = NULL;
        *compress_tail = job;
        compress_tail = &(job->next);
        twist(compress_have, BY, +1);
    } while (more);
    drop_space(next);

    /* wait for the write thread to complete (we leave the compress threads out
       there and waiting in case there is another stream to compress) */
    join(writeth);
    writeth = NULL;
    finish_jobs();
}


// about time to implement this without the awkward callback interface...
int32_t pzwrite_init(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr) {
    // unbuffered, and doesn't need to support Windows
    ps_ptr->outd = open(out_fname, O_WRONLY | (do_append? O_APPEND : (O_CREAT | O_TRUNC)), 0644);
    if (ps_ptr->outd == -1) {
        putc_unlocked('\n', stdout);
	fflush(stdout);
        fprintf(stderr, "Error: Failed to open %s.\n", out_fname);
        return 2; // RET_OPEN_FAIL
    }
    ps_ptr->overflow_buf = overflow_buf;
    return 0;
}

void compressed_pzwrite_init(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr) {
    ps_ptr->outd = -1;
    g.outf = out_fname;
    g.outd = open(g.outf, O_WRONLY | (do_append? O_APPEND : (O_CREAT | O_TRUNC)), 0644);

    /* if first time or after an option change, setup the job lists */
    setup_jobs();

    /* start write thread */
    writeth = launch(write_thread, NULL);

    ps_ptr->overflow_buf = overflow_buf;
    ps_ptr->next = NULL;
}

int32_t flex_pzwrite_init(uint32_t output_gz, char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr) {
    if (!output_gz) {
        return pzwrite_init(out_fname, overflow_buf, do_append, ps_ptr);
    } else {
        compressed_pzwrite_init(out_fname, overflow_buf, do_append, ps_ptr);
        return 0;
    }
}

int32_t force_pzwrite(Pigz_state* ps_ptr, char** writep_ptr, uint32_t write_min) {
    unsigned char* writep = (unsigned char*)(*writep_ptr);
    unsigned char* buf = ps_ptr->overflow_buf;
    uint32_t len = (uintptr_t)(writep - buf);
    ssize_t ret;
    while (len) {
        ret = write(ps_ptr->outd, ps_ptr->overflow_buf, len);
	if (ret < 1) {
	    return 6; // RET_WRITE_FAIL
	}
        buf += ret;
        len -= ret;
    }
    *writep_ptr = (char*)(ps_ptr->overflow_buf);
    return 0;
}

void force_compressed_pzwrite(Pigz_state* ps_ptr, char** writep_ptr, uint32_t write_min) {
    // Caller must not request a length-0 write until it's time to close the
    // file.
    unsigned char* writep = (unsigned char*)(*writep_ptr);
    unsigned char* readp = ps_ptr->overflow_buf;
    uint32_t cur_len = (uintptr_t)(writep - readp);

    struct space* curr;             /* input data to compress */
    struct job *job;                /* job for compress, then write */

    int more;                       /* true if more input to read */
    size_t len;                     /* for various length computations */
    if (!ps_ptr->next) {
        ps_ptr->seq = 0;
        ps_ptr->next = get_space(&in_pool);
        if (cur_len > PIGZ_BLOCK_SIZE) {
	    memcpy(ps_ptr->next->buf, readp, PIGZ_BLOCK_SIZE);
            ps_ptr->next->len = PIGZ_BLOCK_SIZE;
            readp = &(readp[PIGZ_BLOCK_SIZE]);
	    cur_len -= PIGZ_BLOCK_SIZE;
	} else {
	    memcpy(ps_ptr->next->buf, readp, cur_len);
            ps_ptr->next->len = cur_len;
	    readp = writep;
	    cur_len = 0;
	}
        ps_ptr->dict = NULL;
	if ((cur_len <= PIGZ_BLOCK_SIZE) && write_min) {
	    // need more input to handle dict properly
	    if (cur_len) {
		memcpy(ps_ptr->overflow_buf, readp, cur_len);
	    }
	    *writep_ptr = (char*)(&(ps_ptr->overflow_buf[cur_len]));
	    return;
	}
    }

    do {
	// create a new job
	job = (struct job*)malloc(sizeof(struct job));
	if (job == NULL) {
	    bail("not enough memory", "");
	}
	job->calc = new_lock(0);
        curr = ps_ptr->next;
	ps_ptr->next = get_space(&in_pool);
        if (cur_len > PIGZ_BLOCK_SIZE) {
	    memcpy(ps_ptr->next->buf, readp, PIGZ_BLOCK_SIZE);
	    ps_ptr->next->len = PIGZ_BLOCK_SIZE;
            readp = &(readp[PIGZ_BLOCK_SIZE]);
	} else {
	    memcpy(ps_ptr->next->buf, readp, cur_len);
	    ps_ptr->next->len = cur_len;
	    readp = writep;
	}
	job->lens = NULL;
	job->in = curr;
	more = (cur_len != 0);
	job->more = more;
        job->out = ps_ptr->dict;
	if (more) {
	    if (curr->len >= DICT || job->out == NULL) {
	        ps_ptr->dict = curr;
	        use_space(ps_ptr->dict);
	    } else {
	        ps_ptr->dict = get_space(&dict_pool);
                len = DICT - curr->len;
                memcpy(ps_ptr->dict->buf, job->out->buf + (job->out->len - len), len);
                memcpy(ps_ptr->dict->buf + len, curr->buf, curr->len);
                ps_ptr->dict->len = DICT;
	    }
	}
	job->seq = ps_ptr->seq;
        if (++(ps_ptr->seq) < 1) {
	    bail("input too long: ", "");
	}
	if (cthreads < ps_ptr->seq && cthreads < g.procs) {
	    (void)launch(compress_thread, NULL);
            cthreads++;
	}
        possess(compress_have);
        job->next = NULL;
        *compress_tail = job;
        compress_tail = &(job->next);
	twist(compress_have, BY, +1);
	cur_len = (uintptr_t)(writep - readp);
    } while ((cur_len >= write_min) && more);
    if (cur_len) {
        memcpy(ps_ptr->overflow_buf, readp, cur_len);
    }
    *writep_ptr = (char*)(&(ps_ptr->overflow_buf[cur_len]));
}

int32_t flex_pzputs_std(Pigz_state* ps_ptr, char** writep_ptr, char* ss, uint32_t sslen) {
    unsigned char* writep = (unsigned char*)(*writep_ptr);
    unsigned char* readp = (unsigned char*)ss;
    uint32_t cur_write_pos = (uintptr_t)(writep - ps_ptr->overflow_buf);
    uint32_t delta;
    int32_t ii;
    while ((sslen + cur_write_pos) > PIGZ_BLOCK_SIZE) {
        if (cur_write_pos <= PIGZ_BLOCK_SIZE) {
	    delta = PIGZ_BLOCK_SIZE + 1 - cur_write_pos;
	    memcpy(writep, readp, delta);
	    writep = &(writep[delta]);
	    readp = &(readp[delta]);
            sslen -= delta;
        }
	if (is_uncompressed_pzwrite(ps_ptr)) {
	    ii = force_pzwrite(ps_ptr, (char**)(&writep), PIGZ_BLOCK_SIZE + 1);
	    if (ii) {
	        return ii;
	    }
	} else {
	    force_compressed_pzwrite(ps_ptr, (char**)(&writep), PIGZ_BLOCK_SIZE + 1);
	}
        cur_write_pos = (uintptr_t)(writep - ps_ptr->overflow_buf);
    }
    memcpy(writep, readp, sslen);
    *writep_ptr = (char*)(&(writep[sslen]));
    return flex_pzwrite(ps_ptr, writep_ptr);
}

int32_t pzwrite_close_null(Pigz_state* ps_ptr, char* writep) {
    int32_t ii = force_pzwrite(ps_ptr, &writep, 0);
    int32_t jj = close(ps_ptr->outd);
    ps_ptr->overflow_buf = NULL;
    return ii || jj;
}

void compressed_pzwrite_close_null(Pigz_state* ps_ptr, char* writep) {
    force_compressed_pzwrite(ps_ptr, &writep, 0);
    drop_space(ps_ptr->next);
    /* wait for the write thread to complete (we leave the compress threads out
       there and waiting in case there is another stream to compress) */
    join(writeth);
    writeth = NULL;
    finish_jobs();
    ps_ptr->overflow_buf = NULL;
}

int32_t flex_pzwrite_close_null(Pigz_state* ps_ptr, char* writep) {
    if (is_uncompressed_pzwrite(ps_ptr)) {
        return pzwrite_close_null(ps_ptr, writep);
    } else {
        compressed_pzwrite_close_null(ps_ptr, writep);
        return 0;
    }
}
#endif

/* catch termination signal */
local void cut_short(int sig)
{
    (void)sig;
    if (g.outd != -1 && g.outf != NULL)
        unlink(g.outf);
    putc_unlocked('\n', stdout);
    _exit(1);
}

/* set option defaults */
void pigz_init(uint32_t setprocs)
{
    signal(SIGINT, cut_short);
    g.level = Z_DEFAULT_COMPRESSION;
#ifdef NOTHREAD
    g.procs = 1;
#else
    g.procs = setprocs;
    // 1023 threads here failed on the NIH test machine
    if (g.procs > 127) {
      g.procs = 127;
    }
#endif
    yarn_prefix = g.prog;
    yarn_abort = cut_short;
    g.block = PIGZ_BLOCK_SIZE;            /* 128K */
    g.verbosity = 1;                /* normal message level */
}
#endif // _WIN32

// provide identical interface for uncompressed writing, to simplify code that
// can generate either compressed or uncompressed output
int32_t write_uncompressed(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, uint32_t(* emitn)(uint32_t, unsigned char*)) {
  uint32_t overflow_ct = 0;
  // if it's potentially worth compressing, it should be text, hence mode "w"
  // instead of "wb"
  // (er, that actually does the wrong thing on Windows.  Fixed in pzwrite.)
  FILE* outfile = fopen(out_fname, do_append? "a" : "w");
  unsigned char* write_ptr;
  uint32_t last_size;
  if (!outfile) {
    putc_unlocked('\n', stdout);
    fflush(stdout);
    fprintf(stderr, "Error: Failed to open %s.\n", out_fname);
    return 2; // RET_OPEN_FAIL
  }
  do {
    last_size = emitn(overflow_ct, overflow_buf);
    if (last_size > PIGZ_BLOCK_SIZE) {
      overflow_ct = last_size - PIGZ_BLOCK_SIZE;
      last_size = PIGZ_BLOCK_SIZE;
    } else {
      overflow_ct = 0;
    }
    if (last_size) {
      if (!fwrite(overflow_buf, last_size, 1, outfile)) {
	putc_unlocked('\n', stdout);
	fflush(stdout);
	fputs("Error: File write failure.\n", stderr);
	fclose(outfile);
	return 6; // RET_WRITE_FAIL
      }
    }
    if (overflow_ct) {
      write_ptr = &(overflow_buf[PIGZ_BLOCK_SIZE]);
      while (overflow_ct > PIGZ_BLOCK_SIZE) {
	if (!fwrite(write_ptr, PIGZ_BLOCK_SIZE, 1, outfile)) {
	  putc_unlocked('\n', stdout);
	  fflush(stdout);
	  fputs("Error: File write failure.\n", stderr);
	  fclose(outfile);
	  return 6;
	}
	write_ptr = &(write_ptr[PIGZ_BLOCK_SIZE]);
	overflow_ct -= PIGZ_BLOCK_SIZE;
      }
      memcpy(overflow_buf, write_ptr, overflow_ct);
    }
  } while (last_size);
  if (fclose(outfile)) {
    putc_unlocked('\n', stdout);
    fflush(stdout);
    fputs("Error: File write failure.\n", stderr);
    return 6;
  }
  return 0;
}
