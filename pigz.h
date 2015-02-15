#ifndef __PIGZ_H__
#define __PIGZ_H__

#include <stdint.h>
#include <inttypes.h>

#define PIGZ_BLOCK_SIZE 131072

#ifndef _WIN32
/* -- pool of spaces for buffer management -- */

/* These routines manage a pool of spaces.  Each pool specifies a fixed size
   buffer to be contained in each space.  Each space has a use count, which
   when decremented to zero returns the space to the pool.  If a space is
   requested from the pool and the pool is empty, a space is immediately
   created unless a specified limit on the number of spaces has been reached.
   Only if the limit is reached will it wait for a space to be returned to the
   pool.  Each space knows what pool it belongs to, so that it can be returned.
 */

#include "yarn.h"

/* a space (one buffer for each space) */
struct space {
    lock *use;              /* use count -- return to pool when zero */
    unsigned char *buf;     /* buffer of size size */
    size_t size;            /* current size of this buffer */
    size_t len;             /* for application usage (initially zero) */
    struct pool *pool;      /* pool to return to */
    struct space *next;     /* for pool linked list */
};

/* pool of spaces (one pool for each type needed) */
struct pool {
    lock *have;             /* unused spaces available, lock for list */
    struct space *head;     /* linked list of available buffers */
    size_t size;            /* size of new buffers in this pool */
    int limit;              /* number of new spaces allowed, or -1 */
    int made;               /* number of buffers made */
};

// Note that this does NOT actually capture anywhere near all of pigz's state;
// there are plenty of global variables that prevent multiple
// parallel_compress2 instances from running concurrently.  It's just the bare
// minimum to remove parallel_compress's callback requirement.
typedef struct {
    unsigned char* overflow_buf;
    long seq;
    struct space* dict;
    struct space* next;
    int outd; // uncompressed writing
} Pigz_state;

static inline uint32_t is_uncompressed_pzwrite(Pigz_state* ps_ptr) {
    return ps_ptr->outd != -1;
}
#else
typedef struct {
    unsigned char* overflow_buf;
    FILE* outfile;
    gzFile gz_outfile;
} Pigz_state;

static inline uint32_t is_uncompressed_pzwrite(Pigz_state* ps_ptr) {
    return (ps_ptr->outfile != NULL);
}
#endif // _WIN32 / NOTHREAD

// This interface is obsolete; compressed_pzwrite/flex_pzwrite is far easier to
// use.
void parallel_compress(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, uint32_t(* emitn)(uint32_t, unsigned char*));


static inline void pzwrite_init_null(Pigz_state* ps_ptr) {
    ps_ptr->overflow_buf = NULL;
}

int32_t pzwrite_init(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr);

void compressed_pzwrite_init(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr);

int32_t flex_pzwrite_init(uint32_t output_gz, char* out_fname, unsigned char* overflow_buf, uint32_t do_append, Pigz_state* ps_ptr);

int32_t force_pzwrite(Pigz_state* ps_ptr, char** writep_ptr, uint32_t write_min);

void force_compressed_pzwrite(Pigz_state* ps_ptr, char** writep_ptr, uint32_t write_min);

static inline int32_t pzwrite(Pigz_state* ps_ptr, char** writep_ptr) {
    if ((uintptr_t)(((unsigned char*)(*writep_ptr)) - ps_ptr->overflow_buf) >= PIGZ_BLOCK_SIZE + 1) {
        return force_pzwrite(ps_ptr, writep_ptr, PIGZ_BLOCK_SIZE + 1);
    }
    return 0;
}

static inline void compressed_pzwrite(Pigz_state* ps_ptr, char** writep_ptr) {
    if ((uintptr_t)(((unsigned char*)(*writep_ptr)) - ps_ptr->overflow_buf) >= PIGZ_BLOCK_SIZE + 1) {
        force_compressed_pzwrite(ps_ptr, writep_ptr, PIGZ_BLOCK_SIZE + 1);
    }
}

static inline int32_t flex_pzwrite(Pigz_state* ps_ptr, char** writep_ptr) {
    if ((uintptr_t)(((unsigned char*)(*writep_ptr)) - ps_ptr->overflow_buf) >= PIGZ_BLOCK_SIZE + 1) {
        if (is_uncompressed_pzwrite(ps_ptr)) {
	    return force_pzwrite(ps_ptr, writep_ptr, PIGZ_BLOCK_SIZE + 1);
        }
	force_compressed_pzwrite(ps_ptr, writep_ptr, PIGZ_BLOCK_SIZE + 1);
    }
    return 0;
}

// Assumes overflow_buf has size 2 * PIGZ_BLOCK_SIZE.
int32_t flex_pzputs_std(Pigz_state* ps_ptr, char** writep_ptr, char* ss, uint32_t sslen);

// designed to write allele codes, which are usually length-1, but could have
// length in the millions.  Assumes overflow_buf has size 2 * PIGZ_BLOCK_SIZE.
static inline int32_t flex_pzputs_allele(Pigz_state* ps_ptr, char** writep_ptr, char* allele_code, uint32_t allele_len) {
    // optimize the common case
    if (allele_len == 1) {
        **writep_ptr = *allele_code;
	*writep_ptr += 1;
	return flex_pzwrite(ps_ptr, writep_ptr);
    }
    return flex_pzputs_std(ps_ptr, writep_ptr, allele_code, allele_len);
}

int32_t pzwrite_close_null(Pigz_state* ps_ptr, char* writep);

void compressed_pzwrite_close_null(Pigz_state* ps_ptr, char* writep);

int32_t flex_pzwrite_close_null(Pigz_state* ps_ptr, char* writep);

static inline void pzwrite_close_cond(Pigz_state* ps_ptr, char* writep) {
    if (ps_ptr->overflow_buf) {
        pzwrite_close_null(ps_ptr, writep);
    }
}

static inline void compressed_pzwrite_close_cond(Pigz_state* ps_ptr, char* writep) {
    if (ps_ptr->overflow_buf) {
        compressed_pzwrite_close_null(ps_ptr, writep);
    }
}

static inline void flex_pzwrite_close_cond(Pigz_state* ps_ptr, char* writep) {
    if (ps_ptr->overflow_buf) {
        if (is_uncompressed_pzwrite(ps_ptr)) {
	    pzwrite_close_null(ps_ptr, writep);
        } else {
	    compressed_pzwrite_close_null(ps_ptr, writep);
	}
    }
}


void pigz_init(uint32_t setprocs);

int32_t write_uncompressed(char* out_fname, unsigned char* overflow_buf, uint32_t do_append, uint32_t(* emitn)(uint32_t, unsigned char*));

#endif // __PIGZ_H__
