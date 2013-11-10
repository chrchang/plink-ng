#ifndef __PIGZ_H__
#define __PIGZ_H__

#include <stdint.h>
#include <inttypes.h>

#define PIGZ_BLOCK_SIZE 131072

void parallel_compress(char* out_fname, uint32_t do_append, uint32_t(* emitn)(uint32_t, unsigned char*));

void pigz_init(uint32_t setprocs);

int32_t write_uncompressed(char* out_fname, uint32_t do_append, uint32_t(* emitn)(uint32_t, unsigned char*));

#endif // __PIGZ_H__
