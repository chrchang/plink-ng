// Mostly copied from libdeflate/lib/x86/cpu_features.c , except that we're
// interested in different features (SSE4.2, AVX2).
//
// It's currently just used to print a more informative error message when e.g.
// the AVX2 build is run on a machine that doesn't support it.  As a
// consequence, it must be compiled with minimal flags (just -O2,
// -DCPU_CHECK_xxx, and warnings).
//
// Original license text follows.

/*
 * x86/cpu_features.c - feature detection for x86 processors
 *
 * Copyright 2016 Eric Biggers
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#if defined(CPU_CHECK_SSE42) || defined(CPU_CHECK_AVX2)

#  include <stdio.h>
#  include <stdlib.h>
#  include <stdint.h>
#  include <inttypes.h>

/* With old GCC versions we have to manually save and restore the x86_32 PIC
 * register (ebx).  See: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=47602  */
#  if defined(__i386__) && defined(__PIC__)
#    define EBX_CONSTRAINT "=r"
#  else
#    define EBX_CONSTRAINT "=b"
#  endif

/* Execute the CPUID instruction.  */
static inline void cpuid(uint32_t leaf, uint32_t subleaf, uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d) {
  __asm__(".ifnc %%ebx, %1; mov  %%ebx, %1; .endif\n"
	  "cpuid                                  \n"
	  ".ifnc %%ebx, %1; xchg %%ebx, %1; .endif\n"
	  : "=a" (*a), EBX_CONSTRAINT (*b), "=c" (*c), "=d" (*d)
	  : "a" (leaf), "c" (subleaf));
}

#  ifdef CPU_CHECK_AVX2
/* Read an extended control register.  */
static inline uint64_t read_xcr(uint32_t index) {
  uint32_t edx, eax;

  /* Execute the "xgetbv" instruction.  Old versions of binutils do not
   * recognize this instruction, so list the raw bytes instead.  */
  __asm__ (".byte 0x0f, 0x01, 0xd0" : "=d" (edx), "=a" (eax) : "c" (index));
#    ifdef __cplusplus
  return (static_cast<uint64_t>(edx) << 32) | eax;
#    else
  return (((uint64_t)edx) << 32) | eax;
#    endif
}
#  endif

extern int RealMain(int argc, char** argv);

#  define unlikely(expr) __builtin_expect(!!(expr), 0)

int main(int argc, char** argv) {
  uint32_t dummy1, dummy2, dummy3, dummy4;
  uint32_t max_function;
  uint32_t features_2;
#  ifdef CPU_CHECK_AVX2
  uint32_t features_3, features_4;
#  endif

  /* Get maximum supported function  */
  cpuid(0, 0, &max_function, &dummy2, &dummy3, &dummy4);
  if (unlikely(max_function < 1)) {
    // er, is this possible given that the 64-bit build executes at all?  well,
    // leave it in for now.
    goto CpuCheck_ret_SSE42_FAIL;
  }

  /* Standard feature flags  */
  cpuid(1, 0, &dummy1, &dummy2, &features_2, &dummy4);

  // bit 20 = SSE4.2
  if (unlikely(!(features_2 & 0x100000))) {
    goto CpuCheck_ret_SSE42_FAIL;
  }

#  ifdef CPU_CHECK_AVX2
  // os_saves_ymm_regs must be true for AVX2
  if (unlikely(
          (!(features_2 & 0x8000000)) ||
          ((read_xcr(0) & 0x6) != 0x6))) {
    goto CpuCheck_ret_AVX2_FAIL;
  }
  if (unlikely(max_function < 7)) {
    goto CpuCheck_ret_AVX2_FAIL;
  }
  cpuid(7, 0, &dummy1, &features_3, &features_4, &dummy4);
  // bit 3 = BMI
  // bit 5 = AVX2
  // bit 8 = BMI2
  if (unlikely((features_3 & 0x128) != 0x128)) {
    goto CpuCheck_ret_AVX2_FAIL;
  }
#  endif

  return RealMain(argc, argv);
 CpuCheck_ret_SSE42_FAIL:
  fputs("Error: This plink2 build requires a processor which supports SSE4.2\ninstructions.  Try a plain 64-bit build instead.\n", stderr);
  exit(13);  // 13 = kPglRetUnsupportedInstructions
#  ifdef CPU_CHECK_AVX2
 CpuCheck_ret_AVX2_FAIL:
  // SSE4.2 doesn't deliver enough of an advantage to justify more clutter on
  // the main downloads page, but technically sophisticated users should be
  // encouraged to build from source in this case.
  fputs("Error: This plink2 build requires a processor which supports AVX2/Haswell\ninstructions, but only SSE4.2 is available.  Try a plain 64-bit build instead,\nor use the build_dynamic/ Makefile to produce a binary that takes advantage of\nSSE4.2 instructions but not AVX2.\n", stderr);
  exit(13);
#  endif
}
#endif  // CPU_CHECK_SSE42 || CPU_CHECK_AVX2
