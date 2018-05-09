#ifndef __PLINK2_CPU_H__
#define __PLINK2_CPU_H__

#ifdef __cplusplus
namespace plink2 {
#endif

// No return value, just prints error and calls exit() if expected CPU features
// are absent (since caller is likely to quickly hit "Illegal instruction",
// which is exactly what we're trying to avoid).
//
// level == 1: sse4.2 expected
// level == 2: avx2 package expected (includes bmi, bmi2, lzcnt)
void VerifyPlink2CpuFeatures(unsigned int level);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_CPU_H__
