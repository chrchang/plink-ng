// Force-included (via -include) when building the glm accuracy-harness
// binary: makes the --glm writers print 17 significant digits instead of 6, so
// errors far below the normal output precision can be measured.
#ifndef GLM_ACCURACY_FULL_PRECISION_H
#define GLM_ACCURACY_FULL_PRECISION_H

#include <stdio.h>
#include "../../include/plink2_string.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static inline char* GlmAccuracyDtoa(double dxx, char* start) {
  return &(start[snprintf(start, 32, "%.17g", dxx)]);
}

#ifdef __cplusplus
}
#endif

#define dtoa_g(dxx, start) GlmAccuracyDtoa(dxx, start)

#endif
