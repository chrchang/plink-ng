#include <zstd.h>

#if ZSTD_VERSION_NUMBER < 10502
#  error zstd version too old
#endif
