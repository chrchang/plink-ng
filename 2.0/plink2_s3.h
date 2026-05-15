#ifndef __PLINK2_S3_H__
#define __PLINK2_S3_H__

// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Optional S3 support for PLINK 2.0.  Enable by building with USE_S3=1,
// which requires the AWS C++ SDK (aws-sdk-cpp) with the s3 and core
// components installed.
//
// When USE_S3 is enabled, S3 objects are streamed on-demand via HTTP range
// requests — no temporary files are written to disk.

#include <stdio.h>

#include "include/plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Returns 1 if path begins with "s3://", 0 otherwise.
// This check is always available regardless of whether USE_S3 is defined,
// so that a helpful error message can be printed when S3 URIs are detected
// but S3 support was not compiled in.
uint32_t IsS3Uri(const char* path);

#ifdef USE_S3

// Initialize the AWS SDK.  Must be called once before any S3 file opens,
// and before any threads that use the SDK are spawned.
void S3Init();

// Shut down the AWS SDK.  Must be called once at program exit, after all S3
// FILE* handles have been closed.
void S3Shutdown();

// Open an S3 object for reading and return a FILE* that streams data on
// demand using S3 range requests.  Behaves like fopen(path, "rb") but the
// data is fetched from S3 in chunks rather than from a local file.
//
// Returns nullptr on failure (e.g. object not found, access denied, network
// error).  An error message is printed to stderr before returning nullptr.
FILE* OpenMaybeS3(const char* path);

#else

// Without USE_S3, OpenMaybeS3 is a thin wrapper around fopen, allowing
// callers in plink2_text.cc / pgenlib_read.cc to compile unchanged
// regardless of whether USE_S3 is set.
static inline FILE* OpenMaybeS3(const char* path) {
  return fopen(path, FOPEN_RB);
}

#endif  // USE_S3

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_S3_H__
