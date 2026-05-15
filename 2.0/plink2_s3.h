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

// Initialize the AWS SDK.  Must be called once before any S3 operations,
// and before any threads that use the SDK are spawned.
void S3Init();

// Shut down the AWS SDK.  Must be called once at program exit, after all S3
// operations have completed.
void S3Shutdown();

// Download the S3 object identified by s3_uri (e.g. "s3://bucket/key") to a
// newly-created temporary file.  On success, local_path_buf is filled with
// the path to the temporary file and kPglRetSuccess is returned.
// local_path_buf must be at least kPglFnamesize bytes.
// On failure an appropriate kPglRet* error code is returned.
PglErr S3DownloadToTemp(const char* s3_uri, char* local_path_buf);

#endif  // USE_S3

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_S3_H__
