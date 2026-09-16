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
// which requires libcurl (>= 7.75.0, for its SigV4 request signing) and
// nothing else.
//
// S3 objects and presigned HTTP(S) URLs are streamed on demand via HTTP range
// requests -- no temporary files are written to disk.  The implementation
// lives in s3stream/, which is self-contained and knows nothing about plink2;
// this header is the plink2-facing adapter.  See S3_README.md.

#include <stdio.h>

#include "include/plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Returns 1 if path names a remote object: "s3://", "http://" or "https://".
// This check is always available regardless of whether USE_S3 is defined, so
// that a helpful error message can be printed when a remote path is detected
// but S3 support was not compiled in.
uint32_t IsS3Uri(const char* path);

// Idempotent, thread-safe, and always callable regardless of whether USE_S3
// was compiled in (a no-op in that case).  Intended for hosts without an
// obvious "before main()" hook to call S3Init() from, e.g. the Python
// bindings, where object construction is the earliest place S3 support can
// be wired in.
void EnsureS3Ready();

// Send unsigned requests only, for public buckets (the equivalent of
// `aws s3 --no-sign-request`).  Always callable regardless of whether USE_S3
// was compiled in (a no-op in that case).
void S3SetNoSignRequest(uint32_t no_sign);

// Explicit, per-open S3 credentials, letting different files in the same
// process use completely different accounts/buckets/endpoints -- unlike
// OpenMaybeS3(), which resolves credentials from the ambient environment and
// profile.  Any field left null/0 falls back to the usual default for that
// setting (the standard credential chain if no keys are given and
// no_sign_request is 0).
struct S3Credentials {
  const char* access_key_id;
  const char* secret_access_key;
  const char* session_token;
  const char* endpoint_url;
  const char* region;
  uint32_t no_sign_request;
  uint32_t force_path_style;
};

// Opens a remote object using explicit credentials rather than the ambient
// chain, so it has zero effect on any other file opened in the same process
// (compare to environment variables, which are process-global).  Otherwise
// behaves like OpenMaybeS3(): streams via range requests, returns nullptr on
// failure.  `creds` must be non-null.  Always callable regardless of whether
// USE_S3 was compiled in; without it, prints an error and returns nullptr.
FILE* OpenS3WithCredentials(const char* path, const S3Credentials* creds);

// Opens a remote object for reading and returns a FILE* that streams data on
// demand using range requests; local paths fall through to
// fopen(path, "rb").  Returns nullptr on failure, after printing the reason
// to stderr.
FILE* OpenMaybeS3(const char* path);

#ifdef USE_S3

// Initializes libcurl.  Should be called once before any S3 file is opened,
// and before spawning threads that open S3 files.  EnsureS3Ready() does the
// same thing for hosts without a convenient startup hook.
void S3Init();

// Counterpart to S3Init(), called once at program exit after all S3 FILE*
// handles have been closed.
void S3Shutdown();

#endif  // USE_S3

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_S3_H__
