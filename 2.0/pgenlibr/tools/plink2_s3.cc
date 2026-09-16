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

// plink2-facing adapter over s3stream/.  All the S3 machinery -- credential
// resolution, SigV4 signing, range requests, retries -- lives in s3stream and
// has no plink2 dependencies, so it can be reused elsewhere unchanged.

#include "plink2_s3.h"

#include "s3stream/s3stream.h"

namespace plink2 {

uint32_t IsS3Uri(const char* path) {
  return s3stream_is_remote(path)? 1 : 0;
}

void EnsureS3Ready() {
  s3stream_init();
}

void S3SetNoSignRequest(uint32_t no_sign) {
  s3stream_set_no_sign_request(no_sign? 1 : 0);
}

namespace {

FILE* ReportFailure(const char* path) {
  fprintf(stderr, "Error: %s\n", s3stream_last_error());
  (void)path;
  return nullptr;
}

}  // namespace

FILE* OpenMaybeS3(const char* path) {
  if (!s3stream_is_remote(path)) {
    return fopen(path, FOPEN_RB);
  }
  FILE* result = s3stream_open(path);
  if (!result) {
    return ReportFailure(path);
  }
  return result;
}

FILE* OpenS3WithCredentials(const char* path, const S3Credentials* creds) {
  if (!creds) {
    fprintf(stderr, "Error: OpenS3WithCredentials() called with no credentials.\n");
    return nullptr;
  }
  if (!s3stream_is_remote(path)) {
    return fopen(path, FOPEN_RB);
  }
  // The two structs are kept separate so s3stream stays independent of
  // plink2's headers; they are copied field by field rather than cast.
  s3stream_credentials converted;
  converted.access_key_id = creds->access_key_id;
  converted.secret_access_key = creds->secret_access_key;
  converted.session_token = creds->session_token;
  converted.endpoint_url = creds->endpoint_url;
  converted.region = creds->region;
  converted.no_sign_request = creds->no_sign_request;
  converted.force_path_style = creds->force_path_style;

  FILE* result = s3stream_open_with_credentials(path, &converted);
  if (!result) {
    return ReportFailure(path);
  }
  return result;
}

#ifdef USE_S3

void S3Init() {
  s3stream_init();
}

void S3Shutdown() {
  // s3stream tears libcurl down through its own atexit-free path; there is no
  // global state left to release once every FILE* has been closed.
}

#endif  // USE_S3

}  // namespace plink2
