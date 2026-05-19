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

#include "plink2_s3.h"

#include <string.h>  // strchr

#ifdef USE_S3
#  include <aws/core/Aws.h>
#  include <aws/core/auth/AWSCredentialsProvider.h>
#  include <aws/core/utils/logging/LogLevel.h>
#  include <aws/s3/S3Client.h>
#  include <aws/s3/S3ClientConfiguration.h>
#  include <aws/s3/S3Errors.h>
#  include <aws/s3/model/GetObjectRequest.h>
#  include <algorithm>
#  include <cerrno>
#  include <cinttypes>
#  include <cstdio>
#  include <vector>
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

uint32_t IsS3Uri(const char* path) {
  return (path[0] == 's') && (path[1] == '3') && (path[2] == ':') &&
         (path[3] == '/') && (path[4] == '/');
}

#ifdef USE_S3

static Aws::SDKOptions g_s3_sdk_options;
static Aws::S3::S3Client* g_s3_client = nullptr;

void S3Init() {
  g_s3_sdk_options.loggingOptions.logLevel = Aws::Utils::Logging::LogLevel::Off;
  Aws::InitAPI(g_s3_sdk_options);
  S3InitClientOnly();
}

void S3Shutdown() {
  S3ShutdownClientOnly();
  Aws::ShutdownAPI(g_s3_sdk_options);
}

void S3InitClientOnly() {
  // S3ClientConfiguration defaults to virtual-hosted-style addressing, which
  // avoids the legacy path-style URL redirects (HTTP 301) on cross-region
  // requests that caused GetObject/HeadObject to fail.
  Aws::S3::S3ClientConfiguration config;
  g_s3_client = new Aws::S3::S3Client(config);
}

void S3ShutdownClientOnly() {
  delete g_s3_client;
  g_s3_client = nullptr;
}

// Parse "s3://bucket/key" into bucket and key components.
static void ParseS3Uri(const char* s3_uri, Aws::String* bucket,
                       Aws::String* key) {
  const char* rest = s3_uri + 5;  // skip "s3://"
  const char* slash = strchr(rest, '/');
  if (slash) {
    *bucket = Aws::String(rest, static_cast<size_t>(slash - rest));
    *key = Aws::String(slash + 1);
  } else {
    *bucket = Aws::String(rest);
    *key = Aws::String();
  }
}

// Number of bytes fetched per S3 range request.  8 MiB is a good balance
// between API call overhead and memory usage for large sequential reads.
static const int64_t kS3ChunkSize = 8 * 1024 * 1024;

// Per-FILE* state for an S3-backed stream.
struct S3FileState {
  Aws::String bucket;
  Aws::String key;
  int64_t file_size;   // Content-Length from HeadObject; -1 if unknown
  int64_t pos;         // current read position

  // Read-ahead buffer
  std::vector<uint8_t> buf;
  int64_t buf_start;   // file offset of buf[0]
  int64_t buf_end;     // file offset past last valid byte in buf

  // Non-null when anonymous credentials are needed for this file.
  Aws::S3::S3Client* anon_client;

  S3FileState() : file_size(-1), pos(0), buf_start(0), buf_end(0), anon_client(nullptr) {}
  ~S3FileState() { delete anon_client; }
};

// Fetch a chunk of data starting at `offset` via an S3 range request.
// Returns true on success (including EOF where buf becomes empty).
static bool S3FetchChunk(S3FileState* state, int64_t offset) {
  if (state->file_size >= 0 && offset >= state->file_size) {
    // EOF: empty the buffer
    state->buf.clear();
    state->buf_start = offset;
    state->buf_end = offset;
    return true;
  }

  int64_t fetch_end = offset + kS3ChunkSize - 1;
  if (state->file_size >= 0 && fetch_end >= state->file_size) {
    fetch_end = state->file_size - 1;
  }

  char range_str[64];
  snprintf(range_str, sizeof(range_str),
           "bytes=%" PRId64 "-%" PRId64, offset, fetch_end);

  Aws::S3::Model::GetObjectRequest req;
  req.SetBucket(state->bucket);
  req.SetKey(state->key);
  req.SetRange(range_str);

  Aws::S3::S3Client* client = state->anon_client ? state->anon_client : g_s3_client;
  auto outcome = client->GetObject(req);
  if (!outcome.IsSuccess()) {
    const auto& err = outcome.GetError();
    fprintf(stderr, "Error: S3 range read failed for s3://%s/%s: %s\n",
            state->bucket.c_str(), state->key.c_str(),
            err.GetMessage().c_str());
    return false;
  }

  int64_t expected = fetch_end - offset + 1;
  state->buf.resize(static_cast<size_t>(expected));
  auto& body = outcome.GetResult().GetBody();
  body.read(reinterpret_cast<char*>(state->buf.data()), expected);
  int64_t got = static_cast<int64_t>(body.gcount());
  state->buf.resize(static_cast<size_t>(got));
  state->buf_start = offset;
  state->buf_end = offset + got;
  return true;
}

// Core read implementation shared by both platform backends.
static ssize_t S3FileReadImpl(S3FileState* state, char* buf, size_t size) {
  if (size == 0) return 0;
  if (state->file_size >= 0 && state->pos >= state->file_size) {
    return 0;  // EOF
  }

  // If position is outside the buffer, fetch a new chunk.
  if (state->pos < state->buf_start || state->pos >= state->buf_end) {
    if (!S3FetchChunk(state, state->pos)) {
      errno = EIO;
      return -1;
    }
    if (state->buf_start == state->buf_end) {
      return 0;  // EOF
    }
  }

  size_t buf_offset = static_cast<size_t>(state->pos - state->buf_start);
  size_t available = static_cast<size_t>(state->buf_end - state->pos);
  size_t to_copy = size < available ? size : available;
  memcpy(buf, state->buf.data() + buf_offset, to_copy);
  state->pos += static_cast<int64_t>(to_copy);
  return static_cast<ssize_t>(to_copy);
}

// Core seek implementation shared by both platform backends.
// Returns the new position on success, -1 on error.
static int64_t S3FileSeekImpl(S3FileState* state, int64_t offset, int whence) {
  int64_t new_pos;
  switch (whence) {
    case SEEK_SET:
      new_pos = offset;
      break;
    case SEEK_CUR:
      new_pos = state->pos + offset;
      break;
    case SEEK_END:
      if (state->file_size < 0) {
        errno = ENOTSUP;
        return -1;
      }
      new_pos = state->file_size + offset;
      break;
    default:
      errno = EINVAL;
      return -1;
  }
  if (new_pos < 0) {
    errno = EINVAL;
    return -1;
  }
  state->pos = new_pos;
  return new_pos;
}

static int S3FileCloseImpl(S3FileState* state) {
  delete state;  // destructor frees anon_client
  return 0;
}

// ---- Platform-specific FILE* wrappers ----

#ifdef __linux__

static ssize_t S3FileRead(void* cookie, char* buf, size_t size) {
  return S3FileReadImpl(static_cast<S3FileState*>(cookie), buf, size);
}

static int S3FileSeek(void* cookie, off64_t* offset, int whence) {
  int64_t result = S3FileSeekImpl(
      static_cast<S3FileState*>(cookie),
      static_cast<int64_t>(*offset), whence);
  if (result < 0) return -1;
  *offset = static_cast<off64_t>(result);
  return 0;
}

static int S3FileClose(void* cookie) {
  return S3FileCloseImpl(static_cast<S3FileState*>(cookie));
}

static FILE* S3FileOpenPlatform(S3FileState* state) {
  static const cookie_io_functions_t kS3IoFuncs = {
    S3FileRead, nullptr, S3FileSeek, S3FileClose
  };
  return fopencookie(state, "rb", kS3IoFuncs);
}

#else  // macOS / BSD: use funopen

static int S3FileReadFunopen(void* cookie, char* buf, int size) {
  if (size <= 0) return 0;
  ssize_t r = S3FileReadImpl(static_cast<S3FileState*>(cookie), buf,
                              static_cast<size_t>(size));
  return static_cast<int>(r);
}

static fpos_t S3FileSeekFunopen(void* cookie, fpos_t offset, int whence) {
  int64_t result = S3FileSeekImpl(static_cast<S3FileState*>(cookie),
                                   static_cast<int64_t>(offset), whence);
  return static_cast<fpos_t>(result);
}

static int S3FileCloseFunopen(void* cookie) {
  return S3FileCloseImpl(static_cast<S3FileState*>(cookie));
}

static FILE* S3FileOpenPlatform(S3FileState* state) {
  return funopen(state, S3FileReadFunopen, nullptr, S3FileSeekFunopen,
                 S3FileCloseFunopen);
}

#endif  // __linux__

// Probe an S3 object to get its total size by fetching Range: bytes=0-0.
// Uses GetObject instead of HeadObject: cross-region 301 responses carry
// their redirect XML in the body (GET), whereas HEAD 301 has no body and
// the SDK cannot parse the error.
// Returns the total file size on success, or -1 on failure.
// On failure, *out_is_not_found is set to true only for 404-class errors
// (object/bucket does not exist); *out_error_msg receives the SDK message.
static int64_t S3ProbeFileSize(Aws::S3::S3Client* client,
                                const Aws::String& bucket,
                                const Aws::String& key,
                                bool* out_is_not_found,
                                Aws::String* out_error_msg) {
  if (out_is_not_found) *out_is_not_found = false;
  if (out_error_msg) *out_error_msg = "";

  Aws::S3::Model::GetObjectRequest req;
  req.SetBucket(bucket);
  req.SetKey(key);
  req.SetRange("bytes=0-0");

  auto outcome = client->GetObject(req);
  if (!outcome.IsSuccess()) {
    const auto& err = outcome.GetError();
    if (out_error_msg) {
      *out_error_msg = err.GetMessage();
    }
    if (out_is_not_found) {
      const auto err_type = err.GetErrorType();
      *out_is_not_found = (err_type == Aws::S3::S3Errors::NO_SUCH_KEY ||
                           err_type == Aws::S3::S3Errors::NO_SUCH_BUCKET);
    }
    return -1;
  }

  // Parse total size from Content-Range: "bytes 0-0/<total>"
  const auto& content_range = outcome.GetResult().GetContentRange();
  if (!content_range.empty()) {
    const char* slash = strchr(content_range.c_str(), '/');
    if (slash && slash[1] != '*') {
      char* end;
      int64_t total = static_cast<int64_t>(strtoll(slash + 1, &end, 10));
      if (end != slash + 1) {
        return total;
      }
    }
  }
  return -1;
}

// Open an S3 object as a streaming FILE*.  Uses GetObject(range=0-0) to
// determine the file size (enabling SEEK_END and EOF detection), then creates
// a custom FILE* backed by range requests fetched kS3ChunkSize bytes at a
// time.  Falls back to anonymous credentials for public buckets when the
// credentialed probe fails with an auth/other error (but not a 404).
static FILE* S3FileOpenInternal(const char* s3_uri) {
  S3FileState* state = new S3FileState();
  ParseS3Uri(s3_uri, &state->bucket, &state->key);

  // Try with configured credentials first (env vars, instance profile, etc.).
  bool is_not_found = false;
  Aws::String error_msg;
  state->file_size = S3ProbeFileSize(g_s3_client, state->bucket, state->key,
                                      &is_not_found, &error_msg);

  if (state->file_size < 0) {
    if (is_not_found) {
      // 404: the object/bucket does not exist; anonymous access won't help.
      fprintf(stderr, "Error: Cannot access %s: %s\n", s3_uri,
              error_msg.c_str());
      delete state;
      errno = ENOENT;
      return nullptr;
    }

    // Auth or other error — fall back to anonymous credentials for public
    // datasets (e.g. 1000 Genomes, gnomAD, nf-core test data).
    // This is the equivalent of aws --no-sign-request.
    auto anon_creds = Aws::MakeShared<Aws::Auth::AnonymousAWSCredentialsProvider>("S3Anon");
    state->anon_client = new Aws::S3::S3Client(anon_creds, nullptr);
    state->file_size = S3ProbeFileSize(state->anon_client, state->bucket,
                                        state->key, nullptr, &error_msg);
    if (state->file_size < 0) {
      fprintf(stderr, "Error: Cannot access %s: %s\n", s3_uri,
              error_msg.empty() ? "cannot determine file size"
                                : error_msg.c_str());
      delete state;  // destructor frees anon_client
      errno = ENOENT;
      return nullptr;
    }
  }

  FILE* fp = S3FileOpenPlatform(state);
  if (!fp) {
    delete state;  // destructor frees anon_client
  }
  return fp;
}




FILE* OpenMaybeS3(const char* path) {
  if (IsS3Uri(path)) {
    return S3FileOpenInternal(path);
  }
  return fopen(path, FOPEN_RB);
}

#endif  // USE_S3

#ifdef __cplusplus
}  // namespace plink2
#endif
