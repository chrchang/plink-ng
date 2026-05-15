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
#  include <aws/core/utils/logging/LogLevel.h>
#  include <aws/s3/S3Client.h>
#  include <aws/s3/model/GetObjectRequest.h>
#  include <aws/s3/model/HeadObjectRequest.h>
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
  Aws::Client::ClientConfiguration config;
  g_s3_client = new Aws::S3::S3Client(config);
}

void S3Shutdown() {
  delete g_s3_client;
  g_s3_client = nullptr;
  Aws::ShutdownAPI(g_s3_sdk_options);
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

  S3FileState() : file_size(-1), pos(0), buf_start(0), buf_end(0) {}
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

  auto outcome = g_s3_client->GetObject(req);
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
  delete state;
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

// Open an S3 object as a streaming FILE*.  Uses HeadObject to determine the
// file size (enabling SEEK_END and EOF detection), then creates a custom
// FILE* backed by range requests fetched kS3ChunkSize bytes at a time.
static FILE* S3FileOpenInternal(const char* s3_uri) {
  S3FileState* state = new S3FileState();
  ParseS3Uri(s3_uri, &state->bucket, &state->key);

  // HeadObject to get content-length (required for SEEK_END and EOF).
  Aws::S3::Model::HeadObjectRequest head_req;
  head_req.SetBucket(state->bucket);
  head_req.SetKey(state->key);
  auto head_outcome = g_s3_client->HeadObject(head_req);
  if (!head_outcome.IsSuccess()) {
    const auto& err = head_outcome.GetError();
    fprintf(stderr, "Error: Cannot access %s: %s\n", s3_uri,
            err.GetMessage().c_str());
    delete state;
    errno = ENOENT;
    return nullptr;
  }
  state->file_size =
      static_cast<int64_t>(head_outcome.GetResult().GetContentLength());

  FILE* fp = S3FileOpenPlatform(state);
  if (!fp) {
    delete state;
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
