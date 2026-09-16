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
#include <mutex>

#ifdef USE_S3
#  include <aws/core/Aws.h>
#  include <aws/core/auth/AWSCredentialsProvider.h>
#  include <aws/core/platform/Environment.h>
#  include <aws/core/utils/StringUtils.h>
#  include <aws/core/utils/logging/LogLevel.h>
#  include <aws/s3/S3Client.h>
#  include <aws/s3/S3ClientConfiguration.h>
#  include <aws/s3/S3EndpointProvider.h>
#  include <aws/s3/S3Errors.h>
#  include <aws/s3/model/GetObjectRequest.h>
#  include <algorithm>
#  include <cerrno>
#  include <cinttypes>
#  include <cstdio>
#  include <mutex>
#  include <vector>
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

uint32_t IsS3Uri(const char* path) {
  return (path[0] == 's') && (path[1] == '3') && (path[2] == ':') &&
         (path[3] == '/') && (path[4] == '/');
}

static uint32_t g_s3_no_sign_request = 0;

void S3SetNoSignRequest(uint32_t no_sign) {
  g_s3_no_sign_request = no_sign;
}

void EnsureS3Ready() {
#ifdef USE_S3
  static std::once_flag once;
  std::call_once(once, [] { S3Init(); });
#endif  // USE_S3
}

#ifdef USE_S3

static Aws::SDKOptions g_s3_sdk_options;
static Aws::S3::S3ClientConfiguration* g_s3_config = nullptr;
static Aws::S3::S3Client* g_s3_client = nullptr;
static Aws::S3::S3Client* g_s3_anon_client = nullptr;
static std::mutex g_s3_anon_mutex;

void S3Init() {
  g_s3_sdk_options.loggingOptions.logLevel = Aws::Utils::Logging::LogLevel::Off;
  Aws::InitAPI(g_s3_sdk_options);
  S3InitClientOnly();
}

void S3Shutdown() {
  S3ShutdownClientOnly();
  Aws::ShutdownAPI(g_s3_sdk_options);
}

// The AWS C++ SDK, unlike boto3 and the CLI, has no built-in environment
// variable for the S3 addressing style, so honor both the name used by the
// CLI's s3.addressing_style setting and the one used by the JS/Go SDKs.
static bool S3EnvFlagSet(const char* name) {
  const Aws::String value =
      Aws::Utils::StringUtils::ToLower(Aws::Environment::GetEnv(name).c_str());
  return (value == "1") || (value == "true") || (value == "yes");
}

static bool S3UsePathStyle() {
  if (!Aws::Environment::GetEnv("AWS_S3_FORCE_PATH_STYLE").empty()) {
    return S3EnvFlagSet("AWS_S3_FORCE_PATH_STYLE");
  }
  return Aws::Utils::StringUtils::ToLower(
             Aws::Environment::GetEnv("AWS_S3_ADDRESSING_STYLE").c_str()) ==
         "path";
}

void S3InitClientOnly() {
  // S3ClientConfiguration defaults to virtual-hosted-style addressing, which
  // avoids the legacy path-style URL redirects (HTTP 301) on cross-region
  // requests that caused GetObject/HeadObject to fail.  Path-style is still
  // needed for S3-compatible servers (MinIO, Ceph, LocalStack) and for
  // endpoints addressed by IP.
  g_s3_config = new Aws::S3::S3ClientConfiguration();
  if (S3UsePathStyle()) {
    g_s3_config->useVirtualAddressing = false;
  }
  if (S3EnvFlagSet("AWS_NO_SIGN_REQUEST")) {
    g_s3_no_sign_request = 1;
  }
  g_s3_client = new Aws::S3::S3Client(*g_s3_config);
}

void S3ShutdownClientOnly() {
  delete g_s3_anon_client;
  g_s3_anon_client = nullptr;
  delete g_s3_client;
  g_s3_client = nullptr;
  delete g_s3_config;
  g_s3_config = nullptr;
}

// Shared client for public buckets, created on first use.  It reuses the main
// client's configuration so that the endpoint, region, addressing style and
// timeouts stay consistent.
static Aws::S3::S3Client* S3AnonClient() {
  std::lock_guard<std::mutex> guard(g_s3_anon_mutex);
  if (!g_s3_anon_client) {
    auto anon_creds =
        Aws::MakeShared<Aws::Auth::AnonymousAWSCredentialsProvider>("S3Anon");
    g_s3_anon_client = new Aws::S3::S3Client(
        anon_creds, Aws::MakeShared<Aws::S3::S3EndpointProvider>("S3Anon"),
        *g_s3_config);
  }
  return g_s3_anon_client;
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

  // Points at the shared g_s3_client/g_s3_anon_client by default (not
  // owned).  When this file was opened via OpenS3WithCredentials(), it's a
  // dedicated client exclusive to this file, and owns_client is set so it
  // gets freed on close.
  Aws::S3::S3Client* client;
  bool owns_client;

  // ETag seen when the object was opened.  Every subsequent range request is
  // conditioned on it, so an overwrite mid-read fails loudly instead of
  // silently splicing two versions of the object together.
  Aws::String etag;

  S3FileState() : file_size(-1), pos(0), buf_start(0), buf_end(0), client(nullptr), owns_client(false) {}
  ~S3FileState() {
    if (owns_client) {
      delete client;
    }
  }
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
  if (!state->etag.empty()) {
    req.SetIfMatch(state->etag);
  }

  auto outcome = state->client->GetObject(req);
  if (!outcome.IsSuccess()) {
    const auto& err = outcome.GetError();
    if (err.GetResponseCode() ==
        Aws::Http::HttpResponseCode::PRECONDITION_FAILED) {
      fprintf(stderr,
              "Error: s3://%s/%s was modified while it was being read.\n",
              state->bucket.c_str(), state->key.c_str());
      return false;
    }
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
  if (got == 0) {
    // We are before EOF, so an empty body means the response was truncated.
    // Reporting this as EOF would silently feed plink2 a truncated file.
    fprintf(stderr,
            "Error: S3 range read for s3://%s/%s returned no data at offset %" PRId64
            " (object is %" PRId64 " bytes).\n",
            state->bucket.c_str(), state->key.c_str(), offset, state->file_size);
    return false;
  }
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
// Returns the total file size on success, or -1 on failure, in which case
// *out_err receives the SDK error.  On success *out_etag receives the object's
// ETag, which the caller pins subsequent range reads to.
static int64_t S3ProbeFileSize(Aws::S3::S3Client* client,
                                const Aws::String& bucket,
                                const Aws::String& key,
                                Aws::String* out_etag,
                                Aws::S3::S3Error* out_err) {
  Aws::S3::Model::GetObjectRequest req;
  req.SetBucket(bucket);
  req.SetKey(key);
  req.SetRange("bytes=0-0");

  auto outcome = client->GetObject(req);
  if (!outcome.IsSuccess()) {
    const auto& err = outcome.GetError();
    // "bytes=0-0" is unsatisfiable on a zero-byte object, so a 416 means the
    // object exists and is empty.  There is nothing to pin an ETag to.
    if (err.GetResponseCode() ==
        Aws::Http::HttpResponseCode::REQUESTED_RANGE_NOT_SATISFIABLE) {
      return 0;
    }
    *out_err = err;
    return -1;
  }

  *out_etag = outcome.GetResult().GetETag();

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

// Whether a failed credentialed probe is worth retrying without a signature.
// Throttling, 5xx and network failures have nothing to do with credentials,
// and retrying them anonymously only replaces the real error message with a
// misleading one.
static bool S3ShouldRetryAnonymously(const Aws::S3::S3Error& err) {
  const auto err_type = err.GetErrorType();
  if ((err_type == Aws::S3::S3Errors::NO_SUCH_KEY) ||
      (err_type == Aws::S3::S3Errors::NO_SUCH_BUCKET)) {
    return false;
  }
  if (err.ShouldRetry()) {
    return false;
  }
  const int response_code = static_cast<int>(err.GetResponseCode());
  // 0 means the request never went out, which is what an empty or otherwise
  // unusable credential chain looks like.
  return (response_code == 0) || (response_code == 400) ||
         (response_code == 401) || (response_code == 403);
}

// Open an S3 object as a streaming FILE*.  Uses GetObject(range=0-0) to
// determine the file size (enabling SEEK_END and EOF detection), then creates
// a custom FILE* backed by range requests fetched kS3ChunkSize bytes at a
// time.  Retries without a signature for public buckets when the credentialed
// probe fails for a credential-related reason (the equivalent of
// aws --no-sign-request).
static FILE* S3FileOpenInternal(const char* s3_uri) {
  if (!g_s3_client) {
    fprintf(stderr,
            "Error: S3 support was not initialized before opening %s.\n",
            s3_uri);
    errno = EINVAL;
    return nullptr;
  }
  S3FileState* state = new S3FileState();
  ParseS3Uri(s3_uri, &state->bucket, &state->key);

  // --s3-no-sign-request / AWS_NO_SIGN_REQUEST: never send a signed request,
  // which also avoids a pointless credential lookup on every open.
  if (g_s3_no_sign_request) {
    state->client = S3AnonClient();
    Aws::S3::S3Error anon_err;
    state->file_size = S3ProbeFileSize(state->client, state->bucket, state->key,
                                        &state->etag, &anon_err);
    if (state->file_size < 0) {
      fprintf(stderr, "Error: Cannot access %s: %s\n", s3_uri,
              anon_err.GetMessage().c_str());
      const auto err_type = anon_err.GetErrorType();
      const int not_found = (err_type == Aws::S3::S3Errors::NO_SUCH_KEY) ||
                            (err_type == Aws::S3::S3Errors::NO_SUCH_BUCKET);
      delete state;
      errno = not_found? ENOENT : EACCES;
      return nullptr;
    }
    FILE* fp = S3FileOpenPlatform(state);
    if (!fp) {
      delete state;
    }
    return fp;
  }

  // Try with configured credentials first (env vars, instance profile, etc.).
  state->client = g_s3_client;
  Aws::S3::S3Error first_err;
  state->file_size = S3ProbeFileSize(state->client, state->bucket, state->key,
                                      &state->etag, &first_err);

  if (state->file_size < 0) {
    const auto first_err_type = first_err.GetErrorType();
    const int not_found =
        (first_err_type == Aws::S3::S3Errors::NO_SUCH_KEY) ||
        (first_err_type == Aws::S3::S3Errors::NO_SUCH_BUCKET);
    if (!S3ShouldRetryAnonymously(first_err)) {
      fprintf(stderr, "Error: Cannot access %s: %s\n", s3_uri,
              first_err.GetMessage().c_str());
      delete state;
      errno = not_found? ENOENT : EIO;
      return nullptr;
    }

    // Public datasets (e.g. 1000 Genomes, gnomAD, nf-core test data).
    state->client = S3AnonClient();
    Aws::S3::S3Error anon_err;
    state->file_size = S3ProbeFileSize(state->client, state->bucket,
                                        state->key, &state->etag, &anon_err);
    if (state->file_size < 0) {
      // Report the credentialed error: it is the one the user can act on.
      fprintf(stderr, "Error: Cannot access %s: %s\n", s3_uri,
              first_err.GetMessage().c_str());
      delete state;
      errno = EACCES;
      return nullptr;
    }
  }

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

FILE* OpenS3WithCredentials(const char* s3_uri, const S3Credentials* creds) {
  if (!IsS3Uri(s3_uri)) {
    fprintf(stderr, "Error: %s is not an s3:// URI.\n", s3_uri);
    errno = EINVAL;
    return nullptr;
  }
  if (!creds) {
    fprintf(stderr, "Error: OpenS3WithCredentials() called with null credentials for %s.\n", s3_uri);
    errno = EINVAL;
    return nullptr;
  }

  Aws::S3::S3ClientConfiguration config;
  if (creds->region && creds->region[0]) {
    config.region = creds->region;
  }
  if (creds->endpoint_url && creds->endpoint_url[0]) {
    config.endpointOverride = creds->endpoint_url;
  }
  if (creds->force_path_style) {
    config.useVirtualAddressing = false;
  }

  std::shared_ptr<Aws::Auth::AWSCredentialsProvider> provider;
  if (creds->no_sign_request) {
    provider = Aws::MakeShared<Aws::Auth::AnonymousAWSCredentialsProvider>("S3Explicit");
  } else if (creds->access_key_id && creds->access_key_id[0]) {
    const Aws::Auth::AWSCredentials aws_creds(
        creds->access_key_id,
        creds->secret_access_key ? creds->secret_access_key : "",
        creds->session_token ? creds->session_token : "");
    provider = Aws::MakeShared<Aws::Auth::SimpleAWSCredentialsProvider>("S3Explicit", aws_creds);
  }

  S3FileState* state = new S3FileState();
  state->owns_client = true;
  if (provider) {
    state->client = new Aws::S3::S3Client(
        provider, Aws::MakeShared<Aws::S3::S3EndpointProvider>("S3Explicit"), config);
  } else {
    // No explicit credentials given (e.g. only an endpoint/region override
    // was wanted); fall back to the SDK's default credential provider chain,
    // still with a client dedicated to this file.
    state->client = new Aws::S3::S3Client(config);
  }
  ParseS3Uri(s3_uri, &state->bucket, &state->key);

  Aws::S3::S3Error err;
  state->file_size = S3ProbeFileSize(state->client, state->bucket, state->key,
                                      &state->etag, &err);
  if (state->file_size < 0) {
    fprintf(stderr, "Error: Cannot access %s: %s\n", s3_uri,
            err.GetMessage().c_str());
    const auto err_type = err.GetErrorType();
    const int not_found = (err_type == Aws::S3::S3Errors::NO_SUCH_KEY) ||
                          (err_type == Aws::S3::S3Errors::NO_SUCH_BUCKET);
    delete state;  // also frees the owned client
    errno = not_found? ENOENT : EACCES;
    return nullptr;
  }

  FILE* fp = S3FileOpenPlatform(state);
  if (!fp) {
    delete state;
  }
  return fp;
}

#else  // !USE_S3

FILE* OpenS3WithCredentials(const char* s3_uri, const S3Credentials* creds) {
  (void)creds;
  fprintf(stderr,
          "Error: S3 URI detected (%s) but plink2 was not compiled with S3 support.\n"
          "Rebuild with USE_S3=1 to enable S3 support.\n", s3_uri);
  errno = EINVAL;
  return nullptr;
}

#endif  // USE_S3

#ifdef __cplusplus
}  // namespace plink2
#endif
