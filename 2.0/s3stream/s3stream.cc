/* s3stream: endpoint resolution, object streaming, and the public entry
 * points.  See s3stream.h for the API contract. */

#include "s3stream_internal.h"

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef S3STREAM_ENABLE
#  include <inttypes.h>
#  include <mutex>
#  include <curl/curl.h>
#  if defined(_MSC_VER) && !defined(_SSIZE_T_DEFINED)
#    define _SSIZE_T_DEFINED
typedef ptrdiff_t ssize_t;
#  endif
#endif

namespace s3stream {

namespace {

#if defined(_MSC_VER)
#  define S3STREAM_THREAD_LOCAL __declspec(thread)
#else
#  define S3STREAM_THREAD_LOCAL __thread
#endif

/* A plain char array rather than std::string, because a thread_local with a
 * non-trivial destructor would need per-thread teardown in every host. */
S3STREAM_THREAD_LOCAL char g_error[512] = {0};

bool HasPrefix(const char* text, const char* prefix) {
  return text && (strncmp(text, prefix, strlen(prefix)) == 0);
}

}  // namespace

void SetError(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vsnprintf(g_error, sizeof(g_error), fmt, args);
  va_end(args);
}

}  // namespace s3stream

extern "C" int s3stream_is_s3_uri(const char* path) {
  return s3stream::HasPrefix(path, "s3://") ? 1 : 0;
}

extern "C" int s3stream_is_http_url(const char* path) {
  return (s3stream::HasPrefix(path, "http://") ||
          s3stream::HasPrefix(path, "https://"))
             ? 1
             : 0;
}

extern "C" int s3stream_is_remote(const char* path) {
  return (s3stream_is_s3_uri(path) || s3stream_is_http_url(path)) ? 1 : 0;
}

extern "C" const char* s3stream_last_error(void) { return s3stream::g_error; }

#ifndef S3STREAM_ENABLE

/* Stub build: the classification helpers above still work, so callers can
 * detect a remote path and explain that support was not compiled in. */

extern "C" int s3stream_init(void) { return 0; }

extern "C" void s3stream_set_no_sign_request(int no_sign) { (void)no_sign; }

extern "C" FILE* s3stream_open(const char* path) {
  if (s3stream_is_remote(path)) {
    s3stream::SetError(
        "%s: this build has no S3/HTTP support (rebuild with S3STREAM_ENABLE "
        "and libcurl)",
        path);
    return nullptr;
  }
  return fopen(path, "rb");
}

extern "C" FILE* s3stream_open_with_credentials(
    const char* path, const s3stream_credentials* creds) {
  (void)creds;
  return s3stream_open(path);
}

#else  // S3STREAM_ENABLE

namespace s3stream {

namespace {

/* Sized to amortize per-request latency over large sequential reads without
 * making a small read of a large object expensive. */
const int64_t kChunkSize = 8 * 1024 * 1024;

std::once_flag g_init_once;
CURLcode g_init_result = CURLE_OK;
const char* g_too_old_version = nullptr;
bool g_no_sign_request = false;

/* 7.75.0 introduced CURLOPT_AWS_SIGV4.  Checked at run time as well as
 * build time, since the shared library can be older than the headers. */
const unsigned int kMinCurlVersion = 0x074b00;

void GlobalInit() {
  g_init_result = curl_global_init(CURL_GLOBAL_DEFAULT);
  if (g_init_result != CURLE_OK) {
    return;
  }
  const curl_version_info_data* info = curl_version_info(CURLVERSION_NOW);
  if (!info || (info->version_num < kMinCurlVersion)) {
    g_too_old_version = (info && info->version) ? info->version : "unknown";
  }
}

/* ---------------------------------------------------------------------------
 * Endpoint resolution
 * ------------------------------------------------------------------------- */

std::string ResolveRegion(const char* explicit_region) {
  if (explicit_region && *explicit_region) {
    return explicit_region;
  }
  std::string region = GetEnv("AWS_REGION");
  if (region.empty()) {
    region = GetEnv("AWS_DEFAULT_REGION");
  }
  if (region.empty()) {
    region = ConfigValue("region");
  }
  return region.empty() ? std::string("us-east-1") : region;
}

std::string ResolveEndpoint(const char* explicit_endpoint) {
  if (explicit_endpoint && *explicit_endpoint) {
    return explicit_endpoint;
  }
  std::string endpoint = GetEnv("AWS_ENDPOINT_URL_S3");
  if (endpoint.empty()) {
    endpoint = GetEnv("AWS_ENDPOINT_URL");
  }
  if (endpoint.empty()) {
    endpoint = ConfigValue("endpoint_url");
  }
  return endpoint;
}

/* Virtual-hosted addressing puts the bucket in the hostname, which only works
 * when the name is a valid DNS label.  A dot is the common trap: it makes the
 * name fall outside the *.s3.amazonaws.com wildcard certificate, so TLS
 * verification fails even though the request is otherwise correct. */
bool BucketAllowsVirtualHost(const std::string& bucket) {
  if ((bucket.size() < 3) || (bucket.size() > 63)) {
    return false;
  }
  for (size_t i = 0; i < bucket.size(); ++i) {
    const char c = bucket[i];
    const bool ok = ((c >= 'a') && (c <= 'z')) || ((c >= '0') && (c <= '9')) ||
                    (c == '-');
    if (!ok) {
      return false;
    }
  }
  return true;
}

bool ShouldUsePathStyle(const std::string& bucket, bool forced,
                        bool custom_endpoint) {
  if (forced || GetEnvBool("AWS_S3_FORCE_PATH_STYLE")) {
    return true;
  }
  std::string style = GetEnv("AWS_S3_ADDRESSING_STYLE");
  if (style.empty()) {
    style = ConfigValue("s3_addressing_style");
  }
  if (style == "path") {
    return true;
  }
  if (style == "virtual") {
    return false;
  }
  /* S3-compatible servers (MinIO, Ceph, LocalStack) serve every bucket from
   * one hostname, so a custom endpoint means path-style unless the caller
   * said otherwise. */
  if (custom_endpoint) {
    return true;
  }
  return !BucketAllowsVirtualHost(bucket);
}

/* ---------------------------------------------------------------------------
 * Object streaming
 * ------------------------------------------------------------------------- */

struct StreamState {
  std::string url;
  std::string region;
  Credentials creds;
  /* False for anonymous access and for presigned URLs, which carry their own
   * signature in the query string. */
  bool sign;
  /* True when the caller supplied credentials explicitly, which suppresses
   * both the anonymous fallback and refreshes from the ambient chain. */
  bool explicit_creds;
  /* ETag recorded at open.  Every later range request is conditioned on it,
   * so an overwrite mid-read fails loudly instead of silently splicing two
   * versions of the object together. */
  std::string etag;
  /* Human-readable form of the source, for error messages. */
  std::string display;

  int64_t file_size;
  int64_t pos;

  std::string buf;
  int64_t buf_start;
  int64_t buf_end;

  StreamState()
      : sign(true),
        explicit_creds(false),
        file_size(-1),
        pos(0),
        buf_start(0),
        buf_end(0) {}
};

Request MakeRequest(const StreamState* state) {
  Request req;
  req.url = state->url;
  req.region = state->region;
  req.sign = state->sign;
  if (state->sign && !state->creds.Empty()) {
    req.creds = &state->creds;
  }
  return req;
}

/* Re-reads the credential chain when the current ones have lapsed.  Only
 * meaningful for container/IMDS credentials, which are the ones that expire. */
void RefreshIfExpired(StreamState* state) {
  if (state->explicit_creds || !state->sign || !state->creds.Expired()) {
    return;
  }
  Credentials refreshed;
  if (ResolveCredentials(&refreshed)) {
    state->creds = refreshed;
  }
}

void DescribeHttpFailure(const StreamState* state, const Response& response) {
  const char* hint = "";
  switch (response.status) {
    case 403:
      hint = " (access denied: check credentials and bucket policy)";
      break;
    case 404:
      hint = " (not found)";
      break;
    case 412:
      hint = " (object was modified while it was being read)";
      break;
    case 503:
      hint = " (throttled; retries exhausted)";
      break;
    default:
      break;
  }
  SetError("%s: HTTP %ld%s", state->display.c_str(), response.status, hint);
}

/* Fetches the chunk containing `offset` into the read-ahead buffer. */
bool FetchChunk(StreamState* state, int64_t offset) {
  if ((state->file_size >= 0) && (offset >= state->file_size)) {
    state->buf.clear();
    state->buf_start = offset;
    state->buf_end = offset;
    return true;
  }

  int64_t fetch_end = offset + kChunkSize - 1;
  if ((state->file_size >= 0) && (fetch_end >= state->file_size)) {
    fetch_end = state->file_size - 1;
  }
  const int64_t want = fetch_end - offset + 1;

  RefreshIfExpired(state);

  char range[64];
  snprintf(range, sizeof(range), "%" PRId64 "-%" PRId64, offset, fetch_end);

  Request req = MakeRequest(state);
  req.range = range;
  req.if_match = state->etag;
  /* Exactly the requested range and not a byte more.  Without this a server
   * that ignores Range starts streaming the whole object into memory, and the
   * first sign of trouble is the allocator giving up. */
  req.max_body = static_cast<size_t>(want);

  Response response;
  if (!Perform(req, &response)) {
    return false;
  }
  /* 200 means the server served the whole entity instead of the range.  That
   * is only harmless when the range covered the whole object anyway; at any
   * other offset the bytes would land at the wrong file position, so the data
   * would be silently wrong rather than merely oversized. */
  const bool whole_object = (response.status == 200) && (offset == 0) &&
                            (state->file_size >= 0) &&
                            (want == state->file_size) &&
                            !response.body_capped;
  if ((response.status != 206) && !whole_object) {
    if (response.status == 200) {
      SetError("%s: server ignored the Range header at offset %" PRId64,
               state->display.c_str(), offset);
      return false;
    }
    DescribeHttpFailure(state, response);
    return false;
  }

  if (response.body_capped) {
    SetError("%s: server returned more than the %" PRId64
             " bytes requested at offset %" PRId64,
             state->display.c_str(), want, offset);
    return false;
  }

  const int64_t got = static_cast<int64_t>(response.body.size());
  if ((state->file_size >= 0) && (got != want)) {
    /* The range was clipped to the known object size, so a short body is a
     * truncated transfer, never EOF.  Returning it as EOF would hand the
     * caller a silently short file. */
    SetError("%s: short range response at offset %" PRId64 " (%" PRId64
             " of %" PRId64 " bytes; object is %" PRId64 " bytes)",
             state->display.c_str(), offset, got, want, state->file_size);
    return false;
  }
  if (got == 0) {
    SetError("%s: range request at offset %" PRId64 " returned no data",
             state->display.c_str(), offset);
    return false;
  }

  state->buf.swap(response.body);
  state->buf_start = offset;
  state->buf_end = offset + got;
  return true;
}

ssize_t ReadImpl(StreamState* state, char* out, size_t size) {
  if (size == 0) {
    return 0;
  }
  if ((state->file_size >= 0) && (state->pos >= state->file_size)) {
    return 0;
  }
  if ((state->pos < state->buf_start) || (state->pos >= state->buf_end)) {
    if (!FetchChunk(state, state->pos)) {
      errno = EIO;
      return -1;
    }
    if (state->buf_start == state->buf_end) {
      return 0;
    }
  }
  const size_t offset = static_cast<size_t>(state->pos - state->buf_start);
  const size_t available = static_cast<size_t>(state->buf_end - state->pos);
  const size_t to_copy = (size < available) ? size : available;
  memcpy(out, state->buf.data() + offset, to_copy);
  state->pos += static_cast<int64_t>(to_copy);
  return static_cast<ssize_t>(to_copy);
}

int64_t SeekImpl(StreamState* state, int64_t offset, int whence) {
  int64_t target;
  switch (whence) {
    case SEEK_SET:
      target = offset;
      break;
    case SEEK_CUR:
      target = state->pos + offset;
      break;
    case SEEK_END:
      if (state->file_size < 0) {
        errno = ENOTSUP;
        return -1;
      }
      target = state->file_size + offset;
      break;
    default:
      errno = EINVAL;
      return -1;
  }
  if (target < 0) {
    errno = EINVAL;
    return -1;
  }
  state->pos = target;
  return target;
}

}  // namespace

/* ---------------------------------------------------------------------------
 * FILE* adapters
 * ------------------------------------------------------------------------- */

namespace {

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || \
    defined(__NetBSD__)

int FunopenRead(void* cookie, char* buf, int nbytes) {
  if (nbytes < 0) {
    errno = EINVAL;
    return -1;
  }
  const ssize_t got = ReadImpl(static_cast<StreamState*>(cookie), buf,
                               static_cast<size_t>(nbytes));
  return static_cast<int>(got);
}

fpos_t FunopenSeek(void* cookie, fpos_t offset, int whence) {
  return static_cast<fpos_t>(
      SeekImpl(static_cast<StreamState*>(cookie), offset, whence));
}

int FunopenClose(void* cookie) {
  delete static_cast<StreamState*>(cookie);
  return 0;
}

FILE* MakeStreamFile(StreamState* state) {
  return funopen(state, FunopenRead, nullptr, FunopenSeek, FunopenClose);
}

#elif defined(__linux__) || defined(__GLIBC__)

#ifdef __GLIBC__
typedef off64_t cookie_off_t;
#else
typedef off_t cookie_off_t;
#endif

ssize_t CookieRead(void* cookie, char* buf, size_t size) {
  return ReadImpl(static_cast<StreamState*>(cookie), buf, size);
}

int CookieSeek(void* cookie, cookie_off_t* offset, int whence) {
  const int64_t target =
      SeekImpl(static_cast<StreamState*>(cookie),
               static_cast<int64_t>(*offset), whence);
  if (target < 0) {
    return -1;
  }
  *offset = static_cast<cookie_off_t>(target);
  return 0;
}

int CookieClose(void* cookie) {
  delete static_cast<StreamState*>(cookie);
  return 0;
}

FILE* MakeStreamFile(StreamState* state) {
  cookie_io_functions_t funcs;
  memset(&funcs, 0, sizeof(funcs));
  funcs.read = CookieRead;
  funcs.seek = CookieSeek;
  funcs.close = CookieClose;
  return fopencookie(state, "rb", funcs);
}

#else

/* Windows has neither funopen nor fopencookie, and there is no way to hand
 * MSVCRT/UCRT a FILE* backed by our own read/seek callbacks (its FILE
 * struct is opaque and ABI-locked, unlike glibc/BSD libc). The only
 * workaround would be silently downloading the whole object to a temp file
 * before returning a handle, which is unacceptable: it turns a routine
 * `--pfile s3://...` invocation into an unbounded, unannounced download that
 * can be many times larger than expected, with no indication to the user
 * until disk space or bandwidth runs out. So S3/HTTP streaming is disabled
 * on Windows entirely rather than silently staging arbitrarily large
 * objects to local disk. */
FILE* MakeStreamFile(StreamState* state) {
  SetError(
      "%s: S3/HTTP streaming is not supported on Windows (see S3_README.md)",
      state->display.c_str());
  delete state;
  return nullptr;
}

#endif

}  // namespace

/* ---------------------------------------------------------------------------
 * URL construction
 * ------------------------------------------------------------------------- */

bool ParseS3Uri(const char* uri, std::string* bucket, std::string* key) {
  if (!s3stream_is_s3_uri(uri)) {
    return false;
  }
  const char* rest = uri + 5;
  const char* slash = strchr(rest, '/');
  if (!slash || (slash == rest) || (slash[1] == '\0')) {
    return false;
  }
  bucket->assign(rest, static_cast<size_t>(slash - rest));
  key->assign(slash + 1);
  return true;
}

std::string EncodeKey(const std::string& key) {
  static const char kHex[] = "0123456789ABCDEF";
  std::string out;
  out.reserve(key.size());
  size_t start = 0;
  while (true) {
    const size_t slash = key.find('/', start);
    const size_t len =
        (slash == std::string::npos) ? (key.size() - start) : (slash - start);
    const std::string segment = key.substr(start, len);
    /* curl's URL parser collapses "." and ".." segments and SigV4 forbids
     * CURLOPT_PATH_AS_IS, so those two spellings are escaped to survive it.
     * A dot inside a longer segment is left alone: it is unreserved, and S3
     * canonicalizes it unencoded, so escaping it would break the signature. */
    if (segment == ".") {
      out.append("%2E");
    } else if (segment == "..") {
      out.append("%2E%2E");
    } else {
      for (size_t i = 0; i < segment.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(segment[i]);
        const bool unreserved =
            ((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z')) ||
            ((c >= '0') && (c <= '9')) || (c == '-') || (c == '_') ||
            (c == '.') || (c == '~');
        if (unreserved) {
          out.push_back(static_cast<char>(c));
        } else {
          out.push_back('%');
          out.push_back(kHex[c >> 4]);
          out.push_back(kHex[c & 0x0f]);
        }
      }
    }
    if (slash == std::string::npos) {
      break;
    }
    out.push_back('/');
    start = slash + 1;
  }
  return out;
}

std::string BuildObjectUrl(const std::string& endpoint,
                           const std::string& region,
                           const std::string& bucket, const std::string& key,
                           bool force_path_style) {
  const std::string encoded = EncodeKey(key);
  if (!endpoint.empty()) {
    std::string base = endpoint;
    while (!base.empty() && (base[base.size() - 1] == '/')) {
      base.erase(base.size() - 1);
    }
    if (ShouldUsePathStyle(bucket, force_path_style, true)) {
      return base + "/" + bucket + "/" + encoded;
    }
    /* Splice the bucket in front of the custom endpoint's host. */
    const size_t scheme = base.find("://");
    if (scheme == std::string::npos) {
      return base + "/" + bucket + "/" + encoded;
    }
    return base.substr(0, scheme + 3) + bucket + "." +
           base.substr(scheme + 3) + "/" + encoded;
  }
  if (ShouldUsePathStyle(bucket, force_path_style, false)) {
    return "https://s3." + region + ".amazonaws.com/" + bucket + "/" + encoded;
  }
  return "https://" + bucket + ".s3." + region + ".amazonaws.com/" + encoded;
}

namespace {

/* Parses a whole-string non-negative decimal size.  strtoll alone accepts a
 * leading sign and trailing junk, and a negative result would collide with
 * the file_size < 0 "unknown" sentinel. */
bool ParseObjectSize(const std::string& text, int64_t* out) {
  if (text.empty()) {
    return false;
  }
  errno = 0;
  char* end = nullptr;
  const long long value = strtoll(text.c_str(), &end, 10);
  if (errno || !end || (*end != '\0') || (value < 0)) {
    return false;
  }
  *out = static_cast<int64_t>(value);
  return true;
}

/* An ETag is echoed back in If-Match on every subsequent range request, so a
 * value that cannot legally sit in a header must not be stored. */
void RecordEtag(StreamState* state, const Response& response) {
  const std::string* etag = response.Header("etag");
  if (etag && IsSafeHeaderValue(*etag, 512)) {
    state->etag = *etag;
  }
}

/* Determines the object size and current ETag.
 *
 * A HEAD is tried first because it transfers no body.  Presigned URLs are
 * normally signed for GET only and reject it, so the fallback is a one-byte
 * range request whose Content-Range header carries the total size. */
bool ProbeObject(StreamState* state) {
  RefreshIfExpired(state);

  Request head = MakeRequest(state);
  head.head = true;
  head.max_body = 64 * 1024;
  Response response;
  if (!Perform(head, &response)) {
    return false;
  }
  if (RetryInRedirectedRegion(&head, &response)) {
    state->url = head.url;
    state->region = head.region;
  }

  /* An unauthenticated retry rescues the common case of stale ambient
   * credentials pointed at a public bucket.  Explicit credentials are never
   * silently discarded this way. */
  if ((response.status == 403) && state->sign && !state->explicit_creds) {
    Request anon = head;
    anon.sign = false;
    anon.creds = nullptr;
    Response anon_response;
    if (Perform(anon, &anon_response) && (anon_response.status == 200)) {
      state->sign = false;
      response = anon_response;
    }
  }

  if (response.status == 200) {
    const std::string* length = response.Header("content-length");
    if (length && ParseObjectSize(*length, &state->file_size)) {
      RecordEtag(state, response);
      return true;
    }
  }

  Request probe = MakeRequest(state);
  probe.range = "0-0";
  probe.max_body = 64 * 1024;
  Response probe_response;
  if (!Perform(probe, &probe_response)) {
    return false;
  }
  if (RetryInRedirectedRegion(&probe, &probe_response)) {
    state->url = probe.url;
    state->region = probe.region;
  }

  /* A zero-length object has no satisfiable range, so S3 answers 416.  That
   * is a valid empty file, not an error. */
  if (probe_response.status == 416) {
    state->file_size = 0;
    return true;
  }
  if (probe_response.status == 206) {
    const std::string* range = probe_response.Header("content-range");
    if (range) {
      const size_t slash = range->rfind('/');
      if ((slash != std::string::npos) &&
          ParseObjectSize(range->substr(slash + 1), &state->file_size)) {
        RecordEtag(state, probe_response);
        return true;
      }
    }
    SetError("%s: range response carried no usable Content-Range",
             state->display.c_str());
    return false;
  }
  if (probe_response.status == 200) {
    SetError("%s: server ignored the Range header, so it cannot be streamed",
             state->display.c_str());
    return false;
  }

  DescribeHttpFailure(state, probe_response);
  return false;
}

FILE* OpenStream(StreamState* state) {
  if (!ProbeObject(state)) {
    delete state;
    return nullptr;
  }
  FILE* file = MakeStreamFile(state);
  if (!file) {
    /* MakeStreamFile owns the state once it is handed over, including on the
     * failure paths, so there is nothing to free here. */
    return nullptr;
  }
  return file;
}

FILE* OpenRemote(const char* path, const s3stream_credentials* creds) {
  if (s3stream_init() != 0) {
    return nullptr;
  }

  StreamState* state = new StreamState();
  state->display = path;

  if (s3stream_is_http_url(path)) {
    /* Presigned URLs authenticate through their query string; signing them
     * again would invalidate that signature. */
    state->url = path;
    state->sign = false;
    return OpenStream(state);
  }

  std::string bucket;
  std::string key;
  if (!ParseS3Uri(path, &bucket, &key)) {
    SetError("%s: not a valid S3 URI (expected s3://bucket/key)", path);
    delete state;
    return nullptr;
  }

  const bool no_sign =
      creds ? (creds->no_sign_request != 0) : g_no_sign_request;
  state->region = ResolveRegion(creds ? creds->region : nullptr);
  const std::string endpoint =
      ResolveEndpoint(creds ? creds->endpoint_url : nullptr);
  state->url = BuildObjectUrl(endpoint, state->region, bucket, key,
                              creds && creds->force_path_style);

  if (no_sign) {
    state->sign = false;
    return OpenStream(state);
  }

  if (creds && creds->access_key_id && *creds->access_key_id &&
      creds->secret_access_key && *creds->secret_access_key) {
    state->creds.access_key = creds->access_key_id;
    state->creds.secret_key = creds->secret_access_key;
    if (creds->session_token) {
      state->creds.session_token = creds->session_token;
    }
    state->explicit_creds = true;
  } else if (!ResolveCredentials(&state->creds)) {
    /* Nothing in the chain: fall through unsigned, which is right for public
     * buckets and produces a clear 403 otherwise. */
    state->sign = false;
  }
  return OpenStream(state);
}

}  // namespace

}  // namespace s3stream

extern "C" int s3stream_init(void) {
  std::call_once(s3stream::g_init_once, s3stream::GlobalInit);
  if (s3stream::g_init_result != CURLE_OK) {
    s3stream::SetError("could not initialize libcurl: %s",
                       curl_easy_strerror(s3stream::g_init_result));
    return -1;
  }
  if (s3stream::g_too_old_version) {
    s3stream::SetError(
        "libcurl %s is too old for S3 support (7.75.0 or newer is required)",
        s3stream::g_too_old_version);
    return -1;
  }
  return 0;
}

extern "C" void s3stream_set_no_sign_request(int no_sign) {
  s3stream::g_no_sign_request = (no_sign != 0);
}

extern "C" FILE* s3stream_open(const char* path) {
  if (!path) {
    s3stream::SetError("null path");
    return nullptr;
  }
  if (!s3stream_is_remote(path)) {
    return fopen(path, "rb");
  }
  return s3stream::OpenRemote(path, nullptr);
}

extern "C" FILE* s3stream_open_with_credentials(
    const char* path, const s3stream_credentials* creds) {
  if (!path) {
    s3stream::SetError("null path");
    return nullptr;
  }
  if (!creds) {
    s3stream::SetError("%s: null credentials", path);
    return nullptr;
  }
  if (!s3stream_is_remote(path)) {
    return fopen(path, "rb");
  }
  return s3stream::OpenRemote(path, creds);
}

#endif  // S3STREAM_ENABLE
