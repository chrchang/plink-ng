/* HTTP transport for s3stream: request execution, SigV4 signing, retry with
 * backoff, and cross-region redirect recovery. */

#include "s3stream_internal.h"

#ifdef S3STREAM_ENABLE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#  include <windows.h>
#endif

#include <curl/curl.h>

namespace s3stream {

namespace {

/* S3 answers 503 SlowDown and sporadic 500s under normal operation, so a
 * retry loop is required rather than merely nice to have. */
const int kMaxAttempts = 5;
const long kBaseBackoffMs = 200;
const long kMaxBackoffMs = 5000;

const long kConnectTimeoutSeconds = 30;
/* No total-transfer timeout: chunks are large and links are sometimes slow.
 * Stalls are caught by the low-speed limit instead. */
const long kLowSpeedLimitBytes = 64;
const long kLowSpeedTimeSeconds = 120;

/* Header lines longer than this are not worth inspecting; none of the four
 * headers s3stream reads comes close. */
const size_t kMaxHeaderLineBytes = 8192;

size_t WriteToString(char* ptr, size_t size, size_t nmemb, void* userdata) {
  const size_t n = size * nmemb;
  static_cast<std::string*>(userdata)->append(ptr, n);
  return n;
}

struct CappedSink {
  std::string* body;
  size_t limit;
};

/* Returning short of the offered count makes curl abort with
 * CURLE_WRITE_ERROR, which is how the cap stops an oversized response. */
size_t WriteCapped(char* ptr, size_t size, size_t nmemb, void* userdata) {
  CappedSink* sink = static_cast<CappedSink*>(userdata);
  const size_t n = size * nmemb;
  if (sink->body->size() + n > sink->limit) {
    return 0;
  }
  sink->body->append(ptr, n);
  return n;
}

size_t CollectHeader(char* ptr, size_t size, size_t nmemb, void* userdata) {
  const size_t n = size * nmemb;
  /* Only the headers s3stream actually acts on are kept.  curl will happily
   * deliver as many headers as the server cares to send, and storing them all
   * hands a hostile endpoint unbounded memory growth for free. */
  static const char* const kWanted[] = {"content-length", "etag",
                                        "content-range",
                                        "x-amz-bucket-region"};
  if (n > kMaxHeaderLineBytes) {
    return n;
  }
  std::string line(ptr, n);
  const size_t colon = line.find(':');
  if (colon == std::string::npos) {
    return n;
  }
  std::string name = line.substr(0, colon);
  for (size_t i = 0; i < name.size(); ++i) {
    if (name[i] >= 'A' && name[i] <= 'Z') {
      name[i] = static_cast<char>(name[i] - 'A' + 'a');
    }
  }
  bool wanted = false;
  for (size_t i = 0; i < sizeof(kWanted) / sizeof(kWanted[0]); ++i) {
    if (name == kWanted[i]) {
      wanted = true;
      break;
    }
  }
  if (!wanted) {
    return n;
  }
  std::string value = line.substr(colon + 1);
  size_t begin = 0;
  while (begin < value.size() && (value[begin] == ' ' || value[begin] == '\t')) {
    ++begin;
  }
  size_t end = value.size();
  while (end > begin && (value[end - 1] == '\r' || value[end - 1] == '\n' ||
                         value[end - 1] == ' ' || value[end - 1] == '\t')) {
    --end;
  }
  (*static_cast<std::map<std::string, std::string>*>(userdata))[name] =
      value.substr(begin, end - begin);
  return n;
}

bool StatusIsRetryable(long status) {
  return (status == 429) || (status == 500) || (status == 502) ||
         (status == 503) || (status == 504);
}

bool CurlErrorIsRetryable(CURLcode code) {
  switch (code) {
    case CURLE_COULDNT_CONNECT:
    case CURLE_COULDNT_RESOLVE_HOST:
    case CURLE_OPERATION_TIMEDOUT:
    case CURLE_PARTIAL_FILE:
    case CURLE_RECV_ERROR:
    case CURLE_SEND_ERROR:
    case CURLE_GOT_NOTHING:
      return true;
    default:
      return false;
  }
}

/* Full jitter: a uniform draw from [0, capped exponential backoff).  This
 * avoids the synchronized retry storms that plain exponential backoff causes
 * when many readers hit the same throttled prefix.
 *
 * Uses a thread-local xorshift rather than rand(), which is not thread-safe
 * and whose global state belongs to the host program. */
unsigned long NextJitter() {
#if defined(_MSC_VER)
  static __declspec(thread) unsigned long state = 0;
#else
  static __thread unsigned long state = 0;
#endif
  if (state == 0) {
    state = static_cast<unsigned long>(time(nullptr)) ^
            (static_cast<unsigned long>(reinterpret_cast<uintptr_t>(&state)) *
             2654435761UL);
    if (state == 0) {
      state = 1;
    }
  }
  state ^= state << 13;
  state ^= state >> 17;
  state ^= state << 5;
  return state;
}

void SleepBeforeRetry(int attempt) {
  long ceiling = kBaseBackoffMs << attempt;
  if (ceiling > kMaxBackoffMs) {
    ceiling = kMaxBackoffMs;
  }
  const long delay_ms =
      (ceiling > 0) ? static_cast<long>(NextJitter() % static_cast<unsigned long>(ceiling)) : 0;
#ifdef _WIN32
  Sleep(static_cast<DWORD>(delay_ms));
#else
  struct timespec ts;
  ts.tv_sec = delay_ms / 1000;
  ts.tv_nsec = (delay_ms % 1000) * 1000000L;
  nanosleep(&ts, nullptr);
#endif
}

/* S3 reports the real location of a misaddressed bucket in this header,
 * on both the 301 and the 400 AuthorizationHeaderMalformed variants. */
std::string RedirectRegion(const Response& response) {
  const std::string* header = response.Header("x-amz-bucket-region");
  return header ? *header : std::string();
}

bool HasPrefix(const std::string& text, const char* prefix) {
  return text.compare(0, strlen(prefix), prefix) == 0;
}

/* Fails a request before it is sent.  The message goes into error_buffer
 * rather than SetError(), because Perform() reports that buffer. */
CURLcode AbortRequest(CURL* curl, struct curl_slist* headers,
                      char* error_buffer, const char* message) {
  snprintf(error_buffer, CURL_ERROR_SIZE, "%s", message);
  curl_slist_free_all(headers);
  curl_easy_cleanup(curl);
  return CURLE_BAD_FUNCTION_ARGUMENT;
}

CURLcode PerformOnce(const Request& req, Response* out, char* error_buffer) {
  CURL* curl = curl_easy_init();
  if (!curl) {
    return CURLE_FAILED_INIT;
  }
  struct curl_slist* headers = nullptr;

  out->body.clear();
  out->headers.clear();
  out->status = 0;
  out->body_capped = false;

  CappedSink sink;
  sink.body = &out->body;
  sink.limit = req.max_body;

  curl_easy_setopt(curl, CURLOPT_URL, req.url.c_str());
  if (req.max_body) {
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCapped);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &sink);
  } else {
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteToString);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &out->body);
  }
  if (req.head) {
    curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
  }
  curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, CollectHeader);
  curl_easy_setopt(curl, CURLOPT_HEADERDATA, &out->headers);
  curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, error_buffer);
  curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);
  curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, kConnectTimeoutSeconds);
  curl_easy_setopt(curl, CURLOPT_LOW_SPEED_LIMIT, kLowSpeedLimitBytes);
  curl_easy_setopt(curl, CURLOPT_LOW_SPEED_TIME, kLowSpeedTimeSeconds);
  /* A SigV4 signature is bound to the host it was computed for, so following
   * a redirect would send a signature the new host must reject.  Region
   * redirects are handled explicitly instead.
   *
   * The options below are checked because a runtime libcurl older than the
   * headers rejects unknown options at run time rather than at build time,
   * which would silently drop the protection. */
  if (curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 0L) != CURLE_OK) {
    return AbortRequest(curl, headers, error_buffer,
                        "libcurl rejected CURLOPT_FOLLOWLOCATION");
  }
  /* Certificate verification is libcurl's default, but say so explicitly:
   * this is the only thing standing between a signed request and a MITM, and
   * it should not silently depend on how libcurl was built.  Plain-http
   * endpoints are exempt, since a libcurl without TLS may reject these. */
  if (HasPrefix(req.url, "https://") &&
      ((curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L) != CURLE_OK) ||
       (curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L) != CURLE_OK))) {
    return AbortRequest(curl, headers, error_buffer,
                        "libcurl cannot enforce TLS certificate verification");
  }
  if (!RestrictToHttp(curl)) {
    return AbortRequest(curl, headers, error_buffer,
                        "libcurl cannot restrict requests to http/https");
  }

  const std::string ca_bundle = GetEnv("AWS_CA_BUNDLE");
  if (!ca_bundle.empty()) {
    curl_easy_setopt(curl, CURLOPT_CAINFO, ca_bundle.c_str());
  }

  if (!req.range.empty()) {
    curl_easy_setopt(curl, CURLOPT_RANGE, req.range.c_str());
  }
  /* if_match and session_token originate with the server.  Both are validated
   * where they are captured, but this is the single point where they become
   * bytes on the wire, so re-check rather than trust the call chain. */
  if (!req.if_match.empty()) {
    if (!IsSafeHeaderValue(req.if_match, kMaxHeaderLineBytes)) {
      return AbortRequest(
          curl, headers, error_buffer,
          "server supplied an ETag that cannot be sent in a header");
    }
    const std::string header = "If-Match: " + req.if_match;
    headers = curl_slist_append(headers, header.c_str());
  }
  std::string userpwd;
  if (req.sign && req.creds && !req.creds->Empty()) {
    if (!req.creds->session_token.empty()) {
      if (!IsSafeHeaderValue(req.creds->session_token, kMaxHeaderLineBytes)) {
        return AbortRequest(
            curl, headers, error_buffer,
            "session token contains characters that cannot be sent in a header");
      }
      const std::string header =
          "x-amz-security-token: " + req.creds->session_token;
      headers = curl_slist_append(headers, header.c_str());
    }
    /* curl signs every header supplied here, which is what makes the
     * security-token and If-Match headers above acceptable to S3. */
    const std::string sigv4 = "aws:amz:" + req.region + ":s3";
    if (curl_easy_setopt(curl, CURLOPT_AWS_SIGV4, sigv4.c_str()) != CURLE_OK) {
      return AbortRequest(curl, headers, error_buffer,
                          "libcurl was built without AWS SigV4 signing support");
    }
    userpwd = req.creds->access_key + ":" + req.creds->secret_key;
    curl_easy_setopt(curl, CURLOPT_USERPWD, userpwd.c_str());
  }
  if (headers) {
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
  }

  const CURLcode res = curl_easy_perform(curl);
  if ((res == CURLE_OK) || (res == CURLE_WRITE_ERROR)) {
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &out->status);
  }
  if (headers) {
    curl_slist_free_all(headers);
  }
  curl_easy_cleanup(curl);
  SecureZero(&userpwd);
  return res;
}

}  // namespace

bool RestrictToHttp(CURL* curl) {
#if LIBCURL_VERSION_NUM >= 0x075500
  const CURLcode res = curl_easy_setopt(curl, CURLOPT_PROTOCOLS_STR, "http,https");
  if (res != CURLE_UNKNOWN_OPTION) {
    return res == CURLE_OK;
  }
#endif
  /* Deprecated in 7.85 but still honored, and the only spelling that
   * 7.75-7.84 runtimes understand. */
#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
  return curl_easy_setopt(curl, CURLOPT_PROTOCOLS,
                          static_cast<long>(CURLPROTO_HTTP | CURLPROTO_HTTPS)) ==
         CURLE_OK;
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif
}

bool IsSafeHeaderValue(const std::string& value, size_t max_len) {
  if (value.empty() || (value.size() > max_len)) {
    return false;
  }
  for (size_t i = 0; i < value.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(value[i]);
    if ((c < 0x20) || (c == 0x7f)) {
      return false;
    }
  }
  return true;
}

bool IsSafeRegion(const std::string& region) {
  if (region.empty() || (region.size() > 32)) {
    return false;
  }
  for (size_t i = 0; i < region.size(); ++i) {
    const char c = region[i];
    if (!(((c >= 'a') && (c <= 'z')) || ((c >= '0') && (c <= '9')) ||
          (c == '-'))) {
      return false;
    }
  }
  return true;
}

void SecureZero(std::string* secret) {
  if (!secret->empty()) {
    volatile char* p = &(*secret)[0];
    for (size_t i = 0; i < secret->size(); ++i) {
      p[i] = '\0';
    }
  }
  secret->clear();
}

Credentials::~Credentials() {
  SecureZero(&secret_key);
  SecureZero(&session_token);
}

const std::string* Response::Header(const char* name) const {
  const std::map<std::string, std::string>::const_iterator it =
      headers.find(name);
  return (it == headers.end()) ? nullptr : &it->second;
}

bool Perform(const Request& req, Response* out) {
  char error_buffer[CURL_ERROR_SIZE];
  CURLcode res = CURLE_OK;
  for (int attempt = 0; attempt < kMaxAttempts; ++attempt) {
    if (attempt > 0) {
      SleepBeforeRetry(attempt - 1);
    }
    error_buffer[0] = '\0';
    res = PerformOnce(req, out, error_buffer);
    /* The body cap aborts the transfer on purpose; the status line already
     * arrived, so this counts as a usable response -- but the caller has to
     * know the body is incomplete and why. */
    if ((res == CURLE_WRITE_ERROR) && req.max_body && out->status) {
      out->body_capped = true;
      return true;
    }
    if (res != CURLE_OK) {
      if (CurlErrorIsRetryable(res)) {
        continue;
      }
      SetError("%s", error_buffer[0] ? error_buffer : curl_easy_strerror(res));
      return false;
    }
    if (!StatusIsRetryable(out->status)) {
      return true;
    }
  }
  if (res != CURLE_OK) {
    SetError("%s", error_buffer[0] ? error_buffer : curl_easy_strerror(res));
    return false;
  }
  /* Exhausted retries on a retryable status; report it as a normal HTTP
   * outcome and let the caller render the error. */
  return true;
}

bool RetryInRedirectedRegion(Request* req, Response* out) {
  if ((out->status != 301) && (out->status != 400)) {
    return false;
  }
  const std::string region = RedirectRegion(*out);
  if (region.empty() || region == req->region) {
    return false;
  }
  /* The region comes straight out of an error response and is about to be
   * spliced into both the request URL and the SigV4 credential scope, so a
   * value like "evil.example.com/" would redirect the next signed request to
   * a host of the server's choosing. */
  if (!IsSafeRegion(region)) {
    return false;
  }
  /* Both addressing styles spell the host "...s3.<region>.amazonaws.com", so
   * anchoring on the "s3." prefix avoids rewriting a bucket or key that
   * happens to contain the region name. */
  const std::string needle = "s3." + req->region;
  const size_t pos = req->url.find(needle);
  if (pos == std::string::npos) {
    return false;
  }
  req->url.replace(pos, needle.size(), "s3." + region);
  req->region = region;
  return Perform(*req, out);
}

}  // namespace s3stream

#endif  // S3STREAM_ENABLE
