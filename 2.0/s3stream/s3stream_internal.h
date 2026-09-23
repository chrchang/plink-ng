#ifndef S3STREAM_INTERNAL_H
#define S3STREAM_INTERNAL_H

/* Internals shared between the s3stream translation units.  Not installed;
 * nothing outside this directory should include it. */

#include "s3stream.h"

#include <stdint.h>

#include <map>
#include <string>

#ifdef S3STREAM_ENABLE
#  include <curl/curl.h>
#endif

namespace s3stream {

/* Records a message retrievable via s3stream_last_error().  Thread-local, so
 * concurrent opens do not clobber each other's diagnostics. */
void SetError(const char* fmt, ...);

/* ---------------------------------------------------------------------------
 * Environment / profile helpers
 * ------------------------------------------------------------------------- */

std::string GetEnv(const char* name);

/* True for "1", "true", "yes", "on" (any case). */
bool GetEnvBool(const char* name);

/* ---------------------------------------------------------------------------
 * Validation of server-supplied strings
 * ------------------------------------------------------------------------- */

/* True when `value` is safe to splice into an outgoing header.  Several
 * server-controlled strings -- ETag, IMDS session token, STS session token --
 * are echoed back in later requests, so a CR or LF in one of them would be
 * request-splitting injection into our own traffic.  Rejects empty, control
 * characters, and anything longer than max_len. */
bool IsSafeHeaderValue(const std::string& value, size_t max_len);

/* True for a syntactically plausible AWS region ([a-z0-9-], <= 32 chars).
 * Regions arrive from redirect responses and are spliced into both the
 * request URL and the SigV4 credential scope. */
bool IsSafeRegion(const std::string& region);

/* Best-effort overwrite of a secret before its buffer is released.  Earlier
 * reallocations of the same std::string cannot be reached, so this narrows
 * the window rather than closing it. */
void SecureZero(std::string* secret);

/* Looks a key up in the active profile of ~/.aws/config (AWS_CONFIG_FILE),
 * e.g. "region" or "endpoint_url".  Empty string if absent. */
std::string ConfigValue(const char* key);

/* ---------------------------------------------------------------------------
 * Credentials
 * ------------------------------------------------------------------------- */

struct Credentials {
  std::string access_key;
  std::string secret_key;
  std::string session_token;
  /* Unix time at which these stop working, or 0 when they do not expire.
   * Only container/IMDS credentials carry an expiry. */
  int64_t expires_at;

  Credentials() : expires_at(0) {}
  ~Credentials();

  bool Empty() const { return access_key.empty() || secret_key.empty(); }

  /* Treats credentials as expired slightly early, so a long read does not
   * fail on a token that lapses mid-request. */
  bool Expired() const;
};

/* Walks the credential chain.  Returns false when no source yields
 * credentials, which the caller may treat as "send unsigned". */
bool ResolveCredentials(Credentials* out);

/* ---------------------------------------------------------------------------
 * HTTP layer
 * ------------------------------------------------------------------------- */

struct Response {
  long status;
  std::string body;
  /* Response headers, keys lowercased.  Only the few headers s3stream acts on
   * are retained, so a server cannot grow this map without bound. */
  std::map<std::string, std::string> headers;
  /* Set when Request::max_body aborted the transfer, i.e. the server sent
   * more than was asked for.  Distinct from a short read, and the two have
   * opposite diagnoses. */
  bool body_capped;

  Response() : status(0), body_capped(false) {}

  const std::string* Header(const char* name) const;
};

struct Request {
  std::string url;
  /* Signing inputs; ignored when sign is false. */
  std::string region;
  const Credentials* creds;
  /* Inclusive byte range, e.g. "0-1023".  Empty requests the whole object. */
  std::string range;
  /* When set, the object must still have this ETag or the request fails with
   * 412, which is how a mid-read overwrite is detected. */
  std::string if_match;
  bool sign;
  bool head;
  /* Caps how much body is buffered; the transfer is aborted once exceeded and
   * Response::body_capped is set.  Guards every read against a server that
   * ignores Range and starts streaming the whole object.  0 means no limit,
   * which no caller should use for a response it did not size first. */
  size_t max_body;

  Request() : creds(nullptr), sign(true), head(false), max_body(0) {}
};

/* Issues the request, retrying throttling and transient 5xx responses with
 * exponential backoff plus jitter.  Returns false only when no HTTP response
 * could be obtained at all; an HTTP error status is reported through
 * out->status with false never returned for it. */
bool Perform(const Request& req, Response* out);

#ifdef S3STREAM_ENABLE
/* Limits `curl` to http/https.  Returns false if libcurl refuses, which
 * callers must treat as fatal. */
bool RestrictToHttp(CURL* curl);
#endif

/* Retries `req` against the region named in a 301/400 redirect response, if
 * the response carries one.  `url` is rewritten in place on success.  Returns
 * false when the response was not a region redirect. */
bool RetryInRedirectedRegion(Request* req, Response* out);

/* ---------------------------------------------------------------------------
 * Endpoint construction
 * ------------------------------------------------------------------------- */

/* Splits "s3://bucket/key".  Returns false when the URI is malformed. */
bool ParseS3Uri(const char* uri, std::string* bucket, std::string* key);

/* Percent-encodes an object key for use in a URL path, leaving '/' and the
 * RFC 3986 unreserved characters alone. */
std::string EncodeKey(const std::string& key);

/* Builds the object URL, choosing virtual-hosted or path-style addressing.
 * `force_path_style` forces the latter regardless of the usual heuristics. */
std::string BuildObjectUrl(const std::string& endpoint,
                           const std::string& region,
                           const std::string& bucket, const std::string& key,
                           bool force_path_style);

}  // namespace s3stream

#endif  // S3STREAM_INTERNAL_H
