#ifndef S3STREAM_H
#define S3STREAM_H

/* s3stream -- read-only streaming access to S3 objects and HTTP(S) URLs,
 * exposed as an ordinary FILE*.
 *
 * Reads are served on demand through HTTP range requests, so opening a
 * multi-gigabyte object costs one small request and no local disk.  The
 * returned FILE* supports fread/fseek/ftell/feof and is closed with fclose.
 *
 * The only dependency is libcurl (>= 7.75.0, for CURLOPT_AWS_SIGV4, which
 * performs the SigV4 request signing).  Build with -DS3STREAM_ENABLE to
 * compile the remote backends in; without it the entry points remain
 * callable but reject remote URLs with an error, so callers need no #ifdefs.
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Explicit per-open credentials, letting different files in one process use
 * completely different accounts, buckets and endpoints.  Any field left NULL
 * or 0 falls back to the ambient default for that setting: the credential
 * chain for the keys, AWS_REGION/AWS_ENDPOINT_URL for the endpoint, and
 * automatic addressing-style selection.
 *
 * Lifetime: the struct and the strings it points at are only read during the
 * s3stream_open_with_credentials() call and need not outlive it. */
typedef struct s3stream_credentials {
  const char* access_key_id;
  const char* secret_access_key;
  const char* session_token;
  const char* endpoint_url;
  const char* region;
  unsigned int no_sign_request;
  unsigned int force_path_style;
} s3stream_credentials;

/* Path classification.  Always available, including in builds without
 * S3STREAM_ENABLE, so a caller can emit a helpful "rebuild with S3 support"
 * message instead of a confusing file-not-found. */
int s3stream_is_s3_uri(const char* path);   /* "s3://..."            */
int s3stream_is_http_url(const char* path); /* "http://", "https://" */
int s3stream_is_remote(const char* path);   /* either of the above   */

/* Initializes libcurl.  Idempotent and thread-safe; returns 0 on success.
 * Calling this explicitly is optional -- the open functions do it for you --
 * but hosts that spawn threads should call it once up front, since the
 * underlying curl_global_init() is not itself thread-safe. */
int s3stream_init(void);

/* Process-wide default: send unsigned requests, for public buckets.  This is
 * the equivalent of `aws s3 --no-sign-request`.  Per-open credentials
 * override it via s3stream_credentials::no_sign_request. */
void s3stream_set_no_sign_request(int no_sign);

/* Opens a path for reading.  s3:// URIs and http(s):// URLs are streamed;
 * anything else is passed to fopen(path, "rb"), so this is a drop-in
 * replacement for fopen in read paths.
 *
 * S3 credentials come from the standard chain, in order:
 *   1. AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY / AWS_SESSION_TOKEN
 *   2. AWS_SHARED_CREDENTIALS_FILE, else ~/.aws/credentials
 *   3. ~/.aws/config (or AWS_CONFIG_FILE)
 *   4. ECS/EKS container endpoint (AWS_CONTAINER_CREDENTIALS_RELATIVE_URI
 *      or _FULL_URI)
 *   5. EC2 instance metadata (IMDSv2)
 *   6. unsigned, if nothing above yields credentials
 * The profile is selected by AWS_PROFILE (default "default").
 *
 * Returns NULL on failure; s3stream_last_error() explains why. */
FILE* s3stream_open(const char* path);

/* As s3stream_open(), but uses the supplied credentials instead of the
 * ambient chain, with no effect on any other file in the process.  `creds`
 * must be non-NULL; fields left NULL/0 fall back to ambient defaults. */
FILE* s3stream_open_with_credentials(const char* path,
                                     const s3stream_credentials* creds);

/* Describes the most recent failure on the calling thread.  Never NULL;
 * returns "" when nothing has failed yet.  The buffer is owned by s3stream
 * and is overwritten by the next failing call on the same thread. */
const char* s3stream_last_error(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* S3STREAM_H */
