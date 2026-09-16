/* Credential resolution for s3stream: environment, shared config files, ECS
 * container endpoint, EC2 instance metadata.
 *
 * Deliberately unsupported, because each needs a substantial amount of
 * machinery for a case that has a one-line workaround (run the AWS CLI and
 * export the result):
 *   - SSO (`aws sso login`)
 *   - role_arn / source_profile AssumeRole chaining
 *   - credential_process
 *   - web identity / IRSA
 */

#include "s3stream_internal.h"

#ifdef S3STREAM_ENABLE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <curl/curl.h>

namespace s3stream {

namespace {

/* Refresh this long before the advertised expiry, so credentials do not
 * lapse in the middle of a long read. */
const int64_t kExpirySkewSeconds = 300;

/* The link-local metadata endpoints are unreachable off EC2/ECS, so keep the
 * probe short enough that it does not stall an otherwise fine failure path. */
const long kMetadataTimeoutSeconds = 2;

size_t CollectToString(char* ptr, size_t size, size_t nmemb, void* userdata) {
  const size_t n = size * nmemb;
  static_cast<std::string*>(userdata)->append(ptr, n);
  return n;
}

/* Fetches a metadata URL.  Returns false on any transport or HTTP error. */
bool FetchMetadata(const std::string& url, const char* method,
                   const std::string& header, std::string* body) {
  CURL* curl = curl_easy_init();
  if (!curl) {
    return false;
  }
  struct curl_slist* headers = nullptr;
  if (!header.empty()) {
    headers = curl_slist_append(nullptr, header.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
  }
  if (method) {
    curl_easy_setopt(curl, CURLOPT_CUSTOMREQUEST, method);
  }
  curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, CollectToString);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, body);
  curl_easy_setopt(curl, CURLOPT_TIMEOUT, kMetadataTimeoutSeconds);
  curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
  /* Metadata lives on a link-local address; a proxy would break it. */
  curl_easy_setopt(curl, CURLOPT_NOPROXY, "*");
  const CURLcode res = curl_easy_perform(curl);
  if (headers) {
    curl_slist_free_all(headers);
  }
  curl_easy_cleanup(curl);
  return (res == CURLE_OK) && !body->empty();
}

/* Pulls a string value out of a flat JSON object.  The metadata endpoints
 * return four known fields, which does not justify a JSON dependency. */
std::string JsonString(const std::string& json, const char* key) {
  const std::string needle = std::string("\"") + key + "\"";
  size_t pos = json.find(needle);
  if (pos == std::string::npos) {
    return std::string();
  }
  pos = json.find(':', pos + needle.size());
  if (pos == std::string::npos) {
    return std::string();
  }
  pos = json.find('"', pos);
  if (pos == std::string::npos) {
    return std::string();
  }
  ++pos;
  const size_t end = json.find('"', pos);
  if (end == std::string::npos) {
    return std::string();
  }
  return json.substr(pos, end - pos);
}

/* Parses an ISO-8601 instant such as "2026-09-16T10:32:45Z" into Unix time.
 * Returns 0 when the field is missing or unparseable, which the caller reads
 * as "no expiry" -- the safe direction, since a 403 still triggers a
 * refresh. */
int64_t ParseIso8601(const std::string& text) {
  if (text.size() < 19) {
    return 0;
  }
  struct tm tm;
  memset(&tm, 0, sizeof(tm));
  int year, month, day, hour, minute, second;
  if (sscanf(text.c_str(), "%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour,
             &minute, &second) != 6) {
    return 0;
  }
  tm.tm_year = year - 1900;
  tm.tm_mon = month - 1;
  tm.tm_mday = day;
  tm.tm_hour = hour;
  tm.tm_min = minute;
  tm.tm_sec = second;
#ifdef _WIN32
  return static_cast<int64_t>(_mkgmtime(&tm));
#else
  return static_cast<int64_t>(timegm(&tm));
#endif
}

void FillFromJson(const std::string& json, Credentials* out) {
  out->access_key = JsonString(json, "AccessKeyId");
  out->secret_key = JsonString(json, "SecretAccessKey");
  out->session_token = JsonString(json, "Token");
  out->expires_at = ParseIso8601(JsonString(json, "Expiration"));
}

void TrimInPlace(std::string* s) {
  size_t begin = 0;
  while ((begin < s->size()) &&
         ((*s)[begin] == ' ' || (*s)[begin] == '\t' || (*s)[begin] == '\r')) {
    ++begin;
  }
  size_t end = s->size();
  while ((end > begin) && ((*s)[end - 1] == ' ' || (*s)[end - 1] == '\t' ||
                           (*s)[end - 1] == '\r' || (*s)[end - 1] == '\n')) {
    --end;
  }
  *s = s->substr(begin, end - begin);
}

std::string ActiveProfile() {
  std::string profile = GetEnv("AWS_PROFILE");
  if (profile.empty()) {
    profile = GetEnv("AWS_DEFAULT_PROFILE");
  }
  return profile.empty() ? std::string("default") : profile;
}

std::string HomeDirectory() {
  std::string home = GetEnv("HOME");
#ifdef _WIN32
  if (home.empty()) {
    home = GetEnv("USERPROFILE");
  }
#endif
  return home;
}

/* Reads one INI-style AWS config file and returns the key/value pairs of the
 * requested section.
 *
 * ~/.aws/config names non-default profiles "[profile foo]" while
 * ~/.aws/credentials names them "[foo]", so both spellings are accepted. */
std::map<std::string, std::string> ReadIniSection(const std::string& path,
                                                  const std::string& profile) {
  std::map<std::string, std::string> values;
  FILE* file = fopen(path.c_str(), "r");
  if (!file) {
    return values;
  }
  const std::string plain = "[" + profile + "]";
  const std::string prefixed = "[profile " + profile + "]";
  bool in_section = false;
  char line[2048];
  while (fgets(line, static_cast<int>(sizeof(line)), file)) {
    std::string text(line);
    TrimInPlace(&text);
    if (text.empty() || text[0] == '#' || text[0] == ';') {
      continue;
    }
    if (text[0] == '[') {
      in_section = (text == plain) || (text == prefixed);
      continue;
    }
    if (!in_section) {
      continue;
    }
    const size_t eq = text.find('=');
    if (eq == std::string::npos) {
      continue;
    }
    std::string key = text.substr(0, eq);
    std::string value = text.substr(eq + 1);
    TrimInPlace(&key);
    TrimInPlace(&value);
    if (!key.empty()) {
      values[key] = value;
    }
  }
  fclose(file);
  return values;
}

std::string SharedCredentialsPath() {
  const std::string override_path = GetEnv("AWS_SHARED_CREDENTIALS_FILE");
  if (!override_path.empty()) {
    return override_path;
  }
  const std::string home = HomeDirectory();
  return home.empty() ? std::string() : home + "/.aws/credentials";
}

std::string ConfigPath() {
  const std::string override_path = GetEnv("AWS_CONFIG_FILE");
  if (!override_path.empty()) {
    return override_path;
  }
  const std::string home = HomeDirectory();
  return home.empty() ? std::string() : home + "/.aws/config";
}

bool CredentialsFromEnv(Credentials* out) {
  out->access_key = GetEnv("AWS_ACCESS_KEY_ID");
  out->secret_key = GetEnv("AWS_SECRET_ACCESS_KEY");
  out->session_token = GetEnv("AWS_SESSION_TOKEN");
  out->expires_at = 0;
  return !out->Empty();
}

bool CredentialsFromFile(const std::string& path, Credentials* out) {
  if (path.empty()) {
    return false;
  }
  const std::map<std::string, std::string> values =
      ReadIniSection(path, ActiveProfile());
  std::map<std::string, std::string>::const_iterator it;
  if ((it = values.find("aws_access_key_id")) != values.end()) {
    out->access_key = it->second;
  }
  if ((it = values.find("aws_secret_access_key")) != values.end()) {
    out->secret_key = it->second;
  }
  if ((it = values.find("aws_session_token")) != values.end()) {
    out->session_token = it->second;
  }
  out->expires_at = 0;
  return !out->Empty();
}

/* ECS tasks and EKS pods with the credential-provider sidecar expose
 * credentials over a plain HTTP endpoint named by the environment. */
bool CredentialsFromContainer(Credentials* out) {
  const std::string relative =
      GetEnv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI");
  const std::string full = GetEnv("AWS_CONTAINER_CREDENTIALS_FULL_URI");
  std::string url;
  if (!relative.empty()) {
    url = "http://169.254.170.2" + relative;
  } else if (!full.empty()) {
    url = full;
  } else {
    return false;
  }

  std::string auth_header;
  const std::string token_file =
      GetEnv("AWS_CONTAINER_AUTHORIZATION_TOKEN_FILE");
  if (!token_file.empty()) {
    FILE* file = fopen(token_file.c_str(), "r");
    if (file) {
      char buffer[4096];
      const size_t got = fread(buffer, 1, sizeof(buffer) - 1, file);
      fclose(file);
      buffer[got] = '\0';
      std::string token(buffer);
      TrimInPlace(&token);
      if (!token.empty()) {
        auth_header = "Authorization: " + token;
      }
    }
  } else {
    const std::string token = GetEnv("AWS_CONTAINER_AUTHORIZATION_TOKEN");
    if (!token.empty()) {
      auth_header = "Authorization: " + token;
    }
  }

  std::string body;
  if (!FetchMetadata(url, nullptr, auth_header, &body)) {
    return false;
  }
  FillFromJson(body, out);
  return !out->Empty();
}

/* EC2 instance metadata, IMDSv2 only: obtain a session token, discover the
 * attached role, then read that role's credentials. */
bool CredentialsFromImds(Credentials* out) {
  if (GetEnvBool("AWS_EC2_METADATA_DISABLED")) {
    return false;
  }
  std::string base = GetEnv("AWS_EC2_METADATA_SERVICE_ENDPOINT");
  if (base.empty()) {
    base = "http://169.254.169.254";
  }
  while (!base.empty() && base[base.size() - 1] == '/') {
    base.erase(base.size() - 1);
  }

  std::string token;
  if (!FetchMetadata(base + "/latest/api/token", "PUT",
                     "X-aws-ec2-metadata-token-ttl-seconds: 21600", &token)) {
    return false;
  }
  const std::string token_header = "X-aws-ec2-metadata-token: " + token;

  std::string role;
  if (!FetchMetadata(base + "/latest/meta-data/iam/security-credentials/",
                     nullptr, token_header, &role)) {
    return false;
  }
  TrimInPlace(&role);
  /* The listing may hold several roles; the first is the one in effect. */
  const size_t newline = role.find('\n');
  if (newline != std::string::npos) {
    role = role.substr(0, newline);
    TrimInPlace(&role);
  }
  if (role.empty()) {
    return false;
  }

  std::string body;
  if (!FetchMetadata(base + "/latest/meta-data/iam/security-credentials/" + role,
                     nullptr, token_header, &body)) {
    return false;
  }
  FillFromJson(body, out);
  return !out->Empty();
}

}  // namespace

std::string GetEnv(const char* name) {
  const char* value = getenv(name);
  return value ? std::string(value) : std::string();
}

bool GetEnvBool(const char* name) {
  std::string value = GetEnv(name);
  for (size_t i = 0; i < value.size(); ++i) {
    if (value[i] >= 'A' && value[i] <= 'Z') {
      value[i] = static_cast<char>(value[i] - 'A' + 'a');
    }
  }
  return (value == "1") || (value == "true") || (value == "yes") ||
         (value == "on");
}

bool Credentials::Expired() const {
  if (expires_at == 0) {
    return false;
  }
  return static_cast<int64_t>(time(nullptr)) + kExpirySkewSeconds >= expires_at;
}

std::string ConfigValue(const char* key) {
  const std::map<std::string, std::string> values =
      ReadIniSection(ConfigPath(), ActiveProfile());
  const std::map<std::string, std::string>::const_iterator it =
      values.find(key);
  return (it == values.end()) ? std::string() : it->second;
}

bool ResolveCredentials(Credentials* out) {
  *out = Credentials();
  if (CredentialsFromEnv(out)) {
    return true;
  }
  *out = Credentials();
  if (CredentialsFromFile(SharedCredentialsPath(), out)) {
    return true;
  }
  *out = Credentials();
  if (CredentialsFromFile(ConfigPath(), out)) {
    return true;
  }
  *out = Credentials();
  if (CredentialsFromContainer(out)) {
    return true;
  }
  *out = Credentials();
  if (CredentialsFromImds(out)) {
    return true;
  }
  *out = Credentials();
  return false;
}

}  // namespace s3stream

#endif  // S3STREAM_ENABLE
