/* Credential resolution for s3stream: environment, shared config files,
 * credential_process, ECS container endpoint, EC2 instance metadata.
 *
 * Deliberately unsupported, because each needs a substantial amount of
 * machinery for a case that has a one-line workaround (run the AWS CLI and
 * export the result, or point credential_process at
 * `aws configure export-credentials`):
 *   - SSO (`aws sso login`)
 *   - role_arn / source_profile AssumeRole chaining
 *   - web identity / IRSA
 */

#include "s3stream_internal.h"

#ifdef S3STREAM_ENABLE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <curl/curl.h>

#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#else
#  include <errno.h>
#  include <fcntl.h>
#  include <poll.h>
#  include <signal.h>
#  include <spawn.h>
#  include <sys/wait.h>
#  include <unistd.h>
#  ifdef __APPLE__
#    include <crt_externs.h>
#  else
extern char** environ;
#  endif
#endif

namespace s3stream {

namespace {

/* Refresh this long before the advertised expiry, so credentials do not
 * lapse in the middle of a long read. */
const int64_t kExpirySkewSeconds = 300;

/* The link-local metadata endpoints are unreachable off EC2/ECS, so keep the
 * probe short enough that it does not stall an otherwise fine failure path. */
const long kMetadataTimeoutSeconds = 2;

/* Credential documents are a few hundred bytes.  Whatever answers
 * 169.254.169.254 -- or AWS_CONTAINER_CREDENTIALS_FULL_URI, which the
 * environment can point anywhere -- does not get to decide how much memory
 * this process spends. */
const size_t kMaxMetadataBytes = 64 * 1024;

size_t CollectToString(char* ptr, size_t size, size_t nmemb, void* userdata) {
  const size_t n = size * nmemb;
  std::string* body = static_cast<std::string*>(userdata);
  if (body->size() + n > kMaxMetadataBytes) {
    return 0;  /* short return aborts the transfer */
  }
  body->append(ptr, n);
  return n;
}

/* Fetches a metadata URL.  Returns false on any transport or HTTP error. */
bool FetchMetadata(const std::string& url, const char* method,
                   const std::string& header, std::string* body) {
  CURL* curl = curl_easy_init();
  if (!curl) {
    return false;
  }
  /* AWS_CONTAINER_CREDENTIALS_FULL_URI is an arbitrary URL from the
   * environment; it must not be able to name file:// or anything else. */
  if (!RestrictToHttp(curl)) {
    curl_easy_cleanup(curl);
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

/* A credential_process helper may legitimately prompt or contact an identity
 * provider, so this is far longer than the metadata-endpoint timeout. */
const int kProcessTimeoutSeconds = 60;

#ifndef _WIN32
/* Splits like a POSIX shell word list (whitespace, '...', "...", backslash)
 * but performs no expansion; the command is exec'd directly, never through a
 * shell. */
bool SplitCommand(const std::string& command,
                  std::vector<std::string>* args) {
  std::string current;
  bool have = false;
  char quote = 0;
  for (size_t i = 0; i < command.size(); ++i) {
    const char c = command[i];
    if (quote == '\'') {
      if (c == '\'') {
        quote = 0;
      } else {
        current += c;
      }
    } else if (quote == '"') {
      if (c == '"') {
        quote = 0;
      } else if ((c == '\\') && (i + 1 < command.size()) &&
                 ((command[i + 1] == '"') || (command[i + 1] == '\\'))) {
        current += command[++i];
      } else {
        current += c;
      }
    } else if ((c == '\'') || (c == '"')) {
      quote = c;
      have = true;
    } else if ((c == '\\') && (i + 1 < command.size())) {
      current += command[++i];
      have = true;
    } else if ((c == ' ') || (c == '\t')) {
      if (have) {
        args->push_back(current);
        current.clear();
        have = false;
      }
    } else {
      current += c;
      have = true;
    }
  }
  if (quote) {
    return false;
  }
  if (have) {
    args->push_back(current);
  }
  return !args->empty();
}

/* Runs `command` with stdin at /dev/null and stderr inherited, capturing
 * stdout.  True only for a clean exit status 0 within the timeout and size
 * cap. */
bool RunProcess(const std::string& command, std::string* out) {
  std::vector<std::string> args;
  if (!SplitCommand(command, &args)) {
    return false;
  }
  std::vector<char*> argv;
  for (size_t i = 0; i < args.size(); ++i) {
    argv.push_back(const_cast<char*>(args[i].c_str()));
  }
  argv.push_back(nullptr);

  int fds[2];
  if (pipe(fds) != 0) {
    return false;
  }
  fcntl(fds[0], F_SETFD, FD_CLOEXEC);
  fcntl(fds[1], F_SETFD, FD_CLOEXEC);

  posix_spawn_file_actions_t actions;
  if (posix_spawn_file_actions_init(&actions) != 0) {
    close(fds[0]);
    close(fds[1]);
    return false;
  }
  posix_spawn_file_actions_adddup2(&actions, fds[1], STDOUT_FILENO);
  posix_spawn_file_actions_addopen(&actions, STDIN_FILENO, "/dev/null",
                                   O_RDONLY, 0);
#ifdef __APPLE__
  char** const env = *_NSGetEnviron();
#else
  char** const env = environ;
#endif
  pid_t pid = 0;
  const int rc = posix_spawnp(&pid, argv[0], &actions, nullptr, argv.data(),
                              env);
  posix_spawn_file_actions_destroy(&actions);
  close(fds[1]);
  if (rc != 0) {
    close(fds[0]);
    return false;
  }

  const time_t deadline = time(nullptr) + kProcessTimeoutSeconds;
  bool ok = true;
  for (;;) {
    const time_t remaining = deadline - time(nullptr);
    if (remaining <= 0) {
      ok = false;
      break;
    }
    struct pollfd pfd;
    pfd.fd = fds[0];
    pfd.events = POLLIN;
    pfd.revents = 0;
    const int ready = poll(&pfd, 1, static_cast<int>(remaining * 1000));
    if (ready < 0) {
      if (errno == EINTR) {
        continue;
      }
      ok = false;
      break;
    }
    if (ready == 0) {
      continue;
    }
    char buffer[4096];
    const ssize_t got = read(fds[0], buffer, sizeof(buffer));
    if (got < 0) {
      if (errno == EINTR) {
        continue;
      }
      ok = false;
      break;
    }
    if (got == 0) {
      break;
    }
    if (out->size() + static_cast<size_t>(got) > kMaxMetadataBytes) {
      ok = false;
      break;
    }
    out->append(buffer, static_cast<size_t>(got));
  }
  close(fds[0]);
  if (!ok) {
    kill(pid, SIGKILL);
  }
  int wstatus = 0;
  while ((waitpid(pid, &wstatus, 0) < 0) && (errno == EINTR)) {
  }
  return ok && WIFEXITED(wstatus) && (WEXITSTATUS(wstatus) == 0);
}
#else  /* _WIN32 */
/* The command line is handed to CreateProcess as written, as the AWS SDKs do
 * on Windows. */
bool RunProcess(const std::string& command, std::string* out) {
  SECURITY_ATTRIBUTES sa;
  memset(&sa, 0, sizeof(sa));
  sa.nLength = sizeof(sa);
  sa.bInheritHandle = TRUE;
  HANDLE read_end = nullptr;
  HANDLE write_end = nullptr;
  if (!CreatePipe(&read_end, &write_end, &sa, 0)) {
    return false;
  }
  SetHandleInformation(read_end, HANDLE_FLAG_INHERIT, 0);
  HANDLE nul = CreateFileA("NUL", GENERIC_READ,
                           FILE_SHARE_READ | FILE_SHARE_WRITE, &sa,
                           OPEN_EXISTING, 0, nullptr);

  STARTUPINFOA si;
  memset(&si, 0, sizeof(si));
  si.cb = sizeof(si);
  si.dwFlags = STARTF_USESTDHANDLES;
  si.hStdInput = nul;
  si.hStdOutput = write_end;
  si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
  PROCESS_INFORMATION pi;
  memset(&pi, 0, sizeof(pi));
  std::vector<char> line(command.begin(), command.end());
  line.push_back('\0');
  const BOOL created =
      CreateProcessA(nullptr, line.data(), nullptr, nullptr, TRUE,
                     CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi);
  CloseHandle(write_end);
  if (nul != INVALID_HANDLE_VALUE) {
    CloseHandle(nul);
  }
  if (!created) {
    CloseHandle(read_end);
    return false;
  }

  const ULONGLONG deadline =
      GetTickCount64() + static_cast<ULONGLONG>(kProcessTimeoutSeconds) * 1000;
  bool ok = true;
  bool exited = false;
  for (;;) {
    DWORD available = 0;
    if (!PeekNamedPipe(read_end, nullptr, 0, nullptr, &available, nullptr)) {
      break;  /* all write ends closed */
    }
    if (available == 0) {
      if (exited) {
        break;
      }
      if (WaitForSingleObject(pi.hProcess, 0) == WAIT_OBJECT_0) {
        exited = true;  /* drain anything written just before exit */
        continue;
      }
      if (GetTickCount64() >= deadline) {
        ok = false;
        break;
      }
      Sleep(20);
      continue;
    }
    char buffer[4096];
    DWORD want = available < sizeof(buffer) ? available : sizeof(buffer);
    DWORD got = 0;
    if (!ReadFile(read_end, buffer, want, &got, nullptr) || (got == 0)) {
      break;
    }
    if (out->size() + got > kMaxMetadataBytes) {
      ok = false;
      break;
    }
    out->append(buffer, got);
  }
  CloseHandle(read_end);
  if (!ok) {
    TerminateProcess(pi.hProcess, 1);
  }
  DWORD code = 1;
  if (WaitForSingleObject(pi.hProcess, ok ? 5000 : INFINITE) ==
      WAIT_OBJECT_0) {
    GetExitCodeProcess(pi.hProcess, &code);
  } else {
    TerminateProcess(pi.hProcess, 1);
  }
  CloseHandle(pi.hThread);
  CloseHandle(pi.hProcess);
  return ok && (code == 0);
}
#endif

/* Reads an integer field, or -1 when absent/unparseable. */
long JsonInt(const std::string& json, const char* key) {
  const std::string needle = std::string("\"") + key + "\"";
  size_t pos = json.find(needle);
  if (pos == std::string::npos) {
    return -1;
  }
  pos = json.find(':', pos + needle.size());
  if (pos == std::string::npos) {
    return -1;
  }
  ++pos;
  while ((pos < json.size()) && ((json[pos] == ' ') || (json[pos] == '\t') ||
                                 (json[pos] == '\r') || (json[pos] == '\n'))) {
    ++pos;
  }
  char* end = nullptr;
  const long value = strtol(json.c_str() + pos, &end, 10);
  return (end == json.c_str() + pos) ? -1 : value;
}

/* Runs a credential_process command and parses its Version 1 JSON document
 * (AccessKeyId, SecretAccessKey, optional SessionToken and Expiration). */
bool CredentialsFromCommand(const std::string& command, Credentials* out) {
  std::string body;
  const bool ran = RunProcess(command, &body);
  const bool valid = ran && (JsonInt(body, "Version") == 1);
  if (valid) {
    out->access_key = JsonString(body, "AccessKeyId");
    out->secret_key = JsonString(body, "SecretAccessKey");
    out->session_token = JsonString(body, "SessionToken");
    out->expires_at = ParseIso8601(JsonString(body, "Expiration"));
  }
  SecureZero(&body);
  return valid && !out->Empty();
}

/* credential_process may sit in either the credentials file or the config
 * file of the active profile. */
bool CredentialsFromProcess(Credentials* out) {
  const std::string paths[2] = {SharedCredentialsPath(), ConfigPath()};
  for (int i = 0; i < 2; ++i) {
    if (paths[i].empty()) {
      continue;
    }
    const std::map<std::string, std::string> values =
        ReadIniSection(paths[i], ActiveProfile());
    const std::map<std::string, std::string>::const_iterator it =
        values.find("credential_process");
    if ((it == values.end()) || it->second.empty()) {
      continue;
    }
    if (CredentialsFromCommand(it->second, out)) {
      return true;
    }
    *out = Credentials();
  }
  return false;
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
  TrimInPlace(&token);
  /* The token is sent back as a header on the next two requests; a CR or LF
   * in it would splice extra requests into our own connection. */
  if (!IsSafeHeaderValue(token, 4096)) {
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
  /* The role name is appended to a URL path, so restrict it to what an IAM
   * role name can actually contain rather than letting the response steer the
   * next request somewhere else. */
  for (size_t i = 0; i < role.size(); ++i) {
    const char c = role[i];
    const bool ok = ((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z')) ||
                    ((c >= '0') && (c <= '9')) || (c == '+') || (c == '=') ||
                    (c == ',') || (c == '.') || (c == '@') || (c == '-') ||
                    (c == '_');
    if (!ok) {
      return false;
    }
  }
  if (role.size() > 64) {
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
  if (CredentialsFromProcess(out)) {
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
