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
#  include <errno.h>
#  include <unistd.h>  // close, mkstemp, unlink
#  include <fstream>
#  include "plink2_cmdline.h"
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

void S3Init() {
  g_s3_sdk_options.loggingOptions.logLevel = Aws::Utils::Logging::LogLevel::Off;
  Aws::InitAPI(g_s3_sdk_options);
}

void S3Shutdown() {
  Aws::ShutdownAPI(g_s3_sdk_options);
}

// Parse "s3://bucket/key" into bucket and key components.
static void ParseS3Uri(const char* s3_uri, Aws::String* bucket,
                       Aws::String* key) {
  // Skip the leading "s3://"
  const char* rest = s3_uri + 5;
  const char* slash = strchr(rest, '/');
  if (slash) {
    *bucket = Aws::String(rest, static_cast<size_t>(slash - rest));
    *key = Aws::String(slash + 1);
  } else {
    *bucket = Aws::String(rest);
    *key = Aws::String();
  }
}

PglErr S3DownloadToTemp(const char* s3_uri, char* local_path_buf) {
  Aws::String bucket;
  Aws::String key;
  ParseS3Uri(s3_uri, &bucket, &key);

  // Create a temporary file to hold the downloaded object.
  // mkstemp requires a writable template with at least 6 trailing 'X' chars.
  snprintf(local_path_buf, kPglFnamesize, "/tmp/plink2_s3_XXXXXX");
  const int fd = mkstemp(local_path_buf);
  if (fd < 0) {
    logerrprintfww("Error: Failed to create temporary file for S3 download: %s\n",
                   strerror(errno));
    return kPglRetOpenFail;
  }
  close(fd);

  // Set up the S3 client using the default credential provider chain
  // (environment variables, ~/.aws/credentials, IAM role, etc.).
  Aws::Client::ClientConfiguration client_config;
  Aws::S3::S3Client s3_client(client_config);

  Aws::S3::Model::GetObjectRequest request;
  request.SetBucket(bucket);
  request.SetKey(key);

  logprintfww("Downloading %s ...\n", s3_uri);
  auto outcome = s3_client.GetObject(request);
  if (!outcome.IsSuccess()) {
    const auto& error = outcome.GetError();
    logerrprintfww("Error: Failed to download %s: %s\n", s3_uri,
                   error.GetMessage().c_str());
    unlink(local_path_buf);
    return kPglRetOpenFail;
  }

  // Write the response body to the temporary file.
  std::ofstream outfile(local_path_buf, std::ios::binary | std::ios::trunc);
  if (!outfile.is_open()) {
    logerrprintfww("Error: Failed to open temporary file %s for writing: %s\n",
                   local_path_buf, strerror(errno));
    unlink(local_path_buf);
    return kPglRetOpenFail;
  }
  outfile << outcome.GetResult().GetBody().rdbuf();
  if (outfile.fail()) {
    logerrprintfww("Error: Failed to write downloaded S3 object to %s\n",
                   local_path_buf);
    outfile.close();
    unlink(local_path_buf);
    return kPglRetWriteFail;
  }
  outfile.close();

  return kPglRetSuccess;
}

#endif  // USE_S3

#ifdef __cplusplus
}  // namespace plink2
#endif
