#include "../../include/s3stream/s3stream.h"

#include <arpa/inet.h>
#include <poll.h>
#include <signal.h>
#include <sys/socket.h>
#include <unistd.h>

#include <atomic>
#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

class HttpServer {
 public:
  std::atomic<int> connections{0};
  std::atomic<int> metadata_requests{0};
  std::atomic<int> generation{1};
  std::atomic<bool> metadata_failure{false};
  std::atomic<bool> metadata_expired{false};

  HttpServer() {
    listener_ = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in address{};
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    if ((listener_ < 0) ||
        (bind(listener_, reinterpret_cast<sockaddr*>(&address), sizeof(address)) < 0) ||
        (listen(listener_, 16) < 0)) {
      throw std::runtime_error("cannot start loopback HTTP server");
    }
    socklen_t address_size = sizeof(address);
    getsockname(listener_, reinterpret_cast<sockaddr*>(&address), &address_size);
    url = "http://127.0.0.1:" + std::to_string(ntohs(address.sin_port));
    worker_ = std::thread([this] { Run(); });
  }

  ~HttpServer() {
    stopping_ = true;
    worker_.join();
    close(listener_);
  }

  bool SawAuthorization(const std::string& value) {
    std::lock_guard<std::mutex> lock(requests_mutex_);
    for (const auto& request : requests_) {
      if (request.find(value) != std::string::npos) {
        return true;
      }
    }
    return false;
  }

  std::string url;

 private:
  int listener_;
  std::atomic<bool> stopping_{false};
  std::thread worker_;
  std::mutex requests_mutex_;
  std::vector<std::string> requests_;

  std::string Response(const std::string& request) {
    {
      std::lock_guard<std::mutex> lock(requests_mutex_);
      requests_.push_back(request);
    }
    const bool head = request.compare(0, 5, "HEAD ") == 0;
    int status_code = 200;
    std::string body;
    std::string headers;
    size_t length = 1048576;
    if (request.find(" /metadata") != std::string::npos) {
      ++metadata_requests;
      if (metadata_failure) {
        status_code = 404;
      } else {
        body = "{\"AccessKeyId\":\"CACHE" + std::to_string(generation.load()) +
               "\",\"SecretAccessKey\":\"secret\",\"Token\":\"token" +
               std::to_string(generation.load()) + "\",\"Expiration\":\"" +
               (metadata_expired ? "2000-01-01" : "2099-01-01") +
               "T00:00:00Z\"}";
      }
      length = body.size();
    } else if (request.find("/missing") != std::string::npos) {
      status_code = 404;
      length = 0;
    } else if (request.find("/denied") != std::string::npos) {
      status_code = 403;
      length = 0;
    } else if ((request.find("/anonymous") != std::string::npos) &&
               ((request.find("Authorization:") != std::string::npos) ||
                (request.find("x-amz-security-token:") != std::string::npos))) {
      status_code = 403;
      length = 0;
    } else if (!head && (request.find("/refresh") != std::string::npos) &&
               (request.find("Credential=CACHE1/") != std::string::npos)) {
      status_code = 403;
      body = "<Error><Code>ExpiredToken</Code></Error>";
      length = body.size();
    } else {
      headers = "ETag: \"fixture\"\r\n";
      const size_t range = request.find("\r\nRange: bytes=");
      if (!head && (range != std::string::npos)) {
        size_t first = 0;
        size_t last = 0;
        sscanf(request.c_str() + range, "\r\nRange: bytes=%zu-%zu", &first, &last);
        if ((first >= length) || (last < first) || (last >= length)) {
          status_code = 416;
          length = 0;
        } else {
          status_code = 206;
          headers += "Content-Range: bytes " + std::to_string(first) + "-" +
                     std::to_string(last) + "/" + std::to_string(length) + "\r\n";
          length = last - first + 1;
          body.assign(length, 'x');
        }
      }
    }
    return "HTTP/1.1 " + std::to_string(status_code) + " Test\r\n" + headers +
           "Content-Length: " + std::to_string(length) + "\r\n\r\n" +
           (head ? std::string() : body);
  }

  void Run() {
    std::vector<int> clients;
    std::vector<std::string> pending;
    while (!stopping_) {
      std::vector<pollfd> watches{{listener_, POLLIN, 0}};
      for (int client : clients) {
        watches.push_back({client, POLLIN, 0});
      }
      if (poll(watches.data(), watches.size(), 100) <= 0) {
        continue;
      }
      for (size_t client_index = clients.size(); client_index > 0;) {
        --client_index;
        if (!watches[client_index + 1].revents) {
          continue;
        }
        char buffer[4096];
        const ssize_t count = recv(clients[client_index], buffer, sizeof(buffer), 0);
        if (count <= 0) {
          close(clients[client_index]);
          clients.erase(clients.begin() + client_index);
          pending.erase(pending.begin() + client_index);
          continue;
        }
        pending[client_index].append(buffer, static_cast<size_t>(count));
        const size_t end = pending[client_index].find("\r\n\r\n");
        if (end == std::string::npos) {
          continue;
        }
        const std::string response = Response(pending[client_index].substr(0, end + 4));
        pending[client_index].erase(0, end + 4);
        size_t sent = 0;
        while (sent < response.size()) {
          const ssize_t count_sent = send(clients[client_index], response.data() + sent,
                                          response.size() - sent, 0);
          if (count_sent <= 0) {
            break;
          }
          sent += static_cast<size_t>(count_sent);
        }
      }
      if (watches[0].revents & POLLIN) {
        const int client = accept(listener_, nullptr, nullptr);
        if (client >= 0) {
          clients.push_back(client);
          pending.emplace_back();
          ++connections;
        }
      }
    }
    for (int client : clients) {
      close(client);
    }
  }
};

int main() {
  signal(SIGPIPE, SIG_IGN);
  setenv("NO_PROXY", "127.0.0.1", 1);
  setenv("no_proxy", "127.0.0.1", 1);
  assert(s3stream_init() == 0);
  errno = 0;
  assert(s3stream_open(nullptr) == nullptr);
  assert(errno == EINVAL);
  assert(s3stream_open("s3://invalid") == nullptr);
  assert(errno == EINVAL);
  assert(s3stream_open("/nonexistent/plink-fixture") == nullptr);
  assert(errno == ENOENT);

  HttpServer server;
  s3stream_credentials failure_creds{};
  failure_creds.endpoint_url = server.url.c_str();
  failure_creds.region = "us-east-1";
  failure_creds.no_sign_request = 1;
  failure_creds.force_path_style = 1;
  for (const auto& failure : std::vector<std::pair<std::string, int>>{
           {"missing", ENOENT}, {"denied", EACCES}}) {
    const std::string url = server.url + "/" + failure.first;
    const std::string uri = "s3://bucket/" + failure.first;
    errno = 0;
    assert(s3stream_open(url.c_str()) == nullptr);
    assert(errno == failure.second);
    assert(std::string(s3stream_last_error()).find("HTTP ") != std::string::npos);
    errno = 0;
    assert(s3stream_open_with_credentials(uri.c_str(), &failure_creds) == nullptr);
    assert(errno == failure.second);
    assert(std::string(s3stream_last_error()).find("HTTP ") != std::string::npos);
  }

  puts("s3stream errno checks passed.");
  return 0;
}
