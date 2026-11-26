// Unit test for the 64-bit deltalist encoding used by multiallelic dosage
// tracks #5 and #9 (pgenlib_misc.h): round trips, exact byte layout, and
// rejection of malformed input.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "pgenlib_misc.h"

using namespace plink2;

static int g_fail_ct = 0;

#define CHECK(cond) do { if (!(cond)) { fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #cond); ++g_fail_ct; } } while (0)

static uint64_t g_rng = 0x9e3779b97f4a7c15LLU;
static uint64_t NextRand() {
  // splitmix64
  uint64_t z = (g_rng += 0x9e3779b97f4a7c15LLU);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9LLU;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebLLU;
  return z ^ (z >> 31);
}

// Encodes idxs into a buffer of exactly the encoded size (so that ASan can
// catch any read past the end), returning it.
static std::vector<unsigned char> Encode(const std::vector<uint64_t>& idxs, uint64_t idx_end) {
  const uint32_t idx_byte_ct = Deltalist64IdxByteCt(idx_end);
  std::vector<unsigned char> buf(Deltalist64MaxByteCt(idxs.size(), idx_byte_ct));
  unsigned char* end = Deltalist64Append(idxs.data(), idxs.size(), idx_byte_ct, buf.data());
  CHECK(end <= buf.data() + buf.size());
  buf.resize(end - buf.data());
  return buf;
}

// Parses from a heap copy sized exactly to the input.
static BoolErr Parse(const std::vector<unsigned char>& bytes, uint64_t idx_end, uint32_t max_len, std::vector<uint64_t>* out, uint32_t* consumed_ptr) {
  unsigned char* copy = static_cast<unsigned char*>(malloc(bytes.size() + 1));
  memcpy(copy, bytes.data(), bytes.size());
  const unsigned char* iter = copy;
  uint32_t len = 0;
  std::vector<uint64_t> idxs(max_len + 1);
  BoolErr err = ParseDeltalist64(copy + bytes.size(), idx_end, max_len, &iter, idxs.data(), &len);
  if (!err) {
    idxs.resize(len);
    if (out) {
      *out = idxs;
    }
    *consumed_ptr = iter - copy;
    // Count-only mode must agree.
    const unsigned char* iter2 = copy;
    uint32_t len2 = 0;
    CHECK(!ParseDeltalist64(copy + bytes.size(), idx_end, max_len, &iter2, nullptr, &len2));
    CHECK(len2 == len);
    CHECK(iter2 == iter);
  }
  free(copy);
  return err;
}

static void RoundTrip(const std::vector<uint64_t>& idxs, uint64_t idx_end) {
  std::vector<unsigned char> bytes = Encode(idxs, idx_end);
  std::vector<uint64_t> decoded;
  uint32_t consumed = 0;
  CHECK(!Parse(bytes, idx_end, idxs.size(), &decoded, &consumed));
  CHECK(decoded == idxs);
  CHECK(consumed == bytes.size());
  // max_len one too small must fail.
  if (!idxs.empty()) {
    CHECK(Parse(bytes, idx_end, idxs.size() - 1, nullptr, &consumed));
  }
  // Every proper prefix must be rejected.
  for (size_t prefix = 0; prefix < bytes.size(); ++prefix) {
    std::vector<unsigned char> trunc(bytes.begin(), bytes.begin() + prefix);
    CHECK(Parse(trunc, idx_end, idxs.size(), nullptr, &consumed));
  }
}

static std::vector<uint64_t> RandomList(uint32_t len, uint64_t idx_end, uint32_t big_deltas) {
  std::vector<uint64_t> idxs;
  uint64_t cur = NextRand() % 64;
  for (uint32_t i = 0; i != len; ++i) {
    if (cur >= idx_end) {
      break;
    }
    idxs.push_back(cur);
    uint64_t room = idx_end - cur;
    uint64_t step;
    if (big_deltas) {
      // Spread the list over the whole index space, so deltas can need up to
      // 9 varint bytes.
      step = 1 + (room / (len - i + 1));
    } else {
      step = 1 + (NextRand() % 300);
    }
    if (step >= room) {
      break;
    }
    cur += step;
  }
  return idxs;
}

int main() {
  // RareAltIdxBitCt(allele_ct) = ceil(log2(allele_ct - 2)).
  for (uint32_t allele_ct = 3; allele_ct <= 255; ++allele_ct) {
    uint32_t b = 0;
    while ((1U << b) < allele_ct - 2) {
      ++b;
    }
    CHECK(RareAltIdxBitCt(allele_ct) == b);
  }
  CHECK(Deltalist64IdxByteCt(255) == 1);
  CHECK(Deltalist64IdxByteCt(256) == 2);
  CHECK(Deltalist64IdxByteCt(1LLU << 39) == 5);
  CHECK(Deltalist64IdxByteCt((1LLU << 56) - 1) == 7);
  CHECK(Deltalist64IdxByteCt(1LLU << 56) == 8);

  // Exact layout, one group: idx_end = 1000 << 2 = 4000 (2-byte group starts),
  // indexes {5, 9, 300}: len 3, start 5, deltas 4 and 291.
  {
    const std::vector<uint64_t> idxs = {5, 9, 300};
    const std::vector<unsigned char> expected = {0x03, 0x05, 0x00, 0x04, 0xa3, 0x02};
    CHECK(Encode(idxs, 4000) == expected);
  }
  // Exact layout, two groups: indexes 0, 1, ..., 64 with idx_end = 100
  // (1-byte group starts).  Group 0 has 63 one-byte deltas, so its stored
  // uint16 byte count is 63.
  {
    std::vector<uint64_t> idxs;
    for (uint64_t i = 0; i != 65; ++i) {
      idxs.push_back(i);
    }
    std::vector<unsigned char> expected = {65, 0x00, 64, 63, 0x00};
    for (uint32_t i = 0; i != 63; ++i) {
      expected.push_back(0x01);
    }
    CHECK(Encode(idxs, 100) == expected);
  }

  // Round trips.
  const uint64_t idx_ends[] = {1, 8, 300, 65536, 1LLU << 20, 488377LLU << 3, 1LLU << 39, (1LLU << 55) + 12345};
  const uint32_t lens[] = {0, 1, 2, 63, 64, 65, 127, 128, 129, 1000, 4097};
  for (uint64_t idx_end : idx_ends) {
    for (uint32_t len : lens) {
      for (uint32_t big = 0; big != 2; ++big) {
        RoundTrip(RandomList(len, idx_end, big), idx_end);
      }
    }
  }
  // Last valid index; and a full group of 9-byte deltas, whose byte count
  // (567) doesn't fit in the old one-byte field.
  RoundTrip({(1LLU << 55) - 1}, 1LLU << 55);
  {
    std::vector<uint64_t> idxs;
    for (uint64_t i = 0; i != 65; ++i) {
      idxs.push_back(i << 56);
    }
    std::vector<unsigned char> bytes = Encode(idxs, 1LLU << 63);
    CHECK(bytes[1 + 2 * 8] == (567 & 255));
    CHECK(bytes[1 + 2 * 8 + 1] == (567 >> 8));
    RoundTrip(idxs, 1LLU << 63);
  }

  // Malformed input.
  uint32_t consumed;
  {
    // Index == idx_end.
    std::vector<unsigned char> bytes = Encode({5, 9, 300}, 4000);
    CHECK(Parse(bytes, 300, 3, nullptr, &consumed));
    CHECK(!Parse(bytes, 301, 3, nullptr, &consumed));
    // Group start == idx_end.
    CHECK(Parse(Encode({300}, 4000), 300, 1, nullptr, &consumed));
  }
  {
    // Zero delta: {5, 5, 300} hand-encoded.
    const std::vector<unsigned char> bytes = {0x03, 0x05, 0x00, 0x00, 0xa3, 0x02};
    CHECK(Parse(bytes, 4000, 3, nullptr, &consumed));
  }
  {
    // Second group start not greater than the last index of the first group.
    std::vector<uint64_t> idxs;
    for (uint64_t i = 0; i != 65; ++i) {
      idxs.push_back(i);
    }
    std::vector<unsigned char> bytes = Encode(idxs, 100);
    CHECK(!Parse(bytes, 100, 65, nullptr, &consumed));
    bytes[2] = 62;  // second group start: 64 -> 62
    CHECK(Parse(bytes, 100, 65, nullptr, &consumed));
    bytes[2] = 64;
    // Wrong stored group byte count.
    bytes[3] = 62;
    CHECK(Parse(bytes, 100, 65, nullptr, &consumed));
  }
  {
    // Length varint that doesn't fit in 63 bits.
    std::vector<unsigned char> bytes(9, 0xff);
    bytes.push_back(0x01);
    CHECK(Parse(bytes, 4000, 0x7fffffff, nullptr, &consumed));
    // Delta varint that doesn't fit in 63 bits.
    std::vector<unsigned char> bytes2 = {0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    for (uint32_t i = 0; i != 9; ++i) {
      bytes2.push_back(0xff);
    }
    bytes2.push_back(0x01);
    CHECK(Parse(bytes2, 1LLU << 63, 2, nullptr, &consumed));
  }

  if (g_fail_ct) {
    fprintf(stderr, "%d check(s) failed.\n", g_fail_ct);
    return 1;
  }
  printf("deltalist64 tests passed.\n");
  return 0;
}
