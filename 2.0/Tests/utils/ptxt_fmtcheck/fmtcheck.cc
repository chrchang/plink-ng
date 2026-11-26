// Checks FormatNum() in 2.0/utils/ptxt/ptxt.cc, the %.<prec>g formatter
// behind 'ptxt print', against its Ryu-only path (FormatNumRyu) and against
// printf, at every precision from 1 to 17.  The functions are pulled out of
// ptxt.cc at build time (see Makefile), so this always tests the real code.
//
// Usage: fmtcheck [N]   N values per pattern family per precision
//                       (default 3000000; about 1 minute on an M2 Pro).
// Exits 1 on any mismatch.
//
// Comparing against printf assumes a correctly rounding libc, which holds
// for glibc and macOS on non-integral values.  (Integral values are skipped:
// FormatNum prints them as integers, like awk.)

#include "plink2_base.h"
#include "ryu/ryu.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <chrono>
using namespace plink2;
#include "fmt_extracted.h"
static uint64_t s = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
int main(int argc, char** argv) {
  const long per = (argc > 1) ? atol(argv[1]) : 3000000;
  long bad_ryu = 0, bad_pf = 0, n = 0, nonint = 0;
  char a[128], b[128], c[128];
  for (uint32_t prec = 1; prec <= 17; ++prec) {
    for (long i = 0; i < per; ++i) {
      double d;
      switch (i % 6) {
        case 0: { uint64_t bits = rnd(); memcpy(&d, &bits, 8); break; }
        case 1: d = (double)(rnd() % 100000000) / pow(10, rnd() % 12); break;
        case 2: d = ((double)(rnd() % 2000000) + 0.5) / pow(10, rnd() % 10); break;
        case 3: d = (rnd() % 100000000) * ldexp(1.0, (int)(rnd() % 80) - 60); break;
        case 4: {  // decimal ties at the prec-th digit, nudged by 0..2 ulps
          double base = (double)(rnd() % 1000000000000000ULL) + 0.5;
          d = base * pow(10, (int)(rnd() % 40) - 25);
          int nudge = (int)(rnd() % 5) - 2;
          for (; nudge > 0; --nudge) d = nextafter(d, INFINITY);
          for (; nudge < 0; ++nudge) d = nextafter(d, 0);
          break;
        }
        default: d = (double)(rnd() % 1000) / 7.0 * pow(10, (int)(rnd() % 30) - 15); break;
      }
      if (rnd() & 1) d = -d;
      if (!(d == d) || isinf(d)) continue;
      ++n;
      const uint32_t la = FormatNum(d, prec, a); a[la] = 0;
      if (d == floor(d) && fabs(d) < 9.2e18) continue;  // integer path, not under test
      ++nonint;
      const uint32_t lb = FormatNumRyu(d, prec, b); b[lb] = 0;
      snprintf(c, sizeof c, "%.*g", (int)prec, d);
      if (strcmp(a, b)) { if (bad_ryu < 8) printf("vs Ryu: prec %u %.17g fast=%s ryu=%s\n", prec, d, a, b); ++bad_ryu; }
      if (strcmp(a, c)) { if (bad_pf < 8) printf("vs printf: prec %u %.17g fast=%s printf=%s\n", prec, d, a, c); ++bad_pf; }
    }
  }
  printf("%ld values, %ld non-integral: %ld mismatches vs Ryu, %ld vs printf\n", n, nonint, bad_ryu, bad_pf);
  // Timing: %.6g on typical values.
  char buf[64]; double acc = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int i = 1; i < 10000000; ++i) acc += FormatNum(i * 1.37e-3 + 1e-9, 6, buf);
  auto t1 = std::chrono::steady_clock::now();
  for (int i = 1; i < 10000000; ++i) acc += FormatNumRyu(i * 1.37e-3 + 1e-9, 6, buf);
  auto t2 = std::chrono::steady_clock::now();
  for (int i = 1; i < 10000000; ++i) acc += snprintf(buf, 64, "%.6g", i * 1.37e-3 + 1e-9);
  auto t3 = std::chrono::steady_clock::now();
  auto ns = [](auto x, auto y) { return std::chrono::duration<double, std::nano>(y - x).count() / 1e7; };
  printf("%%.6g: fast %.1f ns, Ryu %.1f ns, snprintf %.1f ns (%g)\n", ns(t0, t1), ns(t1, t2), ns(t2, t3), acc);
  return (bad_ryu || bad_pf) ? 1 : 0;
}
