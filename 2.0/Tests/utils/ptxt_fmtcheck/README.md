# ptxt_fmtcheck

Checks `FormatNum()` in `2.0/utils/ptxt/ptxt.cc`, the correctly rounded
`%.<prec>g` formatter behind `ptxt print`. It has a fast path that scales a
value into an integer and rounds it, falling back to Ryu
(`2.0/ryu/`) within 8 ulps of a tie. This harness compares it against the
Ryu-only path and against printf at precisions 1–17, on random bit patterns,
short decimals, exact binary fractions, and values on or within 2 ulps of a
decimal tie. It also times the fast path, Ryu and printf.

```
make check            # about 1 minute on an M2 Pro; exits 1 on any mismatch
make check N=300000   # quicker
```

The Makefile extracts the formatter from `ptxt.cc` at build time, so the
harness always tests the current code. Run it after changing `FormatNum()`,
`FormatNumRyu()` or `LayoutG()`. The printf comparison assumes a correctly
rounding libc for non-integral values (glibc and macOS both qualify;
integral values are skipped, since ptxt prints those as integers).

Last result (M2 Pro, N=3000000): 46,325,458 non-integral values, no
mismatches; `%.6g` fast path 13.8 ns, Ryu 40.0 ns, snprintf 106.8 ns.
