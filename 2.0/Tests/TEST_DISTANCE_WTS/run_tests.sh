#!/bin/bash

set -exo pipefail

# 60 samples, 300 variants.  Two datasets: one complete, for the differential
# against plink 1.9, and one with 5% missing calls, for the oracle and the
# equal-weights identity.
awk 'BEGIN {
  for (i = 0; i < 60; ++i) {
    printf "fam%d per%d 0 0 %d 1\n", i, i, 1 + (i % 2) > "tmp_data.fam"
  }
  for (j = 0; j < 300; ++j) {
    printf "1 var%d 0 %d A B\n", j, (j + 1) * 1000 > "tmp_data.bim"
  }
}'
cp tmp_data.fam tmp_miss.fam
cp tmp_data.bim tmp_miss.bim

# The genotypes come from a plain LCG so every awk produces the same file.  The
# low bits of an LCG have short periods, so values are taken from bits 9-30.
awk -v missrate=0 -v out=tmp_data.bed 'BEGIN {
  s = 5551212
  printf "%c%c%c", 108, 27, 1 > out
  for (j = 0; j < 300; ++j) {
    s = (1103515245 * s + 12345) % 2147483648
    q = 5 + (int(s / 512) % 90)
    for (base = 0; base < 60; base += 4) {
      byte = 0
      for (k = 0; k < 4; ++k) {
        s = (1103515245 * s + 12345) % 2147483648
        if (missrate > 0 && (int(s / 512) % 100) < missrate) {
          g = 1
        } else {
          a = 0
          s = (1103515245 * s + 12345) % 2147483648
          if ((int(s / 512) % 100) < q) { a += 1 }
          s = (1103515245 * s + 12345) % 2147483648
          if ((int(s / 512) % 100) < q) { a += 1 }
          g = (a == 0) ? 0 : ((a == 1) ? 2 : 3)
        }
        byte += g * (4 ^ k)
      }
      printf "%c", byte > out
    }
  }
}'
awk -v missrate=5 -v out=tmp_miss.bed 'BEGIN {
  s = 90210
  printf "%c%c%c", 108, 27, 1 > out
  for (j = 0; j < 300; ++j) {
    s = (1103515245 * s + 12345) % 2147483648
    q = 5 + (int(s / 512) % 90)
    for (base = 0; base < 60; base += 4) {
      byte = 0
      for (k = 0; k < 4; ++k) {
        s = (1103515245 * s + 12345) % 2147483648
        if (missrate > 0 && (int(s / 512) % 100) < missrate) {
          g = 1
        } else {
          a = 0
          s = (1103515245 * s + 12345) % 2147483648
          if ((int(s / 512) % 100) < q) { a += 1 }
          s = (1103515245 * s + 12345) % 2147483648
          if ((int(s / 512) % 100) < q) { a += 1 }
          g = (a == 0) ? 0 : ((a == 1) ? 2 : 3)
        }
        byte += g * (4 ^ k)
      }
      printf "%c", byte > out
    }
  }
}'

# Weight files: one with a header, some omitted variants and some zero weights;
# one without a header; one giving every variant weight 1.
awk 'BEGIN {
  s = 13579
  print "SNP WT" > "wts.txt"
  for (j = 0; j < 300; ++j) {
    # every variant gets weight 1 here, so that the weighted and unweighted
    # runs are comparing the same variant set
    printf "var%d 1\n", j > "ones.txt"
    printf "var%d %g\n", j, 1 + (j % 400) / 100.0 > "wts_nh.txt"
    s = (1103515245 * s + 12345) % 2147483648
    r = int(s / 512) % 10
    # wts.txt omits some variants but has no zero weight, so that the 1.9
    # differential below does not depend on the weight-0 fix in 1.9.
    if (r != 0) {
      s = (1103515245 * s + 12345) % 2147483648
      printf "var%d %g\n", j, 1 + (int(s / 512) % 400) / 100.0 > "wts.txt"
    }
    # zeros.txt is only used for 2.0-internal checks.
    printf "var%d %d\n", j, (r == 1) ? 0 : 1 > "zeros.txt"
    if (r == 1) {
      printf "var%d\n", j > "zero_ids.txt"
    }
  }
}'
sed -i.bak '1i\
SNP WT
' ones.txt

cmp_bin() {
    python3 - "$1" "$2" <<'PYEOF'
import numpy as np, sys
a = np.fromfile(sys.argv[1], dtype=np.float64)
b = np.fromfile(sys.argv[2], dtype=np.float64)
assert a.size == b.size and a.size, 'size %d vs %d' % (a.size, b.size)
den = np.maximum(np.abs(a), np.abs(b))
den[den == 0] = 1.0
r = float((np.abs(a - b) / den).max())
assert r < 1e-9, 'max relative difference %g' % r
PYEOF
}

# Differential against plink 1.9 on complete data, over both weight forms and
# all three report kinds.
for wts in "exp=0.5" "exp=1" "exp=-0.5" "wts.txt" "wts_nh.txt noheader"; do
    for rep in "" "ibs" "1-ibs"; do
        rm -f plink_dw.*.bin plink2_dw.*.bin
        plink --bfile tmp_data --distance square bin $rep --distance-wts $wts --out plink_dw
        $1/plink2 $2 $3 --bfile tmp_data --distance square bin $rep --distance-wts $wts --out plink2_dw
        cmp_bin $(ls plink_dw.*.bin) $(ls plink2_dw.*.bin)
    done
done

# With missing calls, against an independent implementation of the documented
# formula rather than against 1.9, whose weighted path disagrees with its own
# unweighted path there.
for wts in "exp=0.5" "wts.txt"; do
    $1/plink2 $2 $3 --bfile tmp_miss --distance square bin --distance-wts $wts --out plink2_dw
    python3 oracle.py tmp_miss "$wts" oracle.bin
    cmp_bin plink2_dw.dist.bin oracle.bin
done

# Zero-weight variants must drop out entirely: with every other weight equal to
# 1, the answer has to be the unweighted distance over the surviving variants.
sed -i.bak '1i\
SNP WT
' zeros.txt
$1/plink2 $2 $3 --bfile tmp_data --distance square bin --distance-wts zeros.txt --out plink2_zw
$1/plink2 $2 $3 --bfile tmp_data --exclude zero_ids.txt --distance square bin flat-missing --out plink2_zref
cmp plink2_zw.dist.bin plink2_zref.dist.bin

# A constant weight has to cancel out of both the distance and the missing-call
# correction, so the weighted and unweighted results must agree exactly.
$1/plink2 $2 $3 --bfile tmp_miss --distance square bin --out plink2_uw
$1/plink2 $2 $3 --bfile tmp_miss --distance square bin --distance-wts ones.txt --out plink2_ow
cmp plink2_uw.dist.bin plink2_ow.dist.bin

fails() {
    "$@" && false || true
}

# --distance-wts needs --distance, rejects 'flat-missing', and rejects a
# malformed exponent or a second parameter after 'exp='.
fails $1/plink2 $2 $3 --bfile tmp_data --distance-wts exp=0.5 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --distance square flat-missing --distance-wts exp=0.5 --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --distance square --distance-wts exp=notanumber --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --distance square --distance-wts exp=0.5 noheader --out plink2_bad
fails $1/plink2 $2 $3 --bfile tmp_data --distance square --distance-wts wts.txt nonsense --out plink2_bad

# A weight file with no usable entry, and one with a negative weight.
printf 'SNP WT\nnosuchvariant 1\n' > empty.txt
fails $1/plink2 $2 $3 --bfile tmp_data --distance square --distance-wts empty.txt --out plink2_bad
printf 'SNP WT\nvar0 -1\n' > neg.txt
fails $1/plink2 $2 $3 --bfile tmp_data --distance square --distance-wts neg.txt --out plink2_bad
