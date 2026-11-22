#!/bin/bash
# Usage: build.sh [plink2 build dir, default ../../build_dynamic]
# Builds ./plink2_fullprec (next to this script): identical to the plink2 binary in that
# directory except that --glm linear/logistic output and .psam phenotype output
# (used by the harness to recover the exact values plink2 parsed) have 17
# significant digits.
# The build directory must already contain an up-to-date plink2 build.
set -eo pipefail
here=$(cd "$(dirname "$0")" && pwd)
bdir=$(cd "${1:-$here/../../build_dynamic}" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cd "$bdir"
# Reuse the compile and link commands make would use.
dryrun=$(make -n -B 2>/dev/null)
flags=$(echo "$dryrun" | grep -- '-std=c++' | head -n 1 | sed -e 's/ \.\.\/[^ ]*\.cc//g' -e 's/ -c$//')
for f in plink2_glm_linear plink2_glm_logistic plink2_data; do
    eval "$flags -include $here/full_precision.h -Wno-unused-function ../$f.cc -c -o $tmp/$f.o"
done
objs=$(ls *.o | grep -v -e '^plink2_glm_linear\.o$' -e '^plink2_glm_logistic\.o$' -e '^plink2_data\.o$' | sed "s|^|$bdir/|" | tr "\n" " ")
linkcmd=$(echo "$dryrun" | grep -- '-o plink2' | head -n 1)
linkcmd=$(echo "$linkcmd" | sed -e 's/ [A-Za-z0-9_]*\.o//g' -e "s#-o plink2\"\"#-o $here/plink2_fullprec#" -e "s#-o plink2 #-o $here/plink2_fullprec #")
eval "$linkcmd $objs $tmp/plink2_glm_linear.o $tmp/plink2_glm_logistic.o $tmp/plink2_data.o"
echo "built $here/plink2_fullprec"
