#!/bin/bash

# Multipass non-concatenating --pmerge-list.  The hard requirement is that the
# output matches the corresponding single-pass merge byte for byte, so this
# compares the two on randomly generated collections of overlapping filesets;
# see multipass_oracle.py.  The undocumented --pmerge-pass-size flag forces
# small passes, and also makes single-pass merges of more than 20 filesets
# possible as a reference.

set -exo pipefail

python3 multipass_oracle.py $1/plink2 60 1 $2 $3
