#!/bin/bash

# Multipass non-concatenating --pmerge-list.  The hard requirement is that the
# output matches the corresponding single-pass merge byte for byte, so this
# compares the two on randomly generated collections of overlapping filesets;
# see multipass_oracle.py.  '--debug pmerge-pass-size=' forces small
# passes, and also makes single-pass merges of more than 20 filesets
# possible as a reference.

set -exo pipefail

# The weekly workflow (auto_checkcommits_tests_linuxmac.yml) sets
# PLINK2_FULL_TESTS=1 for 60 random configurations; otherwise the first 10
# are run (including two with more than 20 filesets, since every fifth
# configuration has that many).
config_ct=10
[ -n "$PLINK2_FULL_TESTS" ] && config_ct=60
python3 multipass_oracle.py $1/plink2 $config_ct 1 $2 $3
