#!/bin/bash

set -e

# ---------- Cleanup ----------
rm -rf test_data derivatives
make -C 2.0/build_dynamic/ clean 