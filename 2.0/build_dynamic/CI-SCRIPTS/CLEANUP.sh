#!/bin/bash

set -e

# ---------- Cleanup ----------
rm -rf test_data results
make -C 2.0/build_dynamic/ clean 