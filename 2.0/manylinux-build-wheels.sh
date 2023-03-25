#!/bin/bash
set -exuo pipefail
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r Python/dev-requirements.txt
    "${PYBIN}/pip" wheel Python/ --no-deps -w wheelhouse/
done
