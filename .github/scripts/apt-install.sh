#!/bin/bash
# Usage: apt-install.sh <package>...
#
# Runs "apt-get update" and "apt-get install -y" for the given packages,
# retrying when the package mirrors are unreachable.  GitHub's Ubuntu runners
# sometimes lose DNS resolution for archive.ubuntu.com for minutes at a time
# ("Temporary failure resolving ..."), which used to fail whole test runs
# before anything was built.
#
# - Acquire::Retries makes apt retry each failed download a few times, with a
#   short backoff, which covers brief blips.
# - The outer loop retries the whole update + install after a pause, which
#   covers longer outages.  (apt-get update exits 0 even when it cannot reach
#   a mirror, so the install has to be part of the retried command.)
# The job-level timeout-minutes in each workflow bounds the total time.

set -uo pipefail

attempts=3
pause_seconds=60
apt_opts=(-o Acquire::Retries=5 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30)

for attempt in $(seq 1 "$attempts"); do
    if sudo apt-get "${apt_opts[@]}" update &&
       sudo apt-get "${apt_opts[@]}" install -y "$@"; then
        exit 0
    fi
    if [[ $attempt -lt $attempts ]]; then
        echo "apt-install.sh: attempt $attempt of $attempts failed; retrying in ${pause_seconds}s." >&2
        sleep "$pause_seconds"
    fi
done
echo "apt-install.sh: giving up after $attempts attempts." >&2
exit 1
