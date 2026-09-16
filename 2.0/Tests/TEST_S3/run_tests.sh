#!/bin/bash
# S3 input tests for plink2, run against a local MinIO server.
#
# Usage: ./run_tests.sh [dir containing a plink2 built with USE_S3=1]
#
# Runs MinIO in docker when available, otherwise downloads the native minio/mc
# binaries (GitHub's macOS runners have no docker).  Force the latter with
# S3_TEST_NATIVE=1.  The session-token test additionally needs the AWS CLI.
# Set KEEP_MINIO=1 to leave the server running afterwards.

set -euo pipefail

if [[ $# -eq 0 ]]; then
    BUILD_DIR=../../build_dynamic
else
    BUILD_DIR=$1
fi
BUILD_DIR=$(cd "$BUILD_DIR" && pwd)
PLINK2="$BUILD_DIR/plink2"

if [[ ! -x "$PLINK2" ]]; then
    echo "Error: $PLINK2 not found or not executable." 1>&2
    exit 1
fi

MINIO_IMAGE=${MINIO_IMAGE:-quay.io/minio/minio:latest}
MINIO_CONTAINER=${MINIO_CONTAINER:-plink2-s3-test-minio}
MINIO_PORT=${MINIO_PORT:-9000}
MINIO_ROOT_USER=${MINIO_ROOT_USER:-minioadmin}
MINIO_ROOT_PASSWORD=${MINIO_ROOT_PASSWORD:-minioadmin}
# MINIO_DOMAIN makes MinIO accept virtual-hosted-style requests
# (<bucket>.$MINIO_DOMAIN).  plink2 defaults to path-style against a custom
# endpoint, so this is only needed by the virtual-host addressing test.
MINIO_DOMAIN=${MINIO_DOMAIN:-localhost}

PRIVATE_BUCKET=plink2-private
PUBLIC_BUCKET=plink2-public
RO_USER=plink2ro
RO_SECRET=plink2rosecret12

# Endpoint plink2 talks to.  A custom endpoint implies path-style addressing,
# so an IP literal works and no hostname needs to resolve.
S3_ENDPOINT="http://127.0.0.1:${MINIO_PORT}"
# Endpoint the admin tooling talks to (path-style, inside the container).
ADMIN_ENDPOINT="http://127.0.0.1:${MINIO_PORT}"

TEST_DIR=$(pwd)
WORK=$TEST_DIR/tmp_s3
rm -rf "$WORK"
mkdir -p "$WORK"

PASS_CT=0

pass() {
    PASS_CT=$((PASS_CT + 1))
    echo "  ok: $1"
}

fail() {
    echo "  FAILED: $1" 1>&2
    exit 1
}

###########################################################################
echo "=== public bucket on real S3 ==="
###########################################################################

# Reads a public nf-core dataset in eu-west-1, which is the only part of the
# suite that exercises real SigV4 over TLS and a non-default region.  Off by
# default because it needs network access.
if [[ -n "${PLINK_TESTS_S3_SUPPORT:-}" ]]; then
    S3_PREFIX="s3://ngi-igenomes/testdata/nf-core/modules/genomics/homo_sapiens/popgen/plink_simulated"
    cd "$WORK"
    for fileset in pfile bfile; do
        AWS_REGION="${AWS_REGION:-eu-west-1}" AWS_EC2_METADATA_DISABLED=true \
            "$PLINK2" ${2:-} ${3:-} "--$fileset" "$S3_PREFIX" --s3-no-sign-request \
            --geno-counts --out "real_$fileset" --silent > /dev/null ||
            fail "real S3 --$fileset: plink2 exited nonzero"
        diff -q "$TEST_DIR/expected/$fileset.gcount" "real_$fileset.gcount" > /dev/null ||
            fail "real S3 --$fileset: output mismatch"
        pass "real S3 public bucket, --$fileset"
    done
    cd "$TEST_DIR"
else
    echo "  skipped: set PLINK_TESTS_S3_SUPPORT=1 to test against real S3"
fi

if [[ -n "${S3_TEST_NATIVE:-}" ]] || ! command -v docker > /dev/null 2>&1; then
    MODE=native
else
    MODE=docker
fi

# The MinIO half needs either a working docker daemon or native minio/mc.
if [[ $MODE == docker ]] && ! docker info > /dev/null 2>&1; then
    MODE=native
fi
if [[ $MODE == native ]] && { [[ -z "$(type -P minio)" ]] || [[ -z "$(type -P mc)" ]]; }; then
    echo "=== local MinIO tests skipped (no docker daemon, and minio/mc not on PATH) ==="
    echo
    echo "$PASS_CT S3 tests passed."
    exit 0
fi

mc() {
    if [[ $MODE == docker ]]; then
        docker exec "$MINIO_CONTAINER" mc "$@"
    else
        "$MC_BIN" --config-dir "$WORK/mc-config" "$@"
    fi
}

# Makes a local file visible to mc and echoes the path mc should read it from.
stage() {
    if [[ $MODE == docker ]]; then
        docker cp "$WORK/$1" "$MINIO_CONTAINER:/tmp/$1" > /dev/null
        echo "/tmp/$1"
    else
        echo "$WORK/$1"
    fi
}

cleanup() {
    if [[ -n "${KEEP_MINIO:-}" ]]; then
        return
    fi
    if [[ $MODE == docker ]]; then
        docker rm -f "$MINIO_CONTAINER" > /dev/null 2>&1 || true
    elif [[ -n "${MINIO_PID:-}" ]]; then
        kill "$MINIO_PID" > /dev/null 2>&1 || true
    fi
}
trap cleanup EXIT

###########################################################################
echo "=== starting MinIO ($MODE) ==="
###########################################################################

if [[ $MODE == docker ]]; then
    docker rm -f "$MINIO_CONTAINER" > /dev/null 2>&1 || true
    docker run -d --name "$MINIO_CONTAINER" \
        -p "${MINIO_PORT}:9000" \
        -e "MINIO_ROOT_USER=${MINIO_ROOT_USER}" \
        -e "MINIO_ROOT_PASSWORD=${MINIO_ROOT_PASSWORD}" \
        -e "MINIO_DOMAIN=${MINIO_DOMAIN}" \
        "$MINIO_IMAGE" server /data > /dev/null
else
    # MinIO stopped serving standalone binaries from dl.min.io (HTTP 410), so
    # they have to come from the environment: `brew install minio minio-mc` or
    # conda-forge's minio-server / minio-client.
    # type -P searches PATH only; `command -v mc` would find the mc() wrapper
    # defined above and recurse.
    MINIO_BIN=${MINIO_BIN:-$(type -P minio || true)}
    MC_BIN=${MC_BIN:-$(type -P mc || true)}
    if [[ -z "$MINIO_BIN" || -z "$MC_BIN" ]]; then
        fail "native mode needs 'minio' and 'mc' on PATH (brew install minio minio-mc)"
    fi
    mkdir -p "$WORK/minio-data" "$WORK/mc-config"
    MINIO_ROOT_USER="$MINIO_ROOT_USER" \
    MINIO_ROOT_PASSWORD="$MINIO_ROOT_PASSWORD" \
    MINIO_DOMAIN="$MINIO_DOMAIN" \
        "$MINIO_BIN" server --address ":${MINIO_PORT}" "$WORK/minio-data" \
        > "$WORK/minio.log" 2>&1 &
    MINIO_PID=$!
fi

for _ in $(seq 1 60); do
    if curl -fs -o /dev/null "${ADMIN_ENDPOINT}/minio/health/live"; then
        break
    fi
    # Bail out immediately if the server we started has already died, rather
    # than waiting for the timeout or, worse, provisioning a foreign server
    # that happens to hold the port.
    if [[ $MODE == native ]] && ! kill -0 "$MINIO_PID" 2> /dev/null; then
        echo "--- minio.log ---" 1>&2
        cat "$WORK/minio.log" 1>&2 || true
        fail "MinIO exited during startup"
    fi
    sleep 1
done
curl -fs -o /dev/null "${ADMIN_ENDPOINT}/minio/health/live" || fail "MinIO did not become healthy"
if [[ $MODE == native ]] && ! kill -0 "$MINIO_PID" 2> /dev/null; then
    fail "port ${MINIO_PORT} is served by another process; set MINIO_PORT"
fi

# Virtual-hosted-style addressing needs <bucket>.$MINIO_DOMAIN to resolve.
# It is no longer the default, so a name that does not resolve just skips the
# one test that exercises it.  Any HTTP status means the name resolved and
# MinIO answered; 000 means it did not.
VHOST_OK=1
for b in "$PRIVATE_BUCKET" "$PUBLIC_BUCKET"; do
    code=$(curl -s -m 10 -o /dev/null -w '%{http_code}' \
        "http://${b}.${MINIO_DOMAIN}:${MINIO_PORT}/" || true)
    if [[ "$code" == "000" ]]; then
        VHOST_OK=0
    fi
done

###########################################################################
echo "=== generating test data ==="
###########################################################################

cd "$WORK"
"$PLINK2" --dummy 200 2000 --make-pgen --out ref --silent > /dev/null
# ~10 MB .pgen, which forces more than one 8 MiB S3 range request.
"$PLINK2" --dummy 1000 40000 --make-pgen --out big --silent > /dev/null

# Local reference outputs to compare the S3 runs against.
"$PLINK2" --pfile ref --freq --out local_ref --silent > /dev/null
"$PLINK2" --pfile big --freq --out local_big --silent > /dev/null

###########################################################################
echo "=== provisioning MinIO ==="
###########################################################################

mc alias set test "$ADMIN_ENDPOINT" "$MINIO_ROOT_USER" "$MINIO_ROOT_PASSWORD" > /dev/null
mc mb --ignore-existing "test/$PRIVATE_BUCKET" > /dev/null
mc mb --ignore-existing "test/$PUBLIC_BUCKET" > /dev/null

# Fixtures for the non-.pgen input flags and for the zero-byte object case.
head -51 ref.psam > keep50.txt
: > empty.txt
upload_paths=()
for f in ref.pgen ref.pvar ref.psam big.pgen big.pvar big.psam keep50.txt empty.txt; do
    upload_paths+=("$(stage "$f")")
done
for b in "$PRIVATE_BUCKET" "$PUBLIC_BUCKET"; do
    mc cp "${upload_paths[@]}" "test/$b/data/" > /dev/null
done

# Least-privilege user: MinIO's builtin "readonly" policy grants s3:GetObject
# but not s3:ListBucket, matching a typical locked-down AWS setup.
mc admin user add test "$RO_USER" "$RO_SECRET" > /dev/null
mc admin policy attach test readonly --user "$RO_USER" > /dev/null
mc anonymous set download "test/$PUBLIC_BUCKET" > /dev/null

###########################################################################
# Test harness
###########################################################################

# Every test runs plink2 with a clean AWS environment, so that only the
# credentials the test explicitly provides are visible.  HOME points at an
# empty directory unless a test overrides it, so that ~/.aws is never picked up
# by accident (later assignments win in env(1)).
plink2_s3() {
    env -u AWS_ACCESS_KEY_ID -u AWS_SECRET_ACCESS_KEY -u AWS_SESSION_TOKEN \
        -u AWS_PROFILE -u AWS_DEFAULT_PROFILE -u AWS_SHARED_CREDENTIALS_FILE \
        -u AWS_CONFIG_FILE -u AWS_ENDPOINT_URL -u AWS_ENDPOINT_URL_S3 \
        AWS_ENDPOINT_URL_S3="$S3_ENDPOINT" \
        AWS_REGION=us-east-1 \
        AWS_EC2_METADATA_DISABLED=true \
        HOME="$WORK/emptyhome" \
        "$@"
}

run_freq() {
    # usage: run_freq <label> <out prefix> <pfile prefix> [env assignments...]
    local label=$1 out=$2 pfile=$3
    shift 3
    plink2_s3 "$@" "$PLINK2" --pfile "$pfile" --freq --out "$out" --silent > /dev/null ||
        fail "$label: plink2 exited nonzero"
    diff -q local_ref.afreq "$out.afreq" > /dev/null || fail "$label: output mismatch"
    pass "$label"
}

expect_failure() {
    # usage: expect_failure <label> <out prefix> <pfile prefix> [env assignments...]
    local label=$1 out=$2 pfile=$3
    shift 3
    set +e
    plink2_s3 "$@" "$PLINK2" --pfile "$pfile" --freq --out "$out" \
        > "$out.stdout" 2> "$out.stderr"
    local rc=$?
    set -e
    if [[ $rc -eq 0 ]]; then
        fail "$label: expected nonzero exit, got 0"
    fi
    pass "$label (exit $rc)"
}

mkdir -p "$WORK/emptyhome" "$WORK/fakehome"

###########################################################################
echo "=== auth: environment variable credentials ==="
###########################################################################

run_freq "env-var credentials" env_creds "s3://$PRIVATE_BUCKET/data/ref" \
    "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET"

###########################################################################
echo "=== auth: temporary credentials with session token ==="
###########################################################################

if command -v aws > /dev/null; then
    STS_JSON=$(AWS_ACCESS_KEY_ID="$RO_USER" AWS_SECRET_ACCESS_KEY="$RO_SECRET" \
        AWS_DEFAULT_REGION=us-east-1 \
        aws --endpoint-url "$ADMIN_ENDPOINT" sts assume-role \
            --role-arn arn:aws:iam::123456789012:role/plink2-test \
            --role-session-name plink2-test \
            --duration-seconds 900 --output json)
    STS_KEY=$(printf '%s' "$STS_JSON" | sed -n 's/.*"AccessKeyId": "\([^"]*\)".*/\1/p')
    STS_SECRET=$(printf '%s' "$STS_JSON" | sed -n 's/.*"SecretAccessKey": "\([^"]*\)".*/\1/p')
    STS_TOKEN=$(printf '%s' "$STS_JSON" | sed -n 's/.*"SessionToken": "\([^"]*\)".*/\1/p')
    [[ -n "$STS_TOKEN" ]] || fail "could not obtain STS session token"
    run_freq "session-token credentials" sts_creds "s3://$PRIVATE_BUCKET/data/ref" \
        "AWS_ACCESS_KEY_ID=$STS_KEY" "AWS_SECRET_ACCESS_KEY=$STS_SECRET" \
        "AWS_SESSION_TOKEN=$STS_TOKEN"
else
    echo "  skipped: session-token test (no aws CLI)"
fi

###########################################################################
echo "=== auth: shared credentials file ==="
###########################################################################

mkdir -p "$WORK/fakehome/.aws"
cat > "$WORK/fakehome/.aws/credentials" <<EOF
[default]
aws_access_key_id = $RO_USER
aws_secret_access_key = $RO_SECRET

[plink2-named]
aws_access_key_id = $RO_USER
aws_secret_access_key = $RO_SECRET

[plink2-bogus]
aws_access_key_id = AKIAIOSFODNN7EXAMPLE
aws_secret_access_key = wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY
EOF

run_freq "default profile in ~/.aws/credentials" default_profile \
    "s3://$PRIVATE_BUCKET/data/ref" "HOME=$WORK/fakehome"

run_freq "named profile via AWS_PROFILE" named_profile \
    "s3://$PRIVATE_BUCKET/data/ref" "HOME=$WORK/fakehome" "AWS_PROFILE=plink2-named"

cp "$WORK/fakehome/.aws/credentials" "$WORK/alt_credentials"
run_freq "AWS_SHARED_CREDENTIALS_FILE" alt_creds_file \
    "s3://$PRIVATE_BUCKET/data/ref" \
    "AWS_SHARED_CREDENTIALS_FILE=$WORK/alt_credentials" \
    "AWS_PROFILE=plink2-named"

###########################################################################
echo "=== auth: anonymous ==="
###########################################################################

run_freq "anonymous, no credentials anywhere" anon "s3://$PUBLIC_BUCKET/data/ref"

run_freq "anonymous fallback with unusable credentials" anon_fallback \
    "s3://$PUBLIC_BUCKET/data/ref" \
    "AWS_ACCESS_KEY_ID=AKIAIOSFODNN7EXAMPLE" \
    "AWS_SECRET_ACCESS_KEY=wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY"

###########################################################################
echo "=== --s3-no-sign-request ==="
###########################################################################

plink2_s3 "$PLINK2" --pfile "s3://$PUBLIC_BUCKET/data/ref" --s3-no-sign-request \
    --freq --out nosign --silent > /dev/null ||
    fail "--s3-no-sign-request: plink2 exited nonzero"
diff -q local_ref.afreq nosign.afreq > /dev/null || fail "--s3-no-sign-request: output mismatch"
pass "--s3-no-sign-request on a public bucket"

run_freq "AWS_NO_SIGN_REQUEST=1 on a public bucket" nosign_env \
    "s3://$PUBLIC_BUCKET/data/ref" "AWS_NO_SIGN_REQUEST=1"

# Proves signing really is skipped: these credentials would otherwise work.
set +e
plink2_s3 "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "$PLINK2" --pfile "s3://$PRIVATE_BUCKET/data/ref" --s3-no-sign-request \
    --freq --out nosign_private > nosign_private.stdout 2> nosign_private.stderr
nosign_rc=$?
set -e
if [[ $nosign_rc -eq 0 ]]; then
    fail "--s3-no-sign-request still used the supplied credentials"
fi
pass "--s3-no-sign-request ignores valid credentials (exit $nosign_rc)"


###########################################################################
echo "=== mixed local and S3 inputs ==="
###########################################################################

plink2_s3 "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "$PLINK2" --pgen "s3://$PRIVATE_BUCKET/data/ref.pgen" \
    --pvar ref.pvar --psam ref.psam --freq --out mixed --silent > /dev/null ||
    fail "mixed local/S3 inputs: plink2 exited nonzero"
diff -q local_ref.afreq mixed.afreq > /dev/null || fail "mixed local/S3 inputs: output mismatch"
pass "mixed local and S3 inputs"

###########################################################################
echo "=== multi-chunk range reads ==="
###########################################################################

plink2_s3 "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "$PLINK2" --pfile "s3://$PRIVATE_BUCKET/data/big" --freq --out s3_big --silent > /dev/null ||
    fail "multi-chunk read: plink2 exited nonzero"
diff -q local_big.afreq s3_big.afreq > /dev/null || fail "multi-chunk read: output mismatch"
pass "multi-chunk range reads (10 MB .pgen)"

plink2_s3 "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "$PLINK2" --pfile "s3://$PRIVATE_BUCKET/data/big" --make-pgen --out s3_roundtrip --silent > /dev/null ||
    fail "round trip: plink2 exited nonzero"
cmp big.pgen s3_roundtrip.pgen || fail "round trip: .pgen differs"
pass "byte-identical .pgen round trip from S3"

# --make-bed and --make-pgen run realpath() on the input to detect an
# input/output collision, which cannot succeed for an s3:// URI.
"$PLINK2" --pfile ref --make-bed --out local_bed --silent > /dev/null
plink2_s3 "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "$PLINK2" --pfile "s3://$PRIVATE_BUCKET/data/ref" --make-bed --out s3_bed --silent > /dev/null ||
    fail "--make-bed from S3: plink2 exited nonzero"
cmp local_bed.bed s3_bed.bed || fail "--make-bed from S3: .bed differs"
pass "--make-bed from S3 input"

###########################################################################
echo "=== path-style addressing ==="
###########################################################################

# Virtual-hosted style cannot work against an IP literal, so this only passes
# if AWS_S3_FORCE_PATH_STYLE is honored.
run_freq "path-style against an IP endpoint" path_style "s3://$PRIVATE_BUCKET/data/ref" \
    "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "AWS_ENDPOINT_URL_S3=$ADMIN_ENDPOINT" "AWS_S3_FORCE_PATH_STYLE=1"

run_freq "AWS_S3_ADDRESSING_STYLE=path" addressing_style "s3://$PRIVATE_BUCKET/data/ref" \
    "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "AWS_ENDPOINT_URL_S3=$ADMIN_ENDPOINT" "AWS_S3_ADDRESSING_STYLE=path"

# Path-style is the default for a custom endpoint, so virtual-hosted style has
# to be asked for explicitly.  Needs <bucket>.$MINIO_DOMAIN to resolve.
if [[ $VHOST_OK == 1 ]]; then
    run_freq "AWS_S3_ADDRESSING_STYLE=virtual" virtual_style "s3://$PRIVATE_BUCKET/data/ref" \
        "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
        "AWS_ENDPOINT_URL_S3=http://${MINIO_DOMAIN}:${MINIO_PORT}" \
        "AWS_S3_ADDRESSING_STYLE=virtual"
else
    echo "skipping virtual-host addressing: *.${MINIO_DOMAIN} does not resolve"
fi

###########################################################################
echo "=== presigned https:// URLs ==="
###########################################################################

# A presigned URL carries its own signature, so it must be fetched unsigned
# and must work with no credentials in the environment at all.
PRESIGNED=$(AWS_ACCESS_KEY_ID="$MINIO_ROOT_USER" \
    AWS_SECRET_ACCESS_KEY="$MINIO_ROOT_PASSWORD" \
    aws s3 presign "s3://$PRIVATE_BUCKET/data/keep50.txt" \
    --endpoint-url "$ADMIN_ENDPOINT" --expires-in 3600)
[[ -n "$PRESIGNED" ]] || fail "could not generate a presigned URL"

"$PLINK2" --pfile ref --keep keep50.txt --freq --out local_presign --silent > /dev/null
env -u AWS_ACCESS_KEY_ID -u AWS_SECRET_ACCESS_KEY -u AWS_SESSION_TOKEN \
    -u AWS_PROFILE -u AWS_SHARED_CREDENTIALS_FILE \
    "$PLINK2" --pfile ref --keep "$PRESIGNED" \
    --freq --out presign_keep --silent > /dev/null ||
    fail "presigned URL: plink2 exited nonzero"
cmp local_presign.afreq presign_keep.afreq || fail "presigned URL: .afreq differs"
pass "presigned https:// URL as an input file"

###########################################################################
echo "=== S3 inputs other than .pgen/.pvar/.psam ==="
###########################################################################

"$PLINK2" --pfile ref --keep keep50.txt --freq --out local_keep --silent > /dev/null
plink2_s3 "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "$PLINK2" --pfile ref --keep "s3://$PRIVATE_BUCKET/data/keep50.txt" \
    --freq --out s3_keep --silent > /dev/null ||
    fail "--keep from S3: plink2 exited nonzero"
diff -q local_keep.afreq s3_keep.afreq > /dev/null || fail "--keep from S3: output mismatch"
pass "--keep reads an s3:// file with local .pgen input"

###########################################################################
echo "=== zero-byte object ==="
###########################################################################

# "Range: bytes=0-0" is unsatisfiable on an empty object, so without explicit
# 416 handling this fails in the S3 layer instead of being read as empty.
set +e
plink2_s3 "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET" \
    "$PLINK2" --pfile ref --keep "s3://$PRIVATE_BUCKET/data/empty.txt" \
    --freq --out s3_empty > empty.stdout 2> empty.stderr
set -e
if grep -q "Cannot access" empty.stdout empty.stderr; then
    fail "zero-byte object: reported as inaccessible instead of empty"
fi
pass "zero-byte object opens and reads as empty"

###########################################################################
echo "=== error handling ==="
###########################################################################

expect_failure "nonexistent object" no_such_key "s3://$PRIVATE_BUCKET/data/nope" \
    "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET"

expect_failure "nonexistent bucket" no_such_bucket "s3://plink2-does-not-exist/data/ref" \
    "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=$RO_SECRET"

expect_failure "private bucket without credentials" no_creds "s3://$PRIVATE_BUCKET/data/ref"

expect_failure "wrong secret key" bad_secret "s3://$PRIVATE_BUCKET/data/ref" \
    "AWS_ACCESS_KEY_ID=$RO_USER" "AWS_SECRET_ACCESS_KEY=definitely-not-the-secret"

echo
echo "All $PASS_CT S3 tests passed."
