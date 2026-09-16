# S3 and HTTP input support for PLINK 2.0

PLINK 2.0 can read input files directly from S3 (`s3://bucket/key`) and from
HTTP(S) URLs, including presigned URLs. Data is streamed on demand through
HTTP range requests, so no temporary copy is written to disk and a run that
touches only part of a file only transfers that part.

This is an optional feature. It is off by default and adds a single
dependency when enabled: **libcurl 7.75.0 or newer**.

## Building

```sh
# development build
make USE_S3=1

# production build
cd build_dynamic && make USE_S3=1
```

libcurl must be discoverable through the compiler's normal search paths.

| Platform | Command |
| --- | --- |
| Debian/Ubuntu | `apt install libcurl4-openssl-dev` |
| RHEL/Fedora | `dnf install libcurl-devel` |
| macOS | preinstalled, or `brew install curl` |
| conda/pixi | `pixi add libcurl` |

7.75.0 is required because request signing is done by libcurl's
`CURLOPT_AWS_SIGV4`. Anything from 2021 onward qualifies; notably, RHEL 8's
system curl (7.61) does not.

For the Python bindings:

```sh
cd Python && PGENLIB_USE_S3=1 pip install .
```

## Usage

Any input path may be an S3 URI:

```sh
plink2 --pfile s3://my-bucket/cohort/chr1 --freq --out chr1
plink2 --pfile local/chr1 --keep s3://my-bucket/lists/keep.txt --freq
plink2 --bfile s3://my-bucket/cohort/chr1 --make-pgen --out chr1
```

`--pfile`/`--bfile` expand to the usual triple of objects, so
`s3://bucket/cohort/chr1` reads `chr1.pgen`, `chr1.pvar` and `chr1.psam`.

Output always goes to local disk. There is no support for writing to S3.

### Public buckets

```sh
plink2 --pfile s3://1000genomes-dragen/data/chr1 --s3-no-sign-request --freq
```

`AWS_NO_SIGN_REQUEST=1` does the same thing. Unsigned access is also tried
automatically if no credentials can be found, or if the credentials that were
found are rejected.

### Presigned URLs

A presigned URL carries its own signature, so it is fetched unsigned and needs
no credentials:

```sh
url=$(aws s3 presign s3://my-bucket/lists/keep.txt --expires-in 3600)
plink2 --pfile local/chr1 --keep "$url" --freq
```

Because each presigned URL names exactly one object, they cannot be used with
`--pfile`/`--bfile`, which need a shared prefix. Use them for single-file
flags such as `--keep`, `--extract` or `--pheno`.

## Credentials

Sources are tried in this order, stopping at the first that yields a usable
key pair:

1. `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` / `AWS_SESSION_TOKEN`
2. `~/.aws/credentials` (or `AWS_SHARED_CREDENTIALS_FILE`)
3. `~/.aws/config` (or `AWS_CONFIG_FILE`)
4. ECS/EKS container endpoint (`AWS_CONTAINER_CREDENTIALS_RELATIVE_URI`
   or `AWS_CONTAINER_CREDENTIALS_FULL_URI`)
5. EC2 instance metadata, IMDSv2
6. unsigned

`AWS_PROFILE` selects the profile; non-default profiles may be written either
`[name]` or `[profile name]`. Container and instance-metadata credentials
carry an expiry and are re-read automatically before they lapse, so runs
longer than the token lifetime are fine.

### Environment variables

| Variable | Effect |
| --- | --- |
| `AWS_REGION`, `AWS_DEFAULT_REGION` | Signing region (default `us-east-1`) |
| `AWS_ENDPOINT_URL_S3`, `AWS_ENDPOINT_URL` | Custom endpoint, e.g. MinIO |
| `AWS_S3_ADDRESSING_STYLE` | `path` or `virtual` |
| `AWS_S3_FORCE_PATH_STYLE` | Equivalent to `AWS_S3_ADDRESSING_STYLE=path` |
| `AWS_NO_SIGN_REQUEST` | Send unsigned requests |
| `AWS_CA_BUNDLE` | Alternate CA certificate bundle |
| `AWS_EC2_METADATA_DISABLED` | Skip the IMDS probe |
| `HTTPS_PROXY`, `NO_PROXY` | Honored by libcurl |

### S3-compatible services

A custom endpoint implies path-style addressing, which is what MinIO, Ceph and
LocalStack expect:

```sh
export AWS_ENDPOINT_URL_S3=http://localhost:9000
export AWS_ACCESS_KEY_ID=minioadmin AWS_SECRET_ACCESS_KEY=minioadmin
plink2 --pfile s3://my-bucket/chr1 --freq
```

Set `AWS_S3_ADDRESSING_STYLE=virtual` to override. Against real S3, the
addressing style is chosen automatically: virtual-hosted, except for bucket
names that are not valid DNS labels — notably names containing dots, which
fall outside AWS's wildcard certificate.

## Per-file credentials in Python

Environment variables are process-global, so they cannot express "these two
files live in different accounts". The Python API therefore accepts a
`UPath` whose credentials are attached to the object, and applies them to that
file alone:

```python
from upath import UPath
import pgenlib

a = pgenlib.PgenReader(UPath("s3://bucket-a/chr1.pgen",
                             key="AKIA...", secret="..."))
b = pgenlib.PgenReader(UPath("s3://bucket-b/chr1.pgen",
                             key="AKIB...", secret="...",
                             endpoint_url="https://minio.internal"))
```

Recognized `storage_options` keys: `key`, `secret`, `token`, `endpoint_url`,
`anon`, and `client_kwargs={'region_name': ...}`. Plain `str`/`bytes` paths
still use the ambient credential chain.

## Behavior worth knowing

**Consistency.** The object's ETag is recorded when it is opened and every
later range request is conditioned on it. If the object is overwritten
mid-read, the run fails with an explicit error rather than silently splicing
together two versions of the file.

**Truncation.** A range response shorter than requested is an error, not an
early EOF, so a partial transfer cannot be mistaken for a short file.

**Retries.** Throttling (`503 SlowDown`), transient 5xx responses and dropped
connections are retried five times with exponential backoff and full jitter.

**Wrong region.** If a bucket lives in another region, the `x-amz-bucket-region`
hint in the error response is followed automatically.

**Reads are sequential-friendly.** Data is fetched in 8 MiB chunks. `fseek`
within the current chunk is free; seeking outside it costs one request.

## Deliberate limitations

Supporting these would require substantially more machinery than the
workaround costs, so they are out of scope:

- **SSO** (`aws sso login`), **`credential_process`**, **`role_arn`/
  `source_profile` AssumeRole chaining**, and **web identity / IRSA**.
  For all of these, obtain credentials with the AWS CLI and export them:

  ```sh
  eval "$(aws configure export-credentials --format env)"
  ```

  On EKS specifically, the node's instance role is picked up via IMDSv2
  without any extra step, as long as the pod can reach the metadata service.

- **Writing to S3.** Output is always local.
- **S3 Access Points, Multi-Region Access Points, Outposts, S3 Express One
  Zone, dualstack and FIPS endpoints.**
- **SSE-C.** Server-side encryption with S3- or KMS-managed keys is
  transparent and works normally.
- **Requester-pays buckets.**

## Implementation

The S3 machinery lives in [`s3stream/`](s3stream/), which has no plink2
dependencies and can be reused as-is in other projects. `plink2_s3.cc` is only
an adapter that maps plink2's types onto it.

| File | Contents |
| --- | --- |
| `s3stream/s3stream.h` | Public API |
| `s3stream/s3stream.cc` | Endpoint building, streaming, `FILE*` adapters |
| `s3stream/s3stream_creds.cc` | Credential chain |
| `s3stream/s3stream_http.cc` | Requests, signing, retries |

The `FILE*` returned by `s3stream_open()` is backed by `fopencookie` on Linux,
`funopen` on macOS/BSD, and a temporary file on Windows, which has neither.

Build without `S3STREAM_ENABLE` and the module compiles to stubs that reject
remote paths with a clear message, so callers need no `#ifdef`s.

## Tests

```sh
cd Tests/TEST_S3 && ./run_tests.sh ../../build_dynamic
```

24 tests covering every authentication mode, addressing style, presigned URLs,
multi-chunk reads, zero-byte objects and the error paths, against a local
MinIO started via Docker (or natively with `S3_TEST_NATIVE=1`). Two additional
tests run against real S3 when `PLINK_TESTS_S3_SUPPORT=1` is set.
