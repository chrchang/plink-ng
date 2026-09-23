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

**Truncation.** Every range response must be exactly as long as the range
that was asked for. A short body is reported as a truncated transfer rather
than an early EOF, so a partial transfer cannot be mistaken for a short file.
A response *longer* than the range, or a `200` where a `206` was expected,
means the server ignored `Range`; the transfer is aborted at the requested
length and the read fails, because those bytes would otherwise land at the
wrong file offset.

**Retries.** Throttling (`503 SlowDown`), transient 5xx responses and dropped
connections are retried five times with exponential backoff and full jitter.

**Wrong region.** If a bucket lives in another region, the `x-amz-bucket-region`
hint in the error response is followed automatically.

**Reads are sequential-friendly.** Data is fetched in 8 MiB chunks. `fseek`
within the current chunk is free; seeking outside it costs one request.

## Security

S3 support is built so that it does not change what plink2 trusts, and adds
as little new attack surface as possible. Reading a file from S3 carries the
same risk as reading that file from local disk, from NFS, or from a bucket
mounted through FUSE (s3fs, mountpoint-s3). The default binary is not
affected at all.

**The default build is unchanged.** Without `USE_S3=1`, no network code is
compiled in and libcurl is not linked. The S3 entry points become stubs that
reject remote paths. Local paths go straight to `fopen` in both builds.

**plink2's parsers never touch the network.** The `.pgen`/`.pvar`/`.psam`
readers get an ordinary `FILE*` and see the same byte stream a local file
would give them. A malicious or corrupted file is exactly as dangerous from S3
as from NFS or a FUSE mount: the trust boundary is the file's contents, and
that boundary stays where it is. The stream layer also makes sure the bytes
are the right ones:

- Every range response must be exactly the length requested. Short responses,
  oversized responses and responses that ignore `Range` are hard errors, so a
  parser never sees data at the wrong offset or a silently truncated file.
- Every read is pinned to the ETag recorded at open, so the object cannot
  change underneath a running job. That is a stronger guarantee than NFS or
  FUSE give, since there a file can be rewritten between reading its header
  and reading its records.

**The network-facing code is small, isolated and bounded.**

- It lives in [`s3stream/`](s3stream/), includes no plink2 headers, and uses
  `std::string`/`std::map` rather than plink2's manual buffer arithmetic.
- TLS, HTTP parsing and SigV4 signing are handled by libcurl. htslib
  (samtools, bcftools) does its S3/HTTP support the same way, and s3fs-fuse
  is built on the same library.
- Every server response has a fixed memory bound. Data bodies are capped at
  the requested range (at most 8 MiB), only the four response headers
  s3stream uses are kept, and credential-metadata responses are capped at
  64 KiB. A hostile or misbehaving server cannot make plink2 allocate without
  limit.
- Server-supplied strings that are sent back later (ETag, region redirect
  hint, IMDS token and role name) are validated first. A response therefore
  cannot inject headers or send a signed request to another host.
- TLS certificate verification is set explicitly, only `http`/`https` are
  permitted, and HTTP redirects are never followed.
- Nothing is written to disk: no temporary files and no credential caching.
  Requests are signed with SigV4, so the secret key never leaves the process.
  Secrets held in memory are overwritten when released (best-effort).

**Things to keep in mind.**

- A presigned URL works as a bearer credential until it expires. Like any
  other argument, it shows up in the command line that plink2 echoes into its
  `.log`.
- A plain-`http://` custom endpoint (such as a local MinIO) is allowed, but
  then the transfer is not encrypted. Use `https://` for anything that leaves
  the machine.

## Windows is not supported

**S3/HTTP streaming does not work on Windows, full stop, and there is no
planned workaround.** `s3stream_open()` always fails on Windows builds with a
clear error message; local files are unaffected.

The reason is structural, not a missing feature: reads need to be served
lazily and support seeking, which on Linux and macOS/BSD is done by handing
libc a `FILE*` backed by our own read/seek callbacks (`fopencookie` /
`funopen`). MSVCRT/UCRT has no equivalent hook — its `FILE` struct is opaque
and ABI-locked, so there is no supported way to intercept `fread`/`fseek`
with custom logic while still returning a real `FILE*`. The only way to get a
seekable local handle without that hook is to download the entire remote
object to a temporary file before returning it, which a prior version of this
code did. That workaround has been removed: silently turning a routine
`--pfile s3://...` invocation into an unbounded, unannounced download —
potentially far larger than the user expects, with no way to know until disk
space or bandwidth runs out — is not an acceptable default. Genuinely lazy
streaming on Windows would require either depending on a POSIX-emulation
runtime (Cygwin/MSYS2's `newlib`-based CRT, which does support `funopen`, but
pulls in `cygwin1.dll`/`msys-2.0.dll` and drops the dependency-free native
`.exe`) or a virtual-filesystem layer (Windows ProjFS/Cloud Filter API, the
mechanism OneDrive uses for placeholder files) — both far outside the scope
of this change.

## Not planned

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

- **S3 Access Points, Multi-Region Access Points, Outposts, S3 Express One
  Zone, dualstack and FIPS endpoints.**

## Not yet implemented

No architectural blocker; these are just missing plumbing:

- **Writing to S3.** Output is always local. May be added in the future.
- **Requester-pays buckets.** Needs an `x-amz-request-payer: requester`
  header on every request.
- **SSE-C.** Server-side encryption with S3- or KMS-managed keys is
  transparent and works normally; only customer-supplied keys, which require
  sending the key material on each request, are unsupported.
- **`ExpectedBucketOwner`.**

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

The `FILE*` returned by `s3stream_open()` is backed by `fopencookie` on Linux
and `funopen` on macOS/BSD. On Windows it always returns `nullptr` with an
error — see [Windows is not supported](#windows-is-not-supported).

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
