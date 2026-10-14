# syntax=docker/dockerfile:1
#
# Container images for PLINK 2.0 and PLINK 1.9.
#
#   docker build -f docker/Dockerfile --target plink2 -t plink2 .
#   docker build -f docker/Dockerfile --target plink1 -t plink1 .
#   docker run --rm -v "$PWD":/data plink2 plink2 --pfile in --freq --out out
#
# Both are built with the same flags and Debian packages as the functional
# test workflow.  The binaries don't assume AVX2/SSE4.2, so the images run on
# any amd64 or arm64 machine.

ARG DEBIAN_RELEASE=bookworm

FROM debian:${DEBIAN_RELEASE}-slim AS build
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential \
      libopenblas-dev \
      liblapack-dev \
      liblapacke-dev \
      zlib1g-dev \
      libzstd-dev \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .

FROM build AS build-plink2
# Link against the system zstd rather than the bundled copy.
RUN make -C 2.0/build_dynamic -j"$(nproc)" STATIC_ZSTD=

FROM build AS build-plink1
RUN make -C 1.9 -j"$(nproc)" ZLIB=-lz BLASFLAGS="-llapack -lopenblas"

FROM debian:${DEBIAN_RELEASE}-slim AS runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      libopenblas0 \
      liblapack3 \
      liblapacke \
      zlib1g \
      libzstd1 \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /data

FROM runtime AS plink1
COPY --from=build-plink1 /src/1.9/plink /usr/local/bin/
CMD ["plink"]

FROM runtime AS plink2
COPY --from=build-plink2 /src/2.0/build_dynamic/plink2 /src/2.0/build_dynamic/pgen_compress /usr/local/bin/
CMD ["plink2"]
