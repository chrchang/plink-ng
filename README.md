# plink-ng

[![Functional tests](https://img.shields.io/github/actions/workflow/status/chrchang/plink-ng/plink2_functional_tests.yml?branch=master&label=tests&logo=github)](https://github.com/chrchang/plink-ng/actions/workflows/plink2_functional_tests.yml)
[![pgenlib tests](https://img.shields.io/github/actions/workflow/status/chrchang/plink-ng/ci.yaml?branch=master&label=pgenlib&logo=github)](https://github.com/chrchang/plink-ng/actions/workflows/ci.yaml)
[![License](https://img.shields.io/badge/license-GPLv3%20%2F%20LGPLv3-blue)](2.0/COPYING)
[![Paper](https://img.shields.io/badge/DOI-10.1186%2Fs13742--015--0047--8-blue)](https://doi.org/10.1186/s13742-015-0047-8)
[![Forum](https://img.shields.io/badge/support-plink2--users-orange)](https://groups.google.com/g/plink2-users)

[![Homebrew plink1](https://img.shields.io/homebrew/v/plink1?label=brew%20plink1&logo=homebrew)](https://formulae.brew.sh/formula/plink1)
[![Homebrew plink-ng](https://img.shields.io/homebrew/v/plink-ng?label=brew%20plink-ng&logo=homebrew)](https://formulae.brew.sh/formula/plink-ng)
[![Bioconda plink](https://img.shields.io/conda/vn/bioconda/plink?label=bioconda%20plink&logo=anaconda)](https://anaconda.org/bioconda/plink)
[![Bioconda plink2](https://img.shields.io/conda/vn/bioconda/plink2?label=bioconda%20plink2&logo=anaconda)](https://anaconda.org/bioconda/plink2)
[![PyPI pgenlib](https://img.shields.io/pypi/v/pgenlib?label=pypi%20pgenlib&logo=python&logoColor=white)](https://pypi.org/project/pgenlib/)
[![CRAN pgenlibr](https://img.shields.io/cran/v/pgenlibr?label=CRAN%20pgenlibr&logo=r)](https://CRAN.R-project.org/package=pgenlibr)

Source code for PLINK 1.9 and PLINK 2.0, the successors to Shaun Purcell's
PLINK 1.07 (https://zzz.bwh.harvard.edu/plink/), a whole-genome association
analysis toolset.

- User documentation and prebuilt binaries: [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/), [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/)
- Technical support forum: https://groups.google.com/g/plink2-users
- Main methods paper: [Second-generation PLINK](https://academic.oup.com/gigascience/article/4/1/s13742-015-0047-8/2707533) (GigaScience, 2015)

## Which version?

**PLINK 1.9** (`1.9/`) can typically be used as a drop-in replacement for
PLINK 1.07 that scales to much larger datasets. It works with the .bed/.bim/.fam
fileset. It's technically still a beta version because there are a few
rarely-used but possibly-worthwhile PLINK 1.07 commands that are still absent,
but active feature development for it ended in 2016.

**PLINK 2.0** (`2.0/`) is designed to handle VCF files and dosage data, and is
under active development. Its native format is the .pgen/.pvar/.psam fileset,
which stores dosages, phase and multiallelic variants; it also reads and writes
PLINK 1 filesets. Most basic features other than non-concatenating merge are
now in place. See [2.0/README.md](2.0/README.md) for more details.

## Installing

Both versions are available from Homebrew (macOS and Linux):

```
brew install plink1    # PLINK 1.9, installs `plink`
brew install plink-ng  # PLINK 2.0, installs `plink2` and `pgen_compress`
```

Bioconda packages both as well (`conda install -c conda-forge -c bioconda plink plink2`), and
BioContainers builds Docker/Singularity images from those packages
(`quay.io/biocontainers/plink2`, `quay.io/biocontainers/plink`).

Homebrew and Bioconda track tagged releases. The documentation pages above have
the latest builds for Linux, macOS and Windows.

## Building from source

PLINK 2.0 needs a C/C++ compiler, zlib and zstd, plus BLAS/LAPACK (Accelerate
on macOS; e.g. OpenBLAS on Linux):

```
cd 2.0/build_dynamic
make -j8
```

This builds `plink2` and `pgen_compress`. The Makefile header lists the build
options (AVX2, MKL/AOCL, static linking, no LAPACK, ...).

PLINK 1.9:

```
cd 1.9
make -j8 ZLIB=-lz                                          # macOS
make -j8 ZLIB=-lz BLASFLAGS="-llapack -lopenblas"         # Linux
```

On Debian/Ubuntu, the dependencies are `build-essential libopenblas-dev
liblapack-dev liblapacke-dev zlib1g-dev libzstd-dev`.

## Libraries

- **C/C++** (`2.0/include/`): two LGPL-licensed libraries, pgenlib (reads and
  writes .pgen files) and plink2_text (fast line reader with gzip/Zstd
  decompression). The .pgen format is specified in
  [pgen_spec/pgen_spec.pdf](pgen_spec/pgen_spec.pdf).
- **Python**: `pip install pgenlib` ([2.0/Python](2.0/Python)).
- **R**: `install.packages("pgenlibr")` ([CRAN](https://CRAN.R-project.org/package=pgenlibr), source in [2.0/pgenlibr](2.0/pgenlibr)).

## Tests

`2.0/Tests` holds a self-contained suite that compares plink2 against plink 1.9
and checks a number of round trips. It needs `plink` (1.9) on the `PATH`:

```
cd 2.0/Tests
./run_tests.sh ../build_dynamic
```

## Repository layout

| Directory | Contents |
| --- | --- |
| `1.9/` | PLINK 1.9 |
| `2.0/` | PLINK 2.0, with `include/` (pgenlib), `Python/`, `pgenlibr/`, `utils/` and `Tests/` |
| `pgen_spec/` | .pgen format specification |

## License

PLINK 1.9 and PLINK 2.0 are GPLv3 (see [1.9/LICENSE](1.9/LICENSE) and
[2.0/COPYING](2.0/COPYING)). The libraries in `2.0/include/` are LGPLv3
([2.0/COPYING.LESSER](2.0/COPYING.LESSER)).
