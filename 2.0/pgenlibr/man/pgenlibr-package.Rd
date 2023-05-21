\name{pgenlibr-package}
\alias{pgenlibr-package}
\alias{pgenlibr}
\docType{package}
\title{
    PLINK 2 Binary (.pgen) Reader
}
\description{
    A thin wrapper over PLINK 2's core libraries which provides an R interface
    for reading .pgen files.  A minimal .pvar loader is also included.
}
\details{
    \code{NewPvar} and \code{NewPgen} initialize the respective readers.  Then,
    you can either iterate through one variant at a time (\code{Read},
    \code{ReadAlleles}) or perform a multi-variant matrix load
    (\code{ReadIntList}, \read{ReadList}).  When you're done, \code{ClosePgen}
    and \code{ClosePvar} free resources.
}
\author{
    Christopher Chang \email{chrchang@alumni.caltech.edu}
}
\references{
  Chang, C.C. and Chow, C.C. and Tellier, L.C.A.M. and Vattikuti, S. and
  Purcell, S.M. and Lee J.J. (2015) Second-generation PLINK: rising to the
  challenge of larger and richer datasets.  \emph{Gigascience} 4:7.
  \doi{10.1186/s13742-015-0047-8}.
}
\keyword{ package }
\examples{
  \dontrun{
    # This is modified from https://yosuketanigawa.com/posts/2020/09/PLINK2 .
    library(pgenlibr)

    # These files can be downloaded from
    # https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg .
    # The .pgen file must be decompressed before use, but both pgenlibr and
    # the PLINK 2 program can handle compressed .pvar files.
    pfile <- 'chr21_hg38'

    pvar <- pgenlibr::NewPvar(paste0(pfile, '.pvar.zst'))
    pgen <- pgenlibr::NewPgen(paste0(pfile, '.pgen'), pvar=pvar)

    # Check the number of variants and samples.
    pgenlibr::GetVariantCt(pgen)
    pgenlibr::GetRawSampleCt(pgen)

    # Get the ID of the first variant.
    GetVariantId(pvar, 1)

    # Read the first variant.
    buf <- pgenlibr::Buf(pgen)
    pgenlibr::Read(pgen, buf, 1)

    # Get the index of the variant with ID "21:5031194:G:A".
    var_id <- pgenlibr::GetVariantsById(pvar, "21:5031194:G:A")

    # Read that variant.
    pgenlibr::Read(pgen, buf, var_id)

    # Read a matrix with both variants.
    geno_mat <- pgenlibr::ReadList(pgen, c(1, var_id))

    pgenlibr::ClosePgen(pgen)
    pgenlibr::ClosePvar(pvar)
  }
}
