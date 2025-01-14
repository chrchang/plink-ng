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
    (\code{ReadIntList}, \code{ReadList}).  When you're done, \code{ClosePgen}
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
  # This is modified from https://yosuketanigawa.com/posts/2020/09/PLINK2 .
  library(pgenlibr)

  # These files are subsetted from downloads available at
  #   https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg .
  # Note that, after downloading the original files, the .pgen file must be
  # decompressed before use; but both pgenlibr and the PLINK 2 program can
  # handle compressed .pvar files.
  pvar_path <- system.file("extdata", "chr21_phase3_start.pvar.zst", package="pgenlibr")
  pgen_path <- system.file("extdata", "chr21_phase3_start.pgen", package="pgenlibr")

  pvar <- pgenlibr::NewPvar(pvar_path)
  pgen <- pgenlibr::NewPgen(pgen_path, pvar=pvar)

  # Check the number of variants and samples.
  pgenlibr::GetVariantCt(pgen)
  pgenlibr::GetRawSampleCt(pgen)

  # Get the chromosome, position, and ID of the first variant.
  GetVariantChrom(pvar, 1)
  GetVariantPos(pvar, 1)
  GetVariantId(pvar, 1)

  # Read the 14th variant.
  buf <- pgenlibr::Buf(pgen)
  pgenlibr::Read(pgen, buf, 14)

  # Get the index of the variant with ID "rs569225703".
  var_id <- pgenlibr::GetVariantsById(pvar, "rs569225703")

  # Get allele count.
  pgenlibr::GetAlleleCt(pvar, var_id)

  # It has three alleles, i.e. two ALT alleles.
  # Read first-ALT-allele dosages for that variant.
  pgenlibr::Read(pgen, buf, var_id)

  # Read second-ALT-allele dosages.
  pgenlibr::Read(pgen, buf, var_id, allele_num=3)

  # Read a matrix with both variants.  Note that, for the multiallelic variant,
  # the dosages of both ALT alleles are summed here.
  geno_mat <- pgenlibr::ReadList(pgen, c(14, var_id))

  pgenlibr::ClosePgen(pgen)
  pgenlibr::ClosePvar(pvar)
}
