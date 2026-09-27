# R logistf reference values for the two-level Firth check in
# TEST_GLM_MULTINOMIAL (needs the logistf package; not run by the test
# itself, its output, ref_firth_two_levels.txt, is committed).
#
# With two levels, --glm multinomial firth is Firth-penalized logistic
# regression.  For each autosomal variant of mn.vcf.gz, this fits phenotype
# CAT2 (hi versus lo) on the five covariates and the ALT dosage (at plink2's
# 1/16384 resolution, samples with a missing call dropped), with tight
# tolerances, and reports:
# * BETA_hi and SE_hi: the ALT coefficient and logistf's standard error,
#   sqrt(diag(var)), where var is the inverse of the leverage-augmented
#   information matrix;
# * LRT: the penalized likelihood ratio statistic, twice the difference
#   between the penalized log-likelihoods of the full fit and of the fit with
#   the ALT coefficient held at zero (logistf's own restricted fit, whose
#   penalty still involves the ALT column).
#
# Usage: Rscript oracle_logistf.R <output>
suppressMessages(library(logistf))
out_fname <- commandArgs(TRUE)[1]

pheno <- read.table("pheno.txt", header = TRUE, comment.char = "")
covar <- read.table("covar.txt", header = TRUE, comment.char = "")
y <- as.integer(pheno$CAT2 == "hi")
vcf <- readLines(gzfile("mn.vcf.gz"))
vcf <- vcf[!startsWith(vcf, "#")]
ctl <- logistf.control(maxit = 1000, maxstep = 5, lconv = 1e-13, gconv = 1e-13, xconv = 1e-13)
rows <- "#ID\tA1\tOBS_CT\tDF\tLRT\tBETA_hi\tSE_hi"
for (line in vcf) {
  f <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (f[1] == "X") {
    next
  }
  fmt <- strsplit(f[9], ":", fixed = TRUE)[[1]]
  cells <- strsplit(f[10:length(f)], ":", fixed = TRUE)
  gt <- sapply(cells, `[`, 1)
  if ("DS" %in% fmt) {
    g <- round(as.numeric(sapply(cells, `[`, match("DS", fmt))) * 16384) / 16384
  } else {
    g <- sapply(strsplit(gt, "/", fixed = TRUE), function(a) sum(as.integer(a)))
  }
  g[gt %in% c("./.", ".")] <- NA
  d <- data.frame(y = y, covar[, 2:6], g = g)
  d <- d[!is.na(d$g), ]
  full <- logistf(y ~ ., data = d, control = ctl, pl = FALSE)
  restricted <- logistf(y ~ ., data = d, control = ctl, pl = FALSE, modcontrol = logistf.mod.control(terms.fit = 1:6))
  lrt <- 2 * (full$loglik[["full"]] - restricted$loglik[["full"]])
  rows <- c(rows, sprintf("%s\t%s\t%d\t1\t%.10g\t%.10g\t%.10g", f[3], f[5], nrow(d), lrt, full$coefficients[["g"]], sqrt(full$var[7, 7])))
}
writeLines(rows, out_fname)
