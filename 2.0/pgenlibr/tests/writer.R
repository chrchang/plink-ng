# Round-trip tests for the .pgen writer: write with NewPgenWriter() and the
# Append*() functions, then read back with the existing reader.
library(pgenlibr)

set.seed(194)
tmpdir <- tempfile("pgenlibr_writer_")
dir.create(tmpdir)

expect_error <- function(expr) {
  result <- tryCatch({
    expr
    NULL
  }, error = function(e) e)
  stopifnot(inherits(result, "error"))
}

write_pvar <- function(path, alt_alleles) {
  n <- length(alt_alleles)
  lines <- c("#CHROM\tPOS\tID\tREF\tALT",
             paste("1", seq_len(n), paste0("v", seq_len(n)), "A", alt_alleles,
                   sep = "\t"))
  writeLines(lines, path)
}

# 37 samples: not a multiple of the 32-sample word size.
sample_ct <- 37L

# 1. Unphased biallelic hardcalls, from both integer and numeric vectors.
variant_ct <- 60L
geno <- matrix(sample(c(0L, 1L, 2L, NA), sample_ct * variant_ct,
                      replace = TRUE, prob = c(0.5, 0.25, 0.2, 0.05)),
               nrow = sample_ct)
geno[, 2] <- NA_integer_   # all missing
geno[, 3] <- 0L            # monomorphic
geno[, 4] <- 0L            # sparse
geno[c(5, 20), 4] <- c(1L, 2L)
geno[, 5] <- 2L            # sparse relative to hom-alt
geno[7, 5] <- NA_integer_
pgen_path <- file.path(tmpdir, "hardcalls.pgen")
pw <- NewPgenWriter(pgen_path, sample_ct, variant_ct, nonref_flags = FALSE)
for (vidx in seq_len(variant_ct)) {
  if (vidx %% 2 == 0) {
    AppendBiallelic(pw, geno[, vidx])
  } else {
    AppendBiallelic(pw, as.numeric(geno[, vidx]))
  }
}
stopifnot(GetWrittenVariantCt(pw) == variant_ct)
expect_error(AppendBiallelic(pw, geno[, 1]))  # too many variants
ClosePgenWriter(pw)

pgen <- NewPgen(pgen_path)
stopifnot(GetRawSampleCt(pgen) == sample_ct,
          GetVariantCt(pgen) == variant_ct)
stopifnot(identical(ReadIntList(pgen, seq_len(variant_ct)), geno))
ibuf <- IntBuf(pgen)
ReadHardcalls(pgen, ibuf, 4)
stopifnot(identical(ibuf, geno[, 4]))
ClosePgen(pgen)

# -9 is accepted as a missing code.
pgen_path <- file.path(tmpdir, "minus9.pgen")
pw <- NewPgenWriter(pgen_path, 3L, 1L)
AppendBiallelic(pw, c(0L, -9L, 2L))
ClosePgenWriter(pw)
pgen <- NewPgen(pgen_path)
stopifnot(identical(ReadIntList(pgen, 1L)[, 1], c(0L, NA, 2L)))
ClosePgen(pgen)

# 2. Dosages (stored with 1/16384 precision).
variant_ct <- 20L
dosages <- matrix(runif(sample_ct * variant_ct, 0, 2), nrow = sample_ct)
dosages[sample(length(dosages), 30)] <- NA
dosages[, 2] <- sample(c(0, 1, 2), sample_ct, replace = TRUE)  # integral
dosages[1:3, 3] <- c(0, 2, NA)
pgen_path <- file.path(tmpdir, "dosages.pgen")
pw <- NewPgenWriter(pgen_path, sample_ct, variant_ct, dosage_present = TRUE)
for (vidx in seq_len(variant_ct)) {
  AppendDosages(pw, dosages[, vidx])
}
ClosePgenWriter(pw)

pgen <- NewPgen(pgen_path)
buf <- Buf(pgen)
ibuf <- IntBuf(pgen)
for (vidx in seq_len(variant_ct)) {
  expected <- dosages[, vidx]
  Read(pgen, buf, vidx)
  stopifnot(identical(is.na(buf), is.na(expected)),
            all(abs(buf - expected) <= 0.5 / 16384, na.rm = TRUE))
  # Hardcalls are saved when the dosage is within 0.1 of an integer.
  ReadHardcalls(pgen, ibuf, vidx)
  dist <- abs(expected - round(expected))
  near <- !is.na(expected) & (dist < 0.099)
  far <- is.na(expected) | (dist > 0.101)
  stopifnot(identical(ibuf[near], as.integer(round(expected[near]))),
            all(is.na(ibuf[far])))
}
stopifnot(identical(ReadList(pgen, seq_len(variant_ct)) == 0,
                    round(dosages * 16384) == 0))
ClosePgen(pgen)

# Hardcalls and dosages can be mixed in one file.
pgen_path <- file.path(tmpdir, "mixed.pgen")
pw <- NewPgenWriter(pgen_path, 4L, 2L, dosage_present = TRUE)
AppendBiallelic(pw, c(0L, 1L, 2L, NA))
AppendDosages(pw, c(0.25, 1.5, NA, 2))
ClosePgenWriter(pw)
pgen <- NewPgen(pgen_path)
stopifnot(identical(ReadList(pgen, 1:2),
                    matrix(c(0, 1, 2, NA, 0.25, 1.5, NA, 2), nrow = 4)))
ClosePgen(pgen)

# 3. Phased hardcalls via AppendAlleles().
variant_ct <- 30L
pgen_path <- file.path(tmpdir, "phased.pgen")
pw <- NewPgenWriter(pgen_path, sample_ct, variant_ct,
                    hardcall_phase_present = TRUE)
acbufs <- list()
phasepresents <- list()
for (vidx in seq_len(variant_ct)) {
  acbuf <- matrix(sample(0:1, 2 * sample_ct, replace = TRUE), nrow = 2)
  missing <- sample(sample_ct, 2)
  acbuf[, missing] <- NA_integer_
  phasepresent <- sample(c(TRUE, FALSE), sample_ct, replace = TRUE)
  if (vidx == 1) {
    AppendAlleles(pw, acbuf)
    phasepresent[] <- FALSE
  } else if (vidx == 2) {
    AppendAlleles(pw, acbuf, all_phased = TRUE)
    phasepresent[] <- TRUE
  } else if (vidx == 3) {
    # numeric input, -9 missing code
    acbuf_dbl <- acbuf
    storage.mode(acbuf_dbl) <- "double"
    acbuf_dbl[is.na(acbuf_dbl)] <- -9
    AppendAlleles(pw, acbuf_dbl, phasepresent = phasepresent)
  } else {
    AppendAlleles(pw, acbuf, phasepresent = phasepresent)
  }
  acbufs[[vidx]] <- acbuf
  phasepresents[[vidx]] <- phasepresent
}
ClosePgenWriter(pw)

pgen <- NewPgen(pgen_path)
stopifnot(HardcallPhasePresent(pgen))
acbuf_out <- IntAlleleCodeBuf(pgen)
phasepresent_out <- BoolBuf(pgen)
for (vidx in seq_len(variant_ct)) {
  acbuf <- acbufs[[vidx]]
  ReadAlleles(pgen, acbuf_out, vidx, phasepresent_out)
  het <- !is.na(acbuf[1, ]) & (acbuf[1, ] != acbuf[2, ])
  phased_het <- het & phasepresents[[vidx]]
  expected <- acbuf
  # unphased hets are returned in 0/1 order
  expected[, het & !phased_het] <- c(0L, 1L)
  stopifnot(identical(acbuf_out, expected))
  # homozygous calls always have known phase; the reader reports NA for
  # missing calls
  expected_phasepresent <- !het | phased_het
  expected_phasepresent[is.na(acbuf[1, ])] <- NA
  stopifnot(identical(phasepresent_out, expected_phasepresent))
}
ClosePgen(pgen)

# 4. Multiallelic hardcalls, read back through a hand-written .pvar.
variant_ct <- 3L
pgen_path <- file.path(tmpdir, "multiallelic.pgen")
pvar_path <- file.path(tmpdir, "multiallelic.pvar")
write_pvar(pvar_path, c("C", "C,G", "C,G,T"))
acbuf_list <- list(matrix(c(0L, 1L, 1L, 1L, 0L, 0L), nrow = 2),
                   matrix(c(0L, 2L, 1L, 2L, NA, NA), nrow = 2),
                   matrix(c(0L, 0L, 1L, 1L, 0L, 1L), nrow = 2))
pw <- NewPgenWriter(pgen_path, 3L, variant_ct, allele_ct_limit = 4L)
AppendAlleles(pw, acbuf_list[[1]])
AppendAlleles(pw, acbuf_list[[2]])
# The last ALT allele is unobserved, so allele_ct must be given.
AppendAlleles(pw, acbuf_list[[3]], allele_ct = 4L)
ClosePgenWriter(pw)
pvar <- NewPvar(pvar_path)
pgen <- NewPgen(pgen_path, pvar = pvar)
stopifnot(GetAlleleCt(pgen, 2) == 3, GetAlleleCt(pgen, 3) == 4)
ibuf <- IntBuf(pgen)
for (vidx in seq_len(variant_ct)) {
  acbuf <- acbuf_list[[vidx]]
  for (allele_num in seq_len(GetAlleleCt(pgen, vidx))) {
    ReadHardcalls(pgen, ibuf, vidx, allele_num)
    expected <- as.integer(colSums(acbuf == allele_num - 1L))
    stopifnot(identical(ibuf, expected))
  }
}
ClosePgen(pgen)
ClosePvar(pvar)

# 5. Explicit per-variant nonref flags produce a readable file.
pgen_path <- file.path(tmpdir, "nonref.pgen")
pw <- NewPgenWriter(pgen_path, 2L, 3L, nonref_flags = c(TRUE, FALSE, TRUE))
for (vidx in 1:3) {
  AppendBiallelic(pw, c(vidx - 1L, 0L))
}
ClosePgenWriter(pw)
pgen <- NewPgen(pgen_path)
stopifnot(identical(ReadIntList(pgen, 1:3),
                    matrix(c(0:2, 0L, 0L, 0L), nrow = 2, byrow = TRUE)))
ClosePgen(pgen)

# 6. Argument checking.
pgen_path <- file.path(tmpdir, "errors.pgen")
expect_error(NewPgenWriter(pgen_path, 0L, 1L))
expect_error(NewPgenWriter(pgen_path, 2L, 0L))
expect_error(NewPgenWriter(pgen_path, 2L, 3L, nonref_flags = c(TRUE, FALSE)))
expect_error(NewPgenWriter(file.path(tmpdir, "no_such_dir", "x.pgen"), 2L, 1L))
pw <- NewPgenWriter(pgen_path, 3L, 2L)
expect_error(AppendBiallelic(pw, c(0L, 1L)))            # wrong length
expect_error(AppendBiallelic(pw, c(0L, 1L, 3L)))        # invalid value
expect_error(AppendBiallelic(pw, c(0, 0.5, 1)))         # invalid value
expect_error(AppendDosages(pw, c(0, 1, 2)))             # dosage_present FALSE
expect_error(AppendAlleles(pw, matrix(0L, 2, 3), all_phased = TRUE))
expect_error(AppendAlleles(pw, matrix(c(0L, NA, 0L, 0L, 0L, 0L), 2)))  # lone NA
expect_error(AppendAlleles(pw, matrix(c(0L, 2L, 0L, 0L, 0L, 0L), 2)))  # > limit
AppendBiallelic(pw, c(0L, 1L, 2L))
expect_error(ClosePgenWriter(pw))                       # too few variants
stopifnot(GetWrittenVariantCt(pw) == 1)
AppendBiallelic(pw, c(2L, 1L, 0L))
ClosePgenWriter(pw)
expect_error(AppendBiallelic(pw, c(0L, 1L, 2L)))        # closed
ClosePgenWriter(pw)                                     # no-op
pgen <- NewPgen(pgen_path)
expect_error(AppendBiallelic(pgen, c(0L, 1L, 2L)))      # not a writer
stopifnot(identical(ReadIntList(pgen, 1:2), matrix(c(0:2, 2:0), nrow = 3)))
ClosePgen(pgen)

pw <- NewPgenWriter(file.path(tmpdir, "dosage_errors.pgen"), 2L, 1L,
                    dosage_present = TRUE)
expect_error(AppendDosages(pw, c(0, 2.5)))
expect_error(AppendDosages(pw, c(-0.1, 1)))
AppendDosages(pw, c(NA, 1))
ClosePgenWriter(pw)

unlink(tmpdir, recursive = TRUE)
cat("pgenlibr writer tests passed\n")
