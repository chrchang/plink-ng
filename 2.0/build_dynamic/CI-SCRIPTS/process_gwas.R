# =====================================================
# R GWAS pipeline mimicking PLINK 2.0 --glm
# Memory-efficient, batch-enabled, Firth regression safe
# =====================================================

# ------------------- Dependencies -------------------
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

install_if_missing("data.table")
install_if_missing("logistf")
install_if_missing("jsonlite")

library(data.table)
library(logistf)

# ------------------- Prepare data -------------------
prepare_data <- function(effects_file, pheno_file,
                         effects_start_col = 1,
                         iid_col = "IID",
                         pheno_col,
                         covar_cols = NULL,
                         remove_hash_FID = TRUE,
                         pheno_tab_delim = TRUE) {
  
  # Load phenotype/covariates
  sep <- ifelse(pheno_tab_delim, "\t", ",")
  Y_all <- fread(pheno_file, sep = sep, data.table = FALSE)
  colnames(Y_all) <- trimws(colnames(Y_all))
  if (remove_hash_FID) colnames(Y_all) <- gsub("^#", "", colnames(Y_all))
  Y_all[[iid_col]] <- as.character(Y_all[[iid_col]])
  
  stopifnot(iid_col %in% colnames(Y_all))
  stopifnot(pheno_col %in% colnames(Y_all))
  if (!is.null(covar_cols)) stopifnot(all(covar_cols %in% colnames(Y_all)))
  
  merged_Y <- Y_all[, c(iid_col, pheno_col, covar_cols), drop=FALSE]
  
  # Load SNP headers only
  X_all <- fread(effects_file, nrows = 0)
  colnames(X_all) <- trimws(colnames(X_all))
  predictor_cols <- colnames(X_all)[effects_start_col:ncol(X_all)]
  
  cat(sprintf("✅ Loaded phenotype/covariates: %d samples\n", nrow(merged_Y)))
  cat(sprintf("✅ Detected %d SNPs in effects file\n", length(predictor_cols)))
  
  return(list(merged_Y = merged_Y, snp_cols = predictor_cols))
}

# ------------------- GWAS regression -------------------
# =====================================================
# R GWAS pipeline mimicking PLINK 2.0 --glm
# Memory-efficient: batch loading SNPs
# Supports linear, logistic, and Firth logistic regression
# =====================================================

# ------------------- Dependencies -------------------
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

install_if_missing("data.table")
install_if_missing("logistf")

library(data.table)
library(logistf)

# ------------------- Prepare data -------------------
prepare_data <- function(effects_file, pheno_file,
                         effects_start_col = 1,
                         iid_col = "IID",
                         pheno_col,
                         covar_cols = NULL,
                         remove_hash_FID = TRUE,
                         pheno_tab_delim = TRUE) {
  
  sep <- ifelse(pheno_tab_delim, "\t", ",")
  Y_all <- fread(pheno_file, sep = sep, data.table = FALSE)
  colnames(Y_all) <- trimws(colnames(Y_all))
  if (remove_hash_FID) colnames(Y_all) <- gsub("^#", "", colnames(Y_all))
  Y_all[[iid_col]] <- as.character(Y_all[[iid_col]])
  
  stopifnot(iid_col %in% colnames(Y_all))
  stopifnot(pheno_col %in% colnames(Y_all))
  if (!is.null(covar_cols)) stopifnot(all(covar_cols %in% colnames(Y_all)))
  
  merged_Y <- Y_all[, c(iid_col, pheno_col, covar_cols), drop=FALSE]
  
  # Load SNP headers only
  X_all <- fread(effects_file, nrows = 0)
  colnames(X_all) <- trimws(colnames(X_all))
  predictor_cols <- colnames(X_all)[effects_start_col:ncol(X_all)]
  
  cat(sprintf("✅ Loaded phenotype/covariates: %d samples\n", nrow(merged_Y)))
  cat(sprintf("✅ Detected %d SNPs in effects file\n", length(predictor_cols)))
  
  return(list(merged_Y = merged_Y, snp_cols = predictor_cols))
}

# ------------------- GWAS regression -------------------
run_gwas <- function(effects_file,
                     pheno_file,
                     effects_start_col = 7,
                     iid_col = "IID",
                     pheno_col,
                     covar_cols = NULL,
                     regression_type = c("linear", "logistic", "firth"),
                     batch_size = 500,
                     output_prefix = "results/test_run") {
  
  regression_type <- match.arg(regression_type)
  
  # Prepare data
  data_prep <- prepare_data(effects_file, pheno_file,
                            effects_start_col, iid_col, pheno_col, covar_cols)
  merged_Y <- data_prep$merged_Y
  snp_cols <- data_prep$snp_cols
  n_snps <- length(snp_cols)
  
  results <- data.frame(SNP = character(),
                        beta = numeric(),
                        se = numeric(),
                        p = numeric(),
                        stringsAsFactors = FALSE)
  
  # Split SNPs into batches
  if (is.null(batch_size)) {
    batches <- list(snp_cols)
  } else {
    batches <- split(snp_cols, ceiling(seq_along(snp_cols)/batch_size))
  }
  
  snp_counter <- 0
  for (batch in batches) {
    # Load SNP batch
    snp_df <- tryCatch({
      fread(effects_file, select = c(iid_col, batch), data.table = FALSE)
    }, error = function(e) {
      warning(sprintf("Skipping batch %s: %s", paste(batch, collapse=", "), e$message))
      return(NULL)
    })
    if (is.null(snp_df)) next
    snp_df[[iid_col]] <- as.character(snp_df[[iid_col]])
    
    # Merge with phenotype/covariates
    merged_batch <- merge(merged_Y, snp_df, by = iid_col)
    
    for (snp in batch) {
      snp_counter <- snp_counter + 1
      cat(sprintf("Processing SNP %d / %d (%.2f%%): %s\n", 
                  snp_counter, n_snps, 100*snp_counter/n_snps, snp))
      flush.console()
      
      y_sub <- merged_batch[[pheno_col]]
      x_snp <- merged_batch[[snp]]
      X_covars_sub <- if (!is.null(covar_cols)) as.data.frame(merged_batch[, covar_cols, drop=FALSE]) else NULL
      
      # Remove missing
      covar_na <- if (!is.null(X_covars_sub)) apply(is.na(X_covars_sub), 1, any) else rep(FALSE, length(y_sub))
      keep <- !is.na(y_sub) & !is.na(x_snp) & !covar_na
      y_sub <- y_sub[keep]
      x_snp <- x_snp[keep]
      if (!is.null(X_covars_sub)) X_covars_sub <- X_covars_sub[keep, , drop=FALSE]
      if (length(y_sub) < 10) next
      
      # Recode phenotype for logistic/firth
      if (regression_type %in% c("logistic", "firth")) {
        unique_vals <- sort(unique(y_sub))
        if (!all(unique_vals %in% c(0, 1))) {
          max_val <- max(unique_vals)
          y_sub <- ifelse(y_sub == max_val, 1, 0)
        }
      }
      
      # Prepare regression df
      df <- data.frame(y = y_sub, SNP = x_snp)
      if (!is.null(X_covars_sub)) df <- cbind(df, X_covars_sub)
      
      beta <- se <- pval <- NA
      # --- Regression with error handling ---
      tryCatch({
        if (regression_type == "linear") {
          fit <- lm(y ~ ., data = df)
          s <- summary(fit)
          if ("SNP" %in% rownames(s$coefficients)) {
            beta <- s$coefficients["SNP", 1]
            se   <- s$coefficients["SNP", 2]
            pval <- s$coefficients["SNP", 4]
          }
        } else if (regression_type == "logistic") {
          fit <- glm(y ~ ., family = binomial(), data = df)
          s <- summary(fit)
          if ("SNP" %in% rownames(s$coefficients)) {
            beta <- s$coefficients["SNP", 1]
            se   <- s$coefficients["SNP", 2]
            pval <- s$coefficients["SNP", 4]
          }
        } else if (regression_type == "firth") {
  formula_str <- paste("y ~ SNP", 
                       if (!is.null(covar_cols)) paste("+", paste(covar_cols, collapse=" + ")) else "")
  fit <- logistf::logistf(as.formula(formula_str), data = df, pl=TRUE)
  beta <- fit$coefficients["SNP"]
  # Estimate SE from Wald-style CI if var is NA
  se <- sqrt(diag(fit$var))["SNP"]
  if (is.na(se)) {
    se <- (fit$ci.upper["SNP"] - fit$ci.lower["SNP"]) / (2 * 1.96)
  }
  pval <- fit$prob["SNP"]
}
      }, error = function(e) {
        cat("\n❌ Regression failed for SNP:", snp, "\n")
        cat("   Error message:", e$message, "\n")
        beta <<- se <<- pval <<- NA
      })
      
      # Print per-SNP results
      # cat(sprintf("   Result: beta=%.6g, se=%.6g, p=%.6g\n", beta, se, pval))
      
      results <- rbind(results, data.frame(SNP = snp, beta = beta, se = se, p = pval))
    }
  }
  
  # --- Save results ---
  out_dir <- dirname(output_prefix)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_csv <- paste0(output_prefix, "_", regression_type, ".csv")
  write.csv(results, out_csv, row.names = FALSE)
  cat("✅ Saved CSV: ", out_csv, "\n")
  cat("✅ All SNPs processed.\n")
  
  return(results)
}

# ------------------- Example usage -------------------

path <- "test_data/"
fname <- "1kgp3_50k_nomiss_Av_nonintdose"
effects_file <- paste0(path, fname, "_recode_varIDs_A.raw")
pheno_file   <- paste0(path, fname, "_combined_phenocov.csv")
cov_ <- c("COV_1")
model <- "firth"
phenotype <- "ybool"

results <- run_gwas(
  effects_file      = effects_file,
  pheno_file        = pheno_file,
  effects_start_col = 7,
  iid_col           = "IID",
  pheno_col         = phenotype,
  covar_cols        = cov_,
  regression_type   = model,
  batch_size        = 20000,
  output_prefix     = paste0("./results/", fname, "_", phenotype, "_", paste(cov_, collapse="_"), "_glm_", model)
)
