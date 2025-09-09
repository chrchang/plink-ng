# This file is currently incompatible with CRAN submission (causes Windows .dll
# to contain abort() / exit() calls).
files <- Sys.glob(paste0("*", SHLIB_EXT))
files <- c(files, dir(pattern="a$"))
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
