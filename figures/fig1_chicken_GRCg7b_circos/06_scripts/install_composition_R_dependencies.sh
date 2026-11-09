#!/usr/bin/env bash
# Standalone composition dependencies: isolated from the V15 plotting library.
set -euo pipefail
COMPOSE_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
COMPOSE_RSCRIPT="${CHICKEN_CIRCOS_RSCRIPT:-Rscript}"
"${COMPOSE_RSCRIPT}" --vanilla - "${COMPOSE_SCRIPT_DIR}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
lib <- file.path(args[1], "composition_R_library")
archives <- file.path(args[1], "composition_dependencies")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
dir.create(archives, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(lib), .libPaths()))
options(timeout = 180)
versions <- c(png = "0.1-9", rsvg = "2.7.0", cpp11 = "0.5.5",
              systemfonts = "1.0.5", svglite = "2.1.3")
for (pkg in names(versions)) {
  version <- versions[[pkg]]
  found <- tryCatch(packageVersion(pkg), error = function(e) NULL)
  if (!is.null(found) && found == numeric_version(version)) next
  archive <- file.path(archives, paste0(pkg, "_", version, ".tar.gz"))
  if (!file.exists(archive)) {
    success <- FALSE
    for (base in c("https://cloud.r-project.org/src/contrib/",
                   paste0("https://cloud.r-project.org/src/contrib/Archive/", pkg, "/"))) {
      success <- tryCatch({ download.file(paste0(base, basename(archive)), archive,
        mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
      if (success) break
    }
    if (!success) stop("Cannot download ", basename(archive))
  }
  install.packages(archive, lib = lib, repos = NULL, type = "source")
  stopifnot(packageVersion(pkg) == numeric_version(version))
}
cat("Composition library:", normalizePath(lib), "\n")
print(sapply(names(versions), function(pkg) as.character(packageVersion(pkg))))
RS
