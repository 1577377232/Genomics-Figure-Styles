#!/usr/bin/env bash
# 固定兼容版本，安装到项目本地 R_library；不写入系统 R 包目录。
set -euo pipefail
CIRCOS_INSTALL_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CIRCOS_INSTALL_R="${CHICKEN_CIRCOS_RSCRIPT:-Rscript}"
"${CIRCOS_INSTALL_R}" --vanilla - "${CIRCOS_INSTALL_DIR}/R_library" <<'RS'
options(timeout = 240)
lib <- commandArgs(trailingOnly = TRUE)[1]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(lib), .libPaths()))
# GlobalOptions 0.1.4 使用旧版 R 没有的 utils::findMatches；0.1.2 兼容本机 R 4.1.2。
versions <- c(GlobalOptions = "0.1.2", shape = "1.4.6.1", colorspace = "2.1-4", circlize = "0.4.18")
for (pkg in names(versions)) {
  version <- versions[[pkg]]
  found <- tryCatch(utils::packageVersion(pkg), error = function(e) NULL)
  if (!is.null(found) && found == numeric_version(version)) next
  filename <- paste0(pkg, "_", version, ".tar.gz")
  urls <- c(paste0("https://cloud.r-project.org/src/contrib/", filename),
            paste0("https://cloud.r-project.org/src/contrib/Archive/", pkg, "/", filename))
  archive <- tempfile(fileext = ".tar.gz")
  success <- FALSE
  for (url in urls) {
    success <- tryCatch({ utils::download.file(url, archive, mode = "wb", quiet = TRUE); TRUE },
                        error = function(e) FALSE)
    if (success) break
  }
  if (!success) stop("无法下载兼容版本：", filename)
  utils::install.packages(archive, lib = lib, repos = NULL, type = "source")
  unlink(archive)
}
stopifnot(requireNamespace("circlize", quietly = TRUE))
cat("circlize:", as.character(utils::packageVersion("circlize")), "\n")
cat("项目 R library:", normalizePath(lib), "\n")
RS
