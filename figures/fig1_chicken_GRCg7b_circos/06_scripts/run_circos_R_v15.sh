#!/usr/bin/env bash
# V15：真实数据；TE 独立纵轴、平滑轮廓；输出 07_figures/R_v15。
set -euo pipefail
CIRCOS_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CIRCOS_PROJECT_DIR="$(cd -- "${CIRCOS_SCRIPT_DIR}/.." && pwd)"
CIRCOS_RSCRIPT="${CHICKEN_CIRCOS_RSCRIPT:-${R_BIN:-Rscript}}"
if ! command -v "${CIRCOS_RSCRIPT}" >/dev/null 2>&1; then
  printf '找不到 Rscript：%s\n' "${CIRCOS_RSCRIPT}" >&2
  exit 127
fi
exec "${CIRCOS_RSCRIPT}" --vanilla "${CIRCOS_SCRIPT_DIR}/plot_chicken_circos_v15.R" \
  --project-dir "${CIRCOS_PROJECT_DIR}" "$@"
