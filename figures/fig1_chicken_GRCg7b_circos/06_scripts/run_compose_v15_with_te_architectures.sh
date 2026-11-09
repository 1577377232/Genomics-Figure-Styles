#!/usr/bin/env bash
# Compose finalized V15 SVG with original TE architecture SVGs; do not redraw V15.
set -euo pipefail
COMPOSE_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
COMPOSE_PROJECT_DIR="$(cd -- "${COMPOSE_SCRIPT_DIR}/.." && pwd)"
COMPOSE_RSCRIPT="${CHICKEN_CIRCOS_RSCRIPT:-Rscript}"
exec "${COMPOSE_RSCRIPT}" --vanilla "${COMPOSE_SCRIPT_DIR}/compose_v15_with_te_architectures.R" \
  --project-dir "${COMPOSE_PROJECT_DIR}" "$@"
