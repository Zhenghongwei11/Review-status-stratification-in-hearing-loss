#!/usr/bin/env bash
set -euo pipefail

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="logs/hearing_loss_clinvar/smoke_${RUN_ID}"
mkdir -p "${LOG_DIR}"

# Build a small head subset for smoke testing (fast; validates parsing + table writes).
python3 scripts/hearing_loss_clinvar/make_variant_summary_head.py \
  --out data/raw/clinvar/variant_summary.head.tsv.gz \
  --lines 200000 \
  2>&1 | tee "${LOG_DIR}/01_make_head.log"

# Parse only a small number of rows to validate logic and column expectations.
python3 scripts/hearing_loss_clinvar/build_evidence_map.py \
  --variant-summary data/raw/clinvar/variant_summary.head.tsv.gz \
  --outdir results/hearing_loss_clinvar_smoke \
  2>&1 | tee "${LOG_DIR}/02_build_smoke.log"

echo "[done] results/hearing_loss_clinvar_smoke" | tee "${LOG_DIR}/done.txt"
