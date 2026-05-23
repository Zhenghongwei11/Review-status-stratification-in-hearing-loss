#!/usr/bin/env bash
set -euo pipefail

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="logs/hearing_loss_clinvar/sensitivity/${RUN_ID}"
mkdir -p "${LOG_DIR}"

OUT_A="results/hearing_loss_clinvar_hl_only"
OUT_B="results/hearing_loss_clinvar_expanded_synonyms"

python3 scripts/hearing_loss_clinvar/build_evidence_map.py \
  --variant-summary data/raw/clinvar/variant_summary.txt.gz \
  --outdir "${OUT_A}" \
  --phenotype-mode hearing_loss_only \
  2>&1 | tee "${LOG_DIR}/01_run_hl_only.log"

python3 scripts/hearing_loss_clinvar/build_evidence_map.py \
  --variant-summary data/raw/clinvar/variant_summary.txt.gz \
  --outdir "${OUT_B}" \
  --phenotype-mode expanded_synonyms \
  2>&1 | tee "${LOG_DIR}/02_run_expanded_synonyms.log"

python3 scripts/hearing_loss_clinvar/compare_sensitivity_runs.py \
  --run-a "${OUT_A}" \
  --run-b "${OUT_B}" \
  --label-a hearing_loss_only \
  --label-b expanded_synonyms \
  --top-n 50 \
  2>&1 | tee "${LOG_DIR}/03_compare.log"

echo "[done] results/hearing_loss_clinvar_sensitivity/summary.tsv" | tee "${LOG_DIR}/done.txt"

