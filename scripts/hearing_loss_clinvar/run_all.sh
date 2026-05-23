#!/usr/bin/env bash
set -euo pipefail

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="logs/hearing_loss_clinvar/${RUN_ID}"
mkdir -p "${LOG_DIR}"

echo "[run_id] ${RUN_ID}" | tee "${LOG_DIR}/run_id.txt"

if [[ -f data/raw/clinvar/variant_summary.txt.gz.partial && ! -f data/raw/clinvar/variant_summary.txt.gz ]]; then
  echo "[error] Detected in-progress download: data/raw/clinvar/variant_summary.txt.gz.partial" | tee "${LOG_DIR}/01_download_clinvar.log"
  echo "[error] Wait for it to finish, or delete the .partial and rerun." | tee -a "${LOG_DIR}/01_download_clinvar.log"
  exit 2
fi

if [[ -f data/raw/clinvar/variant_summary.txt.gz && -s data/raw/clinvar/variant_summary.txt.gz && -f data/raw/clinvar/variant_summary.txt.gz.sha256 ]]; then
  echo "[skip] ClinVar already downloaded with sha256 present." | tee "${LOG_DIR}/01_download_clinvar.log"
else
  python3 scripts/hearing_loss_clinvar/download_clinvar.py \
    --out data/raw/clinvar/variant_summary.txt.gz \
    --write-sha256 data/raw/clinvar/variant_summary.txt.gz.sha256 \
    2>&1 | tee "${LOG_DIR}/01_download_clinvar.log"
fi

python3 scripts/hearing_loss_clinvar/build_evidence_map.py \
  --variant-summary data/raw/clinvar/variant_summary.txt.gz \
  --outdir results/hearing_loss_clinvar \
  --phenotype-mode hearing_loss_only \
  2>&1 | tee "${LOG_DIR}/02_build_evidence_map.log"

python3 scripts/hearing_loss_clinvar/generate_utility_tables.py 2>&1 | tee "${LOG_DIR}/03_generate_utility_tables.log"

python3 scripts/hearing_loss_clinvar/plot_figures.py 2>&1 | tee "${LOG_DIR}/04_plot_figures.log"

echo "[done] results/hearing_loss_clinvar" | tee "${LOG_DIR}/done.txt"
