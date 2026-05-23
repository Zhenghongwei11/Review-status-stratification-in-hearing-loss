#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

echo "[step] (1/4) Download ClinVar variant_summary export (resume + sha256)"
python3 scripts/hearing_loss_clinvar/download_clinvar.py \
  --out data/raw/clinvar/variant_summary.txt.gz \
  --write-sha256 data/raw/clinvar/variant_summary.txt.gz.sha256

echo "[step] (2/4) Build evidence-map tables (hearing-loss phenotype filter)"
python3 scripts/hearing_loss_clinvar/build_evidence_map.py \
  --variant-summary data/raw/clinvar/variant_summary.txt.gz \
  --outdir results/hearing_loss_clinvar \
  --phenotype-mode hearing_loss_only

echo "[step] (3/4) Generate utility tables"
python3 scripts/hearing_loss_clinvar/generate_utility_tables.py

echo "[step] (4/4) Render figures"
python3 scripts/hearing_loss_clinvar/plot_figures.py

echo "[done] results/hearing_loss_clinvar + plots/hearing_loss_clinvar"
