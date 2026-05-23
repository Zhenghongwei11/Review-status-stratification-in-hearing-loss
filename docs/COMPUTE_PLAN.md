# Compute plan

## Hardware
This pipeline can run on a laptop, but it is I/O-bound and may take substantial time on older machines.

## Disk
- The ClinVar `variant_summary.txt.gz` download is large and stored under `data/raw/clinvar/`.
- The derived outputs are stored under `results/` and `plots/` and are comparatively small.

## Runtime
- `scripts/reproduce_one_click.sh` runs the full pipeline.
- For a fast sanity check of parsing and table writes, run:
  ```bash
  bash scripts/hearing_loss_clinvar/smoke_test.sh
  ```

## Determinism
The pipeline is deterministic given a fixed ClinVar snapshot (download URL + file checksum recorded at download time).
