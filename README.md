# ClinVar evidence mapping for genetic hearing loss (public reproducibility package)

This repository provides a lightweight, reviewer-facing package to reproduce the analysis outputs (tables and figures) for a ClinVar-first evidence map of genetic hearing loss.

## What this repo contains
- `scripts/`: end-to-end pipeline entrypoints
- `results/`: derived tables used to generate the reported figures/tables
- `plots/`: exported figures (`.png` and `.pdf`)
- `docs/`: data manifest, figure provenance, and declared decision rules

## Quick start (one click)
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

bash scripts/reproduce_one_click.sh
```

## What to expect (compute)
- ClinVar `variant_summary.txt.gz` is a large public download; ensure ample disk space under `data/raw/`.
- The pipeline processes the file in a streaming manner, but runtime can be long on laptops.

## Outputs
- Tables: `results/hearing_loss_clinvar/`
- Figures: `plots/hearing_loss_clinvar/`

## License
MIT (see `LICENSE`).
