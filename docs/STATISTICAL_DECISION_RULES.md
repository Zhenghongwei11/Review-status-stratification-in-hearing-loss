# Declared decision rules (analysis invariants)

This document predeclares the core analysis rules used to derive the tables under `results/` from the ClinVar `variant_summary` export.

## Unit of counting
- Unique variant unit: **ClinVar `VariationID`** (used to de-duplicate across submissions/records within the export).

## Hearing-loss phenotype filter
- Inclusion is based on a text-pattern match in the ClinVar phenotype list field, using the curated patterns implemented in `scripts/hearing_loss_clinvar/phenotype_filters.py`.
- A broadened synonym set is used for the sensitivity run under `results/hearing_loss_clinvar_sensitivity/`.

## Review-status tier proxy (0–4)
ClinVar review status wording is mapped to an approximate 0–4 tier:
- Tier 4: practice guideline
- Tier 3: reviewed by expert panel
- Tier 2: criteria provided, multiple submitters, no conflicts
- Tier 1: criteria provided, single submitter
- Tier 0: no assertion criteria provided (and all other/unknown statuses default to tier 0)

## Variant-level clinical-significance bucket
Row-level clinical significance strings are mapped to one of:
- `pathogenic_or_likely_pathogenic`
- `vus`
- `benign_or_likely_benign`
- `conflicting`
- `other` / `unknown`

When multiple row-level buckets exist for the same `VariationID`, the variant-level bucket is aggregated conservatively using:
`conflicting > pathogenic_or_likely_pathogenic > vus > benign_or_likely_benign > other/unknown`.

## Gene assignment for gene-level summaries
- If multiple gene symbols are present for a `VariationID`, the primary gene is the first symbol in lexicographic order after trimming.

## Rounding
- Proportions in figures are rendered as percentages and are rounded by matplotlib defaults; exact counts are preserved in the `results/` tables.
