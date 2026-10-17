#!/usr/bin/env python3
"""Compare the broad ClinGen gene set with a Moderate-or-higher subset."""

import csv
from pathlib import Path

from build_evidence_map import write_tsv
from analysis_inputs import ROOT


def rows(path: Path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def main():
    reference = list(rows(ROOT / "data/references/hearing_genes.tsv"))
    broad = {r["gene_symbol"] for r in reference if r["gene_symbol"] and r["classification"] in {"Definitive", "Strong", "Moderate", "Limited", "Supportive"}}
    moderate = {r["gene_symbol"] for r in reference if r["gene_symbol"] and r["classification"] in {"Definitive", "Strong", "Moderate"}}
    primary = {r["variation_id"] for r in rows(ROOT / "results/primary/variant_aggregate.tsv")}
    counts = {"broad": {"ids": set(), "expert": 0, "plp": 0},
              "moderate_or_higher": {"ids": set(), "expert": 0, "plp": 0}}
    for row in rows(ROOT / "results/curated_gene/variant_aggregate.tsv"):
        symbols = set(row["gene_symbols"].split(";"))
        for label, genes in [("broad", broad), ("moderate_or_higher", moderate)]:
            if not symbols & genes:
                continue
            item = counts[label]
            item["ids"].add(row["variation_id"])
            item["expert"] += int(row["review_tier_proxy_max"]) == 3
            item["plp"] += int(row["plp_review_stars_ge_2_same_row"])
    output = ROOT / "results/clingen_grade_sensitivity.tsv"
    values = []
    for label, genes in [("broad", broad), ("moderate_or_higher", moderate)]:
        item = counts[label]
        values.append({"evidence_threshold": label, "qualifying_gene_symbols": len(genes),
                       "selected_variation_ids": len(item["ids"]),
                       "variation_ids_shared_with_primary_phrase": len(item["ids"] & primary),
                       "expert_panel_records": item["expert"],
                       "plp_review_stars_ge_2_records": item["plp"]})
    write_tsv(output, values)
    print(output)


if __name__ == "__main__":
    main()
