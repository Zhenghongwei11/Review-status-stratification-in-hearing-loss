#!/usr/bin/env python3
"""Describe disease-ID matches among records added by the HPO selection rule.

Counts are distinct VariationIDs per matched disease identifier. A record may
match several identifiers, so rows are not additive and do not assign a
condition-specific classification to the variant.
"""

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

from build_evidence_map import sha256_path, write_tsv
from analysis_inputs import ROOT, RAW_SHA
from run_analysis import canonical_ids


def main():
    source = ROOT / "data/raw/clinvar/variant_summary.txt.gz"
    if sha256_path(source) != RAW_SHA:
        raise RuntimeError("Pinned ClinVar checksum mismatch")
    result_root = ROOT / "results"
    expanded = set()
    with (result_root / "expanded/variant_aggregate.tsv").open(newline="") as handle:
        expanded.update(row["variation_id"] for row in csv.DictReader(handle, delimiter="\t"))
    added = set()
    with (result_root / "disease_supported/variant_aggregate.tsv").open(newline="") as handle:
        added.update(row["variation_id"] for row in csv.DictReader(handle, delimiter="\t") if row["variation_id"] not in expanded)
    names = defaultdict(set)
    with (ROOT / "data/references/hearing_diseases.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for disease_id in canonical_ids(row["disease_id"]):
                names[disease_id].add(row["disease_name"])
    seen = set()
    counts = Counter()
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            vid = row["VariationID"].strip()
            if vid not in added:
                continue
            for disease_id in canonical_ids(row["PhenotypeIDS"]) & names.keys():
                pair = (vid, disease_id)
                if pair in seen:
                    continue
                seen.add(pair)
                counts[disease_id] += 1
    output = ROOT / "results/hpo_disease_contributions.tsv"
    rows = [{"disease_id": did, "disease_name": "; ".join(sorted(names[did])),
             "new_variation_ids_with_matching_disease_id": n,
             "unit": "distinct VariationIDs; disease-ID rows may overlap"}
            for did, n in sorted(counts.items(), key=lambda item: (-item[1], item[0]))]
    write_tsv(output, rows)
    print(f"new_variation_ids={len(added)}; matched_disease_ids={len(counts)}; output={output}")


if __name__ == "__main__":
    main()
