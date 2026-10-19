#!/usr/bin/env python3
"""Export a traceable disease-context display from existing analysis results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = ROOT / "results"
DISEASES = [
    ("ORPHA:154", "Familial isolated dilated cardiomyopathy", 44036),
    ("ORPHA:244", "Primary ciliary dyskinesia", 28729),
    ("ORPHA:84", "Fanconi anemia", 27050),
]
TABLE_COLUMNS = [
    "key", "row_label", "n_variation_ids", "selection_context",
    "classification_label", "review_stars", "observation_type", "source_path",
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def read_tsv(path: Path, key: str, keep=None):
    """Return unique keys, optionally retained rows and retained line numbers."""
    keys, rows, lines = set(), {}, {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or key not in reader.fieldnames:
            raise ValueError(f"Missing {key!r} column in {path}")
        for line, row in enumerate(reader, start=2):
            value = row[key]
            if value in keys:
                raise ValueError(f"Duplicate {key} {value} in {path}:{line}")
            keys.add(value)
            if keep is None or keep(row):
                rows[value], lines[value] = row, line
    return keys, rows, lines


def tokens(value: str, separators: str) -> set[str]:
    for separator in separators[1:]:
        value = value.replace(separator, separators[0])
    return {item.strip() for item in value.split(separators[0]) if item.strip()}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def generate_disease_context_display(results: Path = DEFAULT_RESULTS, output: Path | None = None):
    """Generate the traceable five-row display and return its rows and paths."""
    results = results.resolve()
    output = (output or results).resolve()
    files = {
        "disease_supported": results / "disease_supported/variant_aggregate.tsv",
        "primary": results / "primary/variant_aggregate.tsv",
        "expanded": results / "expanded/variant_aggregate.tsv",
        "coverage_summary": results / "coverage_summary.tsv",
        "hpo_disease_contributions": results / "hpo_disease_contributions.tsv",
        "run_metadata": results / "run_metadata.json",
    }
    for path in files.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    disease_ids, disease_selected, disease_lines = read_tsv(
        files["disease_supported"], "variation_id",
        lambda row: row["variation_id"] == "2353"
        or row["review_status_best"].strip().casefold() == "practice guideline",
    )
    primary_ids, _, _ = read_tsv(files["primary"], "variation_id", lambda _: False)
    expanded_ids, _, _ = read_tsv(files["expanded"], "variation_id", lambda _: False)
    _, coverage, coverage_lines = read_tsv(files["coverage_summary"], "definition")
    _, contributions, contribution_lines = read_tsv(files["hpo_disease_contributions"], "disease_id")
    run_metadata = json.loads(files["run_metadata"].read_text(encoding="utf-8"))

    for name, members in (("disease_supported", disease_ids), ("primary", primary_ids), ("expanded", expanded_ids)):
        require(name in coverage, f"{name} is absent from coverage_summary.tsv")
        require(len(members) == int(coverage[name]["n_variation_ids"]), f"{name} count disagrees with coverage_summary.tsv")
    primary_n = int(coverage["primary"]["n_variation_ids"])
    expanded_n = int(coverage["expanded"]["n_variation_ids"])
    disease_n = int(coverage["disease_supported"]["n_variation_ids"])
    new_ids = disease_ids - expanded_ids
    new_n = len(new_ids)
    require(new_n == disease_n - expanded_n, "Disease-supported minus expanded does not match coverage counts")

    practice = {
        key: row for key, row in disease_selected.items()
        if row["review_status_best"].strip().casefold() == "practice guideline"
    }
    practice_n = int(coverage["disease_supported"]["n_practice_guideline"])
    require(len(practice) == practice_n == 23, "Expected all 23 practice-guideline VariationIDs")
    require(all(r["primary_gene_symbol"] == r["gene_symbols"] == "CFTR" for r in practice.values()), "Practice-guideline group is not all CFTR")
    practice_classifications = {
        classification
        for row in practice.values()
        for classification in row["source_classifications"].split("|")
    }
    practice_stars = {int(row["review_tier_proxy_max"]) for row in practice.values()}
    require(len(practice_classifications) == len(practice_stars) == 1,
            "Practice-guideline records do not share one classification and review tier")
    ush_id = "2353"
    require(ush_id in disease_selected, "USH2A VariationID 2353 is missing")
    ush = disease_selected[ush_id]
    require(ush["primary_gene_symbol"] == ush["gene_symbols"] == "USH2A", "VariationID 2353 is not uniquely annotated to USH2A")
    selected = {ush_id: ush, **practice}

    for variation_id, row in selected.items():
        stars = "3" if variation_id == ush_id else "4"
        require(row["source_classifications"].split("|") == ["Pathogenic"], f"Unexpected classification for {variation_id}")
        require(row["source_classification_review_pairs"].split("|") == [f"Pathogenic@@{stars}"], f"Unexpected joint classification/review pair for {variation_id}")
        require(row["review_tier_proxy_max"] == stars, f"Unexpected review tier for {variation_id}")
        require({"HPO_disease_id", "HPO_disease_name_exact"}.issubset(row["inclusion_reasons"].split("|")), f"Unexpected HPO route for {variation_id}")
        require(variation_id not in primary_ids and variation_id not in expanded_ids, f"{variation_id} occurs in a text-selected set")

    ush_names = tokens(ush["condition_names"], "|;")
    ush_ids = tokens(ush["condition_ids"], "|;,")
    require(any("Usher syndrome" in name for name in ush_names), "VariationID 2353 lacks an Usher syndrome annotation")
    require({"OMIM:276901", "Orphanet:886"}.issubset(ush_ids), "VariationID 2353 lacks expected Usher identifiers")
    for variation_id, row in practice.items():
        require("Cystic fibrosis" in tokens(row["condition_names"], "|;"), f"CFTR {variation_id} lacks Cystic fibrosis annotation")
        require("Orphanet:586" in tokens(row["condition_ids"], "|;,"), f"CFTR {variation_id} lacks Orphanet:586 annotation")

    disease_counts = []
    for disease_id, name, expected in DISEASES:
        require(disease_id in contributions, f"Missing {disease_id} contribution")
        row = contributions[disease_id]
        count = int(row["new_variation_ids_with_matching_disease_id"])
        require(row["disease_name"] == name and count == expected, f"Unexpected contribution for {disease_id}")
        require("distinct VariationIDs" in row["unit"] and "may overlap" in row["unit"], f"Missing unit/nonadditivity note for {disease_id}")
        disease_counts.append((disease_id, name, count, row))

    aggregate_path = relpath(files["disease_supported"])
    rows = [
        {
            "key": "ush2a_2353", "row_label": "USH2A VariationID 2353; Usher syndrome among its annotations",
            "n_variation_ids": len({ush_id} & disease_ids), "selection_context": "HPO-linked disease_supported set; absent from primary and expanded text sets",
            "classification_label": ush["source_classifications"], "review_stars": int(ush["review_tier_proxy_max"]),
            "observation_type": "single_variation_example",
            "source_path": f"{aggregate_path}:{disease_lines[ush_id]}",
        },
        {
            "key": "cftr_practice_guideline_all", "row_label": "CFTR records; cystic fibrosis among every record’s annotations",
            "n_variation_ids": len(practice), "selection_context": "All disease_supported practice-guideline records; absent from primary and expanded text sets",
            "classification_label": next(iter(practice_classifications)), "review_stars": next(iter(practice_stars)),
            "observation_type": "complete_review_status_group",
            "source_path": f"{aggregate_path} (selected source lines in provenance JSON)",
        },
    ]
    for rank, (disease_id, name, count, _) in enumerate(disease_counts, start=1):
        rows.append({
            "key": f"hpo_disease_{rank}", "row_label": f"{name} ({disease_id})", "n_variation_ids": count,
            "selection_context": f"New distinct VariationIDs beyond expanded text with matching disease ID; denominator {new_n}",
            "classification_label": "", "review_stars": "",
            "observation_type": "nonadditive_disease_id_contribution",
            "source_path": f"{relpath(files['hpo_disease_contributions'])}:{contribution_lines[disease_id]}",
        })

    output.mkdir(parents=True, exist_ok=True)
    tsv_path = output / "disease_context_display.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=TABLE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    row_counts = {
        "disease_supported": len(disease_ids), "primary": len(primary_ids), "expanded": len(expanded_ids),
        "coverage_summary": len(coverage), "hpo_disease_contributions": len(contributions),
    }
    input_provenance = {
        name: {"path": relpath(path), "sha256": digest(path), "bytes": path.stat().st_size, "data_rows": row_counts.get(name)}
        for name, path in files.items()
    }
    record_provenance = []
    for variation_id in [ush_id, *sorted(practice, key=int)]:
        row = selected[variation_id]
        names, ids = tokens(row["condition_names"], "|;"), tokens(row["condition_ids"], "|;,")
        record_provenance.append({
            "variation_id": int(variation_id), "source": f"{aggregate_path}:{disease_lines[variation_id]}",
            "gene_symbols": row["gene_symbols"], "classification_review_pair": row["source_classification_review_pairs"],
            "review_status_best": row["review_status_best"], "inclusion_reasons": row["inclusion_reasons"],
            "matched_condition_names": sorted(n for n in names if "Usher syndrome" in n or n == "Cystic fibrosis"),
            "matched_condition_ids": sorted(ids.intersection({"OMIM:276901", "Orphanet:886", "Orphanet:586"})),
            "absent_from_primary": variation_id not in primary_ids, "absent_from_expanded": variation_id not in expanded_ids,
        })
    disease_provenance = [
        {"disease_id": did, "disease_name": name, "n_variation_ids": count, "unit": row["unit"],
         "source": f"{relpath(files['hpo_disease_contributions'])}:{contribution_lines[did]}",
         "classification_label": None, "review_stars": None}
        for did, name, count, row in disease_counts
    ]
    provenance = {
        "schema_version": 1,
        "outputs": {"display_tsv": relpath(tsv_path), "display_tsv_sha256": digest(tsv_path)},
        "inputs": input_provenance,
        "fixed_record_snapshot": {"clinvar_source_sha256": run_metadata.get("source_sha256"),
                                   "clinvar_source_rows": run_metadata.get("total_rows"),
                                   "analysis_started_utc": run_metadata.get("started_utc")},
        "population_counts": {"primary": primary_n, "expanded": expanded_n, "disease_supported": disease_n,
                              "added_beyond_expanded": new_n, "practice_guideline": practice_n,
                              "coverage_source": relpath(files["coverage_summary"]),
                              "coverage_lines": {k: coverage_lines[k] for k in ("primary", "expanded", "disease_supported")}},
        "selected_variation_ids": {"ush2a_example": [2353], "cftr_practice_guideline_group": [int(i) for i in sorted(practice, key=int)]},
        "selected_record_provenance": record_provenance,
        "disease_contribution_provenance": disease_provenance,
        "interpretation_limits": [
            "Classification and review stars are variation-level aggregates, not condition-specific classifications.",
            "Disease-ID contributions may overlap; they are not mutually exclusive and do not establish hearing-loss causality.",
            "Disease-contribution classification fields are blank because condition-specific summaries are not provided.",
        ],
        "validation": {"joint_pairs": True, "24_records_absent_from_primary_and_expanded": True,
                       "23_cftr_records_include_cystic_fibrosis_name_and_orphanet_586": True,
                       "ush2a_2353_has_usher_annotation": True,
                       "new_count_matches_coverage_difference": True, "disease_counts_marked_nonadditive": True},
    }
    json_path = output / "disease_context_display_provenance.json"
    json_path.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {tsv_path}")
    print(f"Wrote {json_path}")
    print(f"Rows: {len(rows)}; selected examples: 1 USH2A + {len(practice)} CFTR; added-beyond-expanded denominator: {new_n}")
    return rows, tsv_path, json_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output-dir", type=Path, help="Defaults to --results-dir")
    args = parser.parse_args()
    generate_disease_context_display(args.results_dir, args.output_dir)


if __name__ == "__main__":
    main()
