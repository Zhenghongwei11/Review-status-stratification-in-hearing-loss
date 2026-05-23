#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_qc(path: Path) -> dict[str, str]:
    m = {}
    with path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            k = (row.get("metric") or "").strip()
            v = (row.get("value") or "").strip()
            if k:
                m[k] = v
    return m


def read_gene_set(path: Path, *, top_n: int) -> set[str]:
    rows = []
    with path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            g = (row.get("gene_symbol") or "").strip()
            if not g:
                continue
            def to_int(x: str) -> int:
                try:
                    return int(x)
                except Exception:
                    return 0
            rows.append(
                (
                    g,
                    to_int((row.get("n_plp") or "").strip()),
                    to_int((row.get("n_vus") or "").strip()),
                    to_int((row.get("n_conflicting") or "").strip()),
                )
            )
    # Sort consistent with "top genes by P/LP count, then VUS, then conflicting".
    rows.sort(key=lambda t: (-t[1], -t[2], -t[3], t[0]))
    return set([t[0] for t in rows[:top_n]])


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def main() -> int:
    ap = argparse.ArgumentParser(description="Compare hearing-loss phenotype sensitivity runs (summary-level).")
    ap.add_argument("--run-a", required=True, help="Output dir A (contains qc_summary.tsv and gene_evidence_map.tsv)")
    ap.add_argument("--run-b", required=True, help="Output dir B (contains qc_summary.tsv and gene_evidence_map.tsv)")
    ap.add_argument("--label-a", default="hearing_loss_only")
    ap.add_argument("--label-b", default="expanded_synonyms")
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--out", default="results/hearing_loss_clinvar_sensitivity/summary.tsv")
    args = ap.parse_args()

    a = Path(args.run_a)
    b = Path(args.run_b)
    qc_a = read_qc(a / "qc_summary.tsv")
    qc_b = read_qc(b / "qc_summary.tsv")

    genes_a = read_gene_set(a / "gene_evidence_map.tsv", top_n=args.top_n)
    genes_b = read_gene_set(b / "gene_evidence_map.tsv", top_n=args.top_n)

    rows = [
        {
            "run_label": args.label_a,
            "hearing_loss_rows_record_level": qc_a.get("hearing_loss_rows_record_level", ""),
            "hearing_loss_variation_ids_unique": qc_a.get("hearing_loss_variation_ids_unique", ""),
            "highconf_plp_tier_ge_2_variation_ids_unique": qc_a.get("highconf_plp_tier_ge_2_variation_ids_unique", ""),
            "unique_genes_hearing_loss": qc_a.get("unique_genes_hearing_loss", ""),
        },
        {
            "run_label": args.label_b,
            "hearing_loss_rows_record_level": qc_b.get("hearing_loss_rows_record_level", ""),
            "hearing_loss_variation_ids_unique": qc_b.get("hearing_loss_variation_ids_unique", ""),
            "highconf_plp_tier_ge_2_variation_ids_unique": qc_b.get("highconf_plp_tier_ge_2_variation_ids_unique", ""),
            "unique_genes_hearing_loss": qc_b.get("unique_genes_hearing_loss", ""),
        },
        {
            "run_label": f"top{args.top_n}_gene_jaccard",
            "hearing_loss_rows_record_level": "",
            "hearing_loss_variation_ids_unique": "",
            "highconf_plp_tier_ge_2_variation_ids_unique": f"{jaccard(genes_a, genes_b):.3f}",
            "unique_genes_hearing_loss": "",
        },
    ]
    write_tsv(Path(args.out), rows)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
