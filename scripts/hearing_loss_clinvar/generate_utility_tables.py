#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        yield from r


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def to_int(x: str) -> int:
    try:
        return int(x)
    except Exception:
        return 0


def safe_div(num: int, den: int) -> float:
    return (num / den) if den else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate descriptive classification and gene membership tables.")
    ap.add_argument("--outdir", default="results/hearing_loss_clinvar")
    ap.add_argument("--variant-aggregate", default="results/hearing_loss_clinvar/variant_aggregate.tsv")
    ap.add_argument("--gene-map", default="results/hearing_loss_clinvar/gene_evidence_map.tsv")
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--single-gene-only", action="store_true", help="Exclude variants annotated to multiple genes.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variants = read_tsv(Path(args.variant_aggregate))
    if args.single_gene_only:
        variants = (v for v in variants if int(v.get("gene_symbol_set_size", "1")) == 1)

    # Descriptive star-by-classification table; the denominator is VariationID records.
    by_tier_bucket: dict[tuple[str, str], int] = defaultdict(int)
    by_tier_total: dict[str, int] = defaultdict(int)
    by_tier_joint_plp: dict[str, int] = defaultdict(int)
    gene_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for v in variants:
        tier_text = (v.get("review_tier_proxy_max") or "").strip()
        bucket_text = (v.get("bucket_variant_level") or "").strip()
        by_tier_bucket[(tier_text, bucket_text)] += 1
        by_tier_total[tier_text] += 1
        by_tier_joint_plp[tier_text] += to_int(v.get("plp_review_stars_ge_2_same_row", "0"))
        for gene in sorted(set((v.get("gene_symbols") or v.get("primary_gene_symbol") or "").split(";"))):
            if not gene: continue
            bucket = (v.get("bucket_variant_level") or "").strip()
            tier = to_int((v.get("review_tier_proxy_max") or "").strip())

            gene_counts[gene]["n_unique_variants"] += 1
            if bucket == "pathogenic_or_likely_pathogenic":
                gene_counts[gene]["n_plp"] += 1
                if to_int(v.get("plp_review_stars_ge_2_same_row", "0")):
                    gene_counts[gene]["n_plp_review_stars_ge_2"] += 1
            elif bucket == "vus":
                gene_counts[gene]["n_vus"] += 1
            elif bucket == "conflicting":
                gene_counts[gene]["n_conflicting"] += 1
            elif bucket == "benign_or_likely_benign":
                gene_counts[gene]["n_benign"] += 1
            else:
                gene_counts[gene]["n_other"] += 1


    tier_bucket_rows: list[dict] = []
    tiers_sorted = sorted(by_tier_total.keys(), key=lambda t: to_int(t), reverse=True)
    buckets = [
        "pathogenic_or_likely_pathogenic",
        "vus",
        "conflicting",
        "benign_or_likely_benign",
        "other",
        "unknown",
    ]
    for tier in tiers_sorted:
        row = {"review_tier_proxy": tier, "n_variants": str(by_tier_total[tier])}
        for b in buckets:
            row[f"n::{b}"] = str(by_tier_bucket.get((tier, b), 0))
        # Descriptive P/LP subset defined by official ClinVar review stars.
        review_status_subset = by_tier_joint_plp.get(tier, 0)
        row["n_plp_review_stars_ge_2"] = str(review_status_subset)
        tier_bucket_rows.append(row)

    write_tsv(
        outdir / "utility_overall_by_tier.tsv",
        tier_bucket_rows,
        fieldnames=list(tier_bucket_rows[0].keys()) if tier_bucket_rows else ["review_tier_proxy", "n_variants"],
    )

    gene_rows: list[dict] = []
    for gene, c in gene_counts.items():
        n_total = c["n_unique_variants"]
        n_plp = c["n_plp"]
        n_review_status = c["n_plp_review_stars_ge_2"]
        n_vus = c["n_vus"]
        n_conf = c["n_conflicting"]
        gene_rows.append(
            {
                "gene_symbol": gene,
                "n_unique_variants": str(n_total),
                "n_plp": str(n_plp),
                "n_plp_review_stars_ge_2": str(n_review_status),
                "n_vus": str(n_vus),
                "n_conflicting": str(n_conf),
                "n_benign": str(c["n_benign"]),
                "n_other": str(c["n_other"]),
                "uncertainty_burden": str(n_vus + n_conf),
                "conflict_rate": f"{safe_div(n_conf, n_total):.4f}",
                "vus_rate": f"{safe_div(n_vus, n_total):.4f}",
                "plp_review_stars_ge_2_rate": f"{safe_div(n_review_status, n_total):.4f}",
            }
        )

    # Sort for two separate tables: actionability-focused and re-review hotspot-focused
    def sort_actionability(x: dict) -> tuple:
        return (
            -to_int(x["n_plp_review_stars_ge_2"]),
            -to_int(x["n_plp"]),
            -to_int(x["n_unique_variants"]),
            x["gene_symbol"],
        )

    def sort_hotspot(x: dict) -> tuple:
        return (-to_int(x["n_conflicting"]), -to_int(x["uncertainty_burden"]), -to_int(x["n_unique_variants"]), x["gene_symbol"])

    gene_rows_action = sorted(gene_rows, key=sort_actionability)[: max(1, args.top_n)]
    gene_rows_hot = sorted(gene_rows, key=sort_hotspot)[: max(1, args.top_n)]

    fields = list(gene_rows[0].keys()) if gene_rows else []
    if gene_rows:
        write_tsv(outdir / "utility_gene_all.tsv", gene_rows, fieldnames=fields)
    if gene_rows_action:
        write_tsv(outdir / f"utility_gene_actionability_top{args.top_n}.tsv", gene_rows_action, fieldnames=fields)
    if gene_rows_hot:
        write_tsv(outdir / f"utility_gene_conflict_hotspots_top{args.top_n}.tsv", gene_rows_hot, fieldnames=fields)

    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
