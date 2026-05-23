#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import sys
import time
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


from phenotype_filters import compile_pattern


DEFAULT_PHENOTYPE_MODE = "hearing_loss_only"


def normalize_review_status(text: str) -> str:
    return (text or "").strip().lower()


def review_tier(review_status: str) -> int:
    s = normalize_review_status(review_status)
    # Approximate star tiers using ClinVar review status wording.
    if "practice guideline" in s:
        return 4
    if "reviewed by expert panel" in s:
        return 3
    if "criteria provided, multiple submitters, no conflicts" in s:
        return 2
    if "criteria provided, single submitter" in s:
        return 1
    if "no assertion criteria provided" in s:
        return 0
    return 0


def clin_sig_bucket(clinsig: str) -> str:
    s = (clinsig or "").strip().lower()
    if not s:
        return "unknown"
    if "pathogenic" in s:
        # includes "likely pathogenic" and combined strings like "Pathogenic/Likely pathogenic"
        if "conflict" in s:
            return "conflicting"
        return "pathogenic_or_likely_pathogenic"
    if "benign" in s:
        if "conflict" in s:
            return "conflicting"
        return "benign_or_likely_benign"
    if "uncertain significance" in s:
        return "vus"
    if "conflicting interpretations" in s:
        return "conflicting"
    return "other"


def phenotype_mentions_hearing_loss(phenotype_list: str) -> bool:
    if not phenotype_list:
        return False
    # Backward-compatible fallback. Prefer the mode-driven pattern.
    return re.search(r"\bhearing\s+loss\b", phenotype_list, re.IGNORECASE) is not None


@dataclass
class GeneAgg:
    gene_id: str
    gene_symbol: str
    n_variants: int = 0
    n_plp: int = 0
    n_benign: int = 0
    n_vus: int = 0
    n_conflicting: int = 0
    n_other: int = 0
    n_tier4_practice_guideline: int = 0
    n_tier3_expert_panel: int = 0
    n_tier2_multi_submitters_no_conflicts: int = 0
    n_tier1_single_submitter: int = 0
    n_tier0_no_criteria: int = 0
    max_review_tier: int = 0

    def update(self, bucket: str, tier: int) -> None:
        self.n_variants += 1
        if bucket == "pathogenic_or_likely_pathogenic":
            self.n_plp += 1
        elif bucket == "benign_or_likely_benign":
            self.n_benign += 1
        elif bucket == "vus":
            self.n_vus += 1
        elif bucket == "conflicting":
            self.n_conflicting += 1
        elif bucket == "other":
            self.n_other += 1
        else:
            self.n_other += 1

        if tier == 4:
            self.n_tier4_practice_guideline += 1
        elif tier == 3:
            self.n_tier3_expert_panel += 1
        elif tier == 2:
            self.n_tier2_multi_submitters_no_conflicts += 1
        elif tier == 1:
            self.n_tier1_single_submitter += 1
        else:
            self.n_tier0_no_criteria += 1
        self.max_review_tier = max(self.max_review_tier, tier)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict], *, fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        raise RuntimeError(f"No rows to write for {path} (and no fieldnames provided).")
    if not fieldnames:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def tsv_writer(f: TextIO, fieldnames: list[str]) -> csv.DictWriter:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    w.writeheader()
    return w


def truncate_text(s: str, max_len: int) -> str:
    if s is None:
        return ""
    s = str(s)
    if len(s) <= max_len:
        return s
    return s[: max(0, max_len - 3)] + "..."


def aggregate_bucket(buckets: set[str]) -> str:
    """
    Aggregate multiple row-level buckets into a single variant-level bucket.

    Conservative ordering:
      conflicting > pathogenic/likely pathogenic > VUS > benign/likely benign > other/unknown
    """
    if "conflicting" in buckets:
        return "conflicting"
    if "pathogenic_or_likely_pathogenic" in buckets:
        return "pathogenic_or_likely_pathogenic"
    if "vus" in buckets:
        return "vus"
    if "benign_or_likely_benign" in buckets:
        return "benign_or_likely_benign"
    if "other" in buckets:
        return "other"
    if "unknown" in buckets:
        return "unknown"
    return "other"


def choose_primary_gene(gene_symbols: set[str]) -> str:
    cleaned = sorted({g.strip() for g in gene_symbols if g and g.strip()})
    return cleaned[0] if cleaned else ""


def main() -> int:
    ap = argparse.ArgumentParser(description="Build hearing-loss ClinVar evidence-map tables from variant_summary.txt.gz.")
    ap.add_argument(
        "--variant-summary",
        default="data/raw/clinvar/variant_summary.txt.gz",
        help="Path to ClinVar variant_summary.txt.gz",
    )
    ap.add_argument(
        "--outdir",
        default="results/hearing_loss_clinvar",
        help="Output directory for tables",
    )
    ap.add_argument(
        "--max-rows",
        type=int,
        default=0,
        help="If >0, only process this many rows (debug/smoke).",
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Emit progress to stderr every N processed rows (0 disables).",
    )
    ap.add_argument(
        "--phenotype-mode",
        default=DEFAULT_PHENOTYPE_MODE,
        help="Phenotype filter mode for ClinVar PhenotypeList (e.g., hearing_loss_only, expanded_synonyms).",
    )
    args = ap.parse_args()

    vs_path = Path(args.variant_summary)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    total_rows = 0
    hl_rows = 0
    hl_genes: set[str] = set()
    all_genes: set[str] = set()

    review_status_counter = Counter()
    bucket_counter = Counter()
    tier_counter = Counter()

    gene_aggs: dict[tuple[str, str], GeneAgg] = {}

    phenotype_re = compile_pattern(args.phenotype_mode)

    # ClinVar variant_summary repeats the same VariationID across assemblies (and other dimensions),
    # so we must de-duplicate at the VariationID level before counting.
    #
    # We only retain per-VariationID aggregates for variants that match the phenotype filter at least once.
    per_vid: dict[str, dict] = {}
    hl_variation_ids: set[str] = set()
    multi_gene_variants = 0
    highconf_n = 0

    with gzip.open(vs_path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"GeneID", "GeneSymbol", "ClinicalSignificance", "PhenotypeList", "ReviewStatus", "VariationID", "RCVaccession"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"variant_summary missing columns: {sorted(missing)}")

        for row in reader:
            total_rows += 1
            if args.max_rows and total_rows > args.max_rows:
                break
            if args.progress_every and (total_rows % args.progress_every == 0):
                elapsed = time.time() - t0
                rate = total_rows / elapsed if elapsed > 0 else 0.0
                print(
                    f"[progress] rows={total_rows:,} hearing_loss_rows_record={hl_rows:,} unique_variation_ids={len(hl_variation_ids):,} rate={rate:,.0f} rows/s",
                    file=sys.stderr,
                    flush=True,
                )

            gene_id = (row.get("GeneID") or "").strip()
            gene_symbol = (row.get("GeneSymbol") or "").strip()
            if gene_symbol:
                all_genes.add(gene_symbol)

            phenotype_list = row.get("PhenotypeList") or ""
            if phenotype_re.search(phenotype_list) is None:
                continue

            hl_rows += 1

            vid = (row.get("VariationID") or "").strip()
            if not vid:
                continue
            hl_variation_ids.add(vid)
            if gene_symbol:
                hl_genes.add(gene_symbol)

            clinsig = row.get("ClinicalSignificance") or ""
            bucket = clin_sig_bucket(clinsig)
            rs = row.get("ReviewStatus") or ""
            tier = review_tier(rs)

            entry = per_vid.get(vid)
            if entry is None:
                entry = {
                    "gene_symbols": set(),
                    "gene_ids": set(),
                    "buckets": set(),
                    "max_tier": tier,
                    "best_review_status": rs,
                    "best_review_tier": tier,
                    "allele_id": (row.get("#AlleleID") or "").strip(),
                    "phenotype_ids": (row.get("PhenotypeIDS") or "").strip(),
                    "phenotype_list": phenotype_list,
                    "rcv_accession": (row.get("RCVaccession") or "").strip(),
                    "clinsig_example": clinsig,
                }
                per_vid[vid] = entry

            if gene_symbol:
                entry["gene_symbols"].add(gene_symbol)
            if gene_id:
                entry["gene_ids"].add(gene_id)
            entry["buckets"].add(bucket)
            if tier > entry["max_tier"]:
                entry["max_tier"] = tier
            if tier > entry["best_review_tier"] or (tier == entry["best_review_tier"] and not entry.get("best_review_status") and rs):
                entry["best_review_tier"] = tier
                entry["best_review_status"] = rs

    # Variant-level aggregation (unique VariationID).
    # We export:
    #  - `variant_aggregate.tsv`: all unique HL-matching VariationIDs (analysis backbone for utility endpoints)
    #  - `variant_high_confidence.tsv`: P/LP AND tier>=2 proxy (clinic-facing "reportable-now" shortlist)
    var_fields = [
        "variation_id",
        "primary_gene_symbol",
        "primary_gene_id",
        "bucket_variant_level",
        "review_tier_proxy_max",
        "review_status_best",
        "allele_id_example",
        "phenotype_list_trunc",
        "phenotype_ids_trunc",
        "rcv_accession_example",
        "clinsig_example",
        "gene_symbol_set_size",
    ]
    agg_path = outdir / "variant_aggregate.tsv"
    highconf_path = outdir / "variant_high_confidence.tsv"
    agg_path.parent.mkdir(parents=True, exist_ok=True)
    highconf_path.parent.mkdir(parents=True, exist_ok=True)

    with agg_path.open("w", newline="", encoding="utf-8") as agg_f, highconf_path.open(
        "w", newline="", encoding="utf-8"
    ) as highconf_f:
        agg_w = tsv_writer(agg_f, var_fields)
        highconf_w = tsv_writer(highconf_f, var_fields)

        for vid, entry in per_vid.items():
            bucket_v = aggregate_bucket(entry["buckets"])
            tier_max = int(entry["max_tier"] or 0)
            rs_best = entry.get("best_review_status") or ""

            gene_symbols: set[str] = entry["gene_symbols"]
            gene_ids: set[str] = entry["gene_ids"]
            primary_symbol = choose_primary_gene(gene_symbols)
            primary_gene_id = sorted({g for g in gene_ids if g and g.strip()})[0] if gene_ids else ""
            if len(gene_symbols) > 1:
                multi_gene_variants += 1

            row_out = {
                "variation_id": vid,
                "primary_gene_symbol": primary_symbol,
                "primary_gene_id": primary_gene_id,
                "bucket_variant_level": bucket_v,
                "review_tier_proxy_max": str(tier_max),
                "review_status_best": rs_best,
                "allele_id_example": entry.get("allele_id", ""),
                "phenotype_list_trunc": truncate_text(entry.get("phenotype_list", ""), 500),
                "phenotype_ids_trunc": truncate_text(entry.get("phenotype_ids", ""), 500),
                "rcv_accession_example": truncate_text(entry.get("rcv_accession", ""), 200),
                "clinsig_example": truncate_text(entry.get("clinsig_example", ""), 200),
                "gene_symbol_set_size": str(len(gene_symbols)),
            }
            agg_w.writerow(row_out)

            bucket_counter[bucket_v] += 1
            tier_counter[str(tier_max)] += 1
            review_status_counter[normalize_review_status(rs_best)] += 1

            key = (primary_gene_id, primary_symbol)
            if primary_symbol:
                if key not in gene_aggs:
                    gene_aggs[key] = GeneAgg(gene_id=primary_gene_id, gene_symbol=primary_symbol)
                gene_aggs[key].update(bucket=bucket_v, tier=tier_max)

            if bucket_v == "pathogenic_or_likely_pathogenic" and tier_max >= 2:
                highconf_w.writerow(row_out)
                highconf_n += 1

    gene_rows: list[dict] = []
    for (gene_id, gene_symbol), agg in sorted(gene_aggs.items(), key=lambda kv: (-(kv[1].n_plp + kv[1].n_vus), kv[0][1])):
        if not gene_symbol:
            continue
        gene_rows.append(
            {
                "gene_symbol": gene_symbol,
                "gene_id": gene_id,
                "n_variants_hearing_loss": str(agg.n_variants),
                "n_plp": str(agg.n_plp),
                "n_benign": str(agg.n_benign),
                "n_vus": str(agg.n_vus),
                "n_conflicting": str(agg.n_conflicting),
                "n_other": str(agg.n_other),
                "max_review_tier_proxy": str(agg.max_review_tier),
                "n_tier4_practice_guideline": str(agg.n_tier4_practice_guideline),
                "n_tier3_expert_panel": str(agg.n_tier3_expert_panel),
                "n_tier2_multi_submitters_no_conflicts": str(agg.n_tier2_multi_submitters_no_conflicts),
                "n_tier1_single_submitter": str(agg.n_tier1_single_submitter),
                "n_tier0_no_criteria_or_other": str(agg.n_tier0_no_criteria),
            }
        )

    qc_rows = [
        {"metric": "variant_summary_sha256", "value": sha256_path(vs_path)},
        {"metric": "total_rows_processed", "value": str(total_rows)},
        {"metric": "hearing_loss_rows_record_level", "value": str(hl_rows)},
        {"metric": "hearing_loss_variation_ids_unique", "value": str(len(hl_variation_ids))},
        {"metric": "highconf_plp_tier_ge_2_variation_ids_unique", "value": str(highconf_n)},
        {"metric": "unique_genes_all", "value": str(len(all_genes))},
        {"metric": "unique_genes_hearing_loss", "value": str(len(hl_genes))},
        {"metric": "hearing_loss_variants_multi_gene_symbol", "value": str(multi_gene_variants)},
        {"metric": "elapsed_seconds", "value": str(int(time.time() - t0))},
    ]
    for k, v in bucket_counter.most_common():
        qc_rows.append({"metric": f"bucket_count::{k}", "value": str(v)})
    for k, v in tier_counter.most_common():
        qc_rows.append({"metric": f"review_tier_proxy_count::{k}", "value": str(v)})

    gene_fields = [
        "gene_symbol",
        "gene_id",
        "n_variants_hearing_loss",
        "n_plp",
        "n_benign",
        "n_vus",
        "n_conflicting",
        "n_other",
        "max_review_tier_proxy",
        "n_tier4_practice_guideline",
        "n_tier3_expert_panel",
        "n_tier2_multi_submitters_no_conflicts",
        "n_tier1_single_submitter",
        "n_tier0_no_criteria_or_other",
    ]
    qc_fields = ["metric", "value"]

    write_tsv(outdir / "gene_evidence_map.tsv", gene_rows, fieldnames=gene_fields)
    write_tsv(outdir / "qc_summary.tsv", qc_rows, fieldnames=qc_fields)

    # Also write a compact review-status distribution (top 50) for debugging.
    rs_rows = [{"review_status_normalized": k, "n": str(v)} for k, v in review_status_counter.most_common(50)]
    write_tsv(outdir / "review_status_top50.tsv", rs_rows, fieldnames=["review_status_normalized", "n"])

    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
