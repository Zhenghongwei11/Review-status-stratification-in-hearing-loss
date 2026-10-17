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
from pathlib import Path
from typing import TextIO


from phenotype_filters import compile_pattern


DEFAULT_PHENOTYPE_MODE = "hearing_loss_only"


def normalize_review_status(text: str) -> str:
    return (text or "").strip().lower()


def review_tier(review_status: str) -> int:
    s = normalize_review_status(review_status)
    # Official germline review statuses, not a clinical confidence score.
    # Unknown statuses remain distinguishable from recognized zero-star statuses.
    return {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, single submitter": 1,
        "criteria provided, conflicting classifications": 1,
        "criteria provided, conflicting interpretations": 1,
        "no assertion criteria provided": 0,
        "no classification provided": 0,
        "no assertion provided": 0,
        "no classifications from unflagged records": 0,
    }.get(s, -1)


GENE_SUMMARY_RE = re.compile(r"^covers\s+\d+\s+genes?\b", re.IGNORECASE)
GENE_SUBSET_RE = re.compile(r"^subset\s+of\s+\d+\s+genes?\s*:\s*(.+)$", re.IGNORECASE)
GENE_SYMBOL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@-]*$")


def parse_gene_annotation_details(symbol_text: str, id_text: str) -> dict:
    """Parse ClinVar GeneSymbol conservatively and retain an audit trail.

    ClinVar structural-variant rows can contain prose summaries such as
    ``covers 17 genes`` or ``subset of 41 genes: PKD1:TSC2``.  Those strings
    are not gene symbols. Explicit symbols after ``subset`` are retained
    independently of whether the source supplies usable, paired GeneIDs.
    The stated total never supplies unnamed members of a structural variant.
    """
    raw = (symbol_text or "").strip()
    ids = [s.strip() for s in re.split(r"[;,]", id_text or "") if s.strip() and s.strip() not in {"-1", "-"}]
    symbols: list[str] = []
    status = "missing"
    if not raw or raw == "-":
        status = "missing"
    elif GENE_SUMMARY_RE.match(raw):
        status = "unresolved_summary_label"
    else:
        subset = GENE_SUBSET_RE.match(raw)
        if subset:
            candidates = [s.strip() for s in subset.group(1).split(":") if s.strip()]
            if candidates and all(GENE_SYMBOL_RE.fullmatch(s) for s in candidates):
                symbols = candidates
                status = "resolved_explicit_subset" if len(candidates) == len(ids) else "resolved_explicit_subset_unpaired_ids"
            else:
                status = "unresolved_non_symbol_subset"
        else:
            candidates = [s.strip() for s in raw.split(";") if s.strip() and s.strip() != "-"]
            if candidates and all(GENE_SYMBOL_RE.fullmatch(s) for s in candidates):
                symbols = candidates
                status = "resolved"
            else:
                status = "unresolved_non_symbol_label"
    # Only a one-to-one source order can support symbol-to-ID pairing.
    pairs = set(zip(symbols, ids)) if symbols and len(symbols) == len(ids) else {(s, "") for s in symbols}
    return {"raw_annotation": raw, "symbols": set(symbols), "ids": set(ids), "pairs": pairs,
            "parse_status": status, "unresolved": bool(raw and status.startswith("unresolved_"))}


def parse_gene_annotations(symbol_text: str, id_text: str) -> tuple[set[str], set[str], set[tuple[str, str]]]:
    details = parse_gene_annotation_details(symbol_text, id_text)
    return details["symbols"], details["ids"], details["pairs"]


def new_entry(row: dict) -> dict:
    return {"gene_symbols": set(), "gene_ids": set(), "gene_pairs": set(), "buckets": set(),
            "max_tier": -1, "best_review_status": "", "best_review_tier": -2,
            "allele_id": row.get("#AlleleID", ""), "phenotype_ids": row.get("PhenotypeIDS", ""),
            "phenotype_list": row.get("PhenotypeList", ""), "rcv_accession": row.get("RCVaccession", ""),
            "clinsig_example": row.get("ClinicalSignificance", ""), "source_rows": 0,
            "assemblies": set(), "classifications": set(), "review_statuses": set(),
            "condition_names": set(), "condition_ids": set(), "scvs": set(), "rcvs": set(),
            "submitter_counts": set(), "inclusion_reasons": set(), "classification_review_pairs": set(),
            "raw_gene_annotations": set(), "gene_parse_statuses": set(), "unresolved_gene_annotations": set()}


def update_entry(entry: dict, row: dict, reasons: set[str] | None = None) -> None:
    details = parse_gene_annotation_details(row.get("GeneSymbol", ""), row.get("GeneID", ""))
    symbols, ids, pairs = details["symbols"], details["ids"], details["pairs"]
    entry["gene_symbols"].update(symbols)
    entry["gene_ids"].update(ids)
    entry["gene_pairs"].update(pairs)
    if details["raw_annotation"]:
        entry["raw_gene_annotations"].add(details["raw_annotation"])
    entry["gene_parse_statuses"].add(details["parse_status"])
    if details["unresolved"]:
        entry["unresolved_gene_annotations"].add((details["raw_annotation"], details["parse_status"]))
    entry["buckets"].add(clin_sig_bucket(row.get("ClinicalSignificance", "")))
    rs = row.get("ReviewStatus", "")
    tier = review_tier(rs)
    entry["classification_review_pairs"].add((row.get("ClinicalSignificance", ""), tier))
    entry["max_tier"] = max(entry["max_tier"], tier)
    if tier > entry["best_review_tier"]:
        entry["best_review_tier"], entry["best_review_status"] = tier, rs
    entry["source_rows"] += 1
    for field, source in [("assemblies", "Assembly"), ("classifications", "ClinicalSignificance"), ("review_statuses", "ReviewStatus")]:
        entry[field].add(row.get(source, ""))
    entry["condition_names"].update(x.strip() for x in row.get("PhenotypeList", "").split("|") if x.strip())
    entry["condition_ids"].update(x.strip() for x in row.get("PhenotypeIDS", "").split("|") if x.strip())
    entry["scvs"].update(re.findall(r"SCV\d+", row.get("SCVsForAggregateGermlineClassification", "")))
    entry["rcvs"].update(re.findall(r"RCV\d+", row.get("RCVaccession", "")))
    n = row.get("NumberSubmitters", "")
    if n.isdigit():
        entry["submitter_counts"].add(int(n))
    entry["inclusion_reasons"].update(reasons or set())


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

    Fixed precedence for repeated preaggregated source rows (not SCV voting):
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


def export_tables(per_vid: dict, outdir: Path, source_sha: str, total_rows: int, matched_rows: int, elapsed: float) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    rows, links, unresolved_rows = [], [], []
    gene_counts = {}
    tiers, buckets, statuses = Counter(), Counter(), Counter()
    unknown = Counter()
    unknown_variants = 0
    for vid in sorted(per_vid, key=int):
        e = per_vid[vid]
        symbols = sorted(e["gene_symbols"])
        bucket, tier = aggregate_bucket(e["buckets"]), e["max_tier"]
        # Legacy primary fields are populated only for unambiguous singleton annotations.
        # All analyses use gene_symbols and the exported many-to-many membership table.
        gene = symbols[0] if len(symbols) == 1 else ""
        known_ids = sorted({i for g, i in e["gene_pairs"] if g == gene and i})
        joint_plp_ge_2 = any(
            clin_sig_bucket(classification) == "pathogenic_or_likely_pathogenic" and stars >= 2
            for classification, stars in e["classification_review_pairs"]
        )
        r = {"variation_id": vid, "primary_gene_symbol": gene,
             "primary_gene_id": ";".join(known_ids), "gene_symbols": ";".join(symbols),
             "gene_ids": ";".join(sorted(e["gene_ids"])), "gene_symbol_set_size": len(symbols),
             "raw_gene_annotations": " || ".join(sorted(e["raw_gene_annotations"])),
             "gene_parse_statuses": "|".join(sorted(e["gene_parse_statuses"])),
             "unresolved_gene_annotation_count": len(e["unresolved_gene_annotations"]),
             "bucket_variant_level": bucket, "review_tier_proxy_max": tier,
             "plp_review_stars_ge_2_same_row": int(joint_plp_ge_2),
             "review_status_best": e["best_review_status"], "allele_id_example": e["allele_id"],
             "phenotype_list_trunc": truncate_text(e["phenotype_list"], 500),
             "phenotype_ids_trunc": truncate_text(e["phenotype_ids"], 500),
             "rcv_accession_example": e["rcv_accession"], "clinsig_example": e["clinsig_example"],
             "source_row_count": e["source_rows"], "assemblies": "|".join(sorted(e["assemblies"])),
             "source_classifications": "|".join(sorted(e["classifications"])),
             "source_classification_count": len(e["classifications"]),
             "source_bucket_set": "|".join(sorted(e["buckets"])), "source_bucket_count": len(e["buckets"]),
             "source_review_statuses": "|".join(sorted(e["review_statuses"])),
             "source_classification_review_pairs": "|".join(f"{classification}@@{stars}" for classification, stars in sorted(e["classification_review_pairs"])),
             "condition_names": "|".join(sorted(e["condition_names"])),
             "condition_ids": "|".join(sorted(e["condition_ids"])),
             "germline_scv_ids": "|".join(sorted(e["scvs"])), "n_germline_scvs": len(e["scvs"]),
             "rcv_ids": "|".join(sorted(e["rcvs"])), "n_rcvs": len(e["rcvs"]),
             "number_submitters_min": min(e["submitter_counts"]) if e["submitter_counts"] else "",
             "number_submitters_max": max(e["submitter_counts"]) if e["submitter_counts"] else "",
             "inclusion_reasons": "|".join(sorted(e["inclusion_reasons"]))}
        rows.append(r)
        for raw_annotation, parse_status in sorted(e["unresolved_gene_annotations"]):
            unresolved_rows.append({"variation_id": vid, "raw_annotation": raw_annotation,
                                    "parse_status": parse_status})
        tiers[tier] += 1; buckets[bucket] += 1; statuses[e["best_review_status"]] += 1
        unknown_labels = [status for status in e["review_statuses"] if review_tier(status) < 0]
        unknown_variants += bool(unknown_labels)
        for status in unknown_labels:
            unknown[status] += 1
        for symbol in symbols:
            ids = sorted({i for g, i in e["gene_pairs"] if g == symbol and i})
            links.append({"variation_id": vid, "gene_symbol": symbol, "gene_id_unambiguous": ";".join(ids),
                          "n_annotated_genes": len(symbols), "bucket": bucket, "review_stars": tier})
            c = gene_counts.setdefault(symbol, Counter())
            c["n_variants"] += 1; c[bucket] += 1; c[f"tier_{tier}"] += 1
            if len(symbols) == 1: c["n_single_gene_variants"] += 1
            if joint_plp_ge_2: c["n_higher_review_plp"] += 1
    fields = list(rows[0]) if rows else ["variation_id", "gene_symbols", "bucket_variant_level", "review_tier_proxy_max"]
    high = [r for r in rows if r["plp_review_stars_ge_2_same_row"]]
    write_tsv(outdir / "variant_aggregate.tsv", rows, fieldnames=fields)
    write_tsv(outdir / "variant_high_confidence.tsv", high, fieldnames=fields)
    write_tsv(outdir / "variant_gene_links.tsv", links, fieldnames=["variation_id", "gene_symbol", "gene_id_unambiguous", "n_annotated_genes", "bucket", "review_stars"])
    write_tsv(outdir / "unresolved_gene_annotations.tsv", unresolved_rows,
              fieldnames=["variation_id", "raw_annotation", "parse_status"])
    gene_rows = []
    for symbol, c in sorted(gene_counts.items()):
        gene_rows.append({"gene_symbol": symbol, "n_variants_hearing_loss": c["n_variants"],
            "n_plp": c["pathogenic_or_likely_pathogenic"], "n_vus": c["vus"], "n_conflicting": c["conflicting"],
            "n_benign": c["benign_or_likely_benign"], "n_other": c["other"] + c["unknown"],
            "n_single_gene_variants": c["n_single_gene_variants"], "n_multi_gene_variants": c["n_variants"] - c["n_single_gene_variants"],
            "n_plp_review_stars_ge_2": c["n_higher_review_plp"], **{f"n_stars_{i}": c[f"tier_{i}"] for i in range(-1, 5)}})
    write_tsv(outdir / "gene_evidence_map.tsv", gene_rows, fieldnames=list(gene_rows[0]) if gene_rows else ["gene_symbol", "n_variants_hearing_loss"])
    metrics = {"variant_summary_sha256": source_sha, "total_rows_processed": total_rows,
        "hearing_loss_rows_record_level": matched_rows, "hearing_loss_variation_ids_unique": len(rows),
        "plp_review_stars_ge_2_variation_ids_unique": len(high), "unique_genes_hearing_loss": len(gene_counts),
        "plp_review_stars_ge_2_unique_genes": len({g for r in high for g in r["gene_symbols"].split(";") if g}),
        "hearing_loss_variants_multi_gene_symbol": sum(r["gene_symbol_set_size"] > 1 for r in rows),
        "gene_variant_memberships": len(links),
        "variants_with_unresolved_gene_annotations": sum(r["unresolved_gene_annotation_count"] > 0 for r in rows),
        "unresolved_gene_annotation_memberships": len(unresolved_rows),
        "variants_differing_source_classifications": sum(r["source_classification_count"] > 1 for r in rows),
        "variants_differing_source_buckets": sum(r["source_bucket_count"] > 1 for r in rows),
        "conflicting_variants_stars_ge_2": sum(r["bucket_variant_level"] == "conflicting" and r["review_tier_proxy_max"] >= 2 for r in rows),
        "unknown_review_status_variants": unknown_variants, "elapsed_seconds": round(elapsed, 2)}
    metrics.update({f"bucket_count::{k}": v for k, v in sorted(buckets.items())})
    metrics.update({f"review_tier_proxy_count::{k}": v for k, v in sorted(tiers.items())})
    write_tsv(outdir / "qc_summary.tsv", [{"metric": k, "value": v} for k, v in metrics.items()], fieldnames=["metric", "value"])
    write_tsv(outdir / "review_status_top50.tsv", [{"review_status_normalized": k, "n": v} for k, v in statuses.most_common()], fieldnames=["review_status_normalized", "n"])
    write_tsv(outdir / "unknown_review_statuses.tsv", [{"review_status": k, "n_variants": v} for k, v in unknown.items()], fieldnames=["review_status", "n_variants"])
    return metrics


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize phenotype-selected ClinVar variants with official stars and all gene annotations.")
    ap.add_argument("--variant-summary", default="data/raw/clinvar/variant_summary.txt.gz")
    ap.add_argument("--outdir", default="results/primary")
    ap.add_argument("--phenotype-mode", default=DEFAULT_PHENOTYPE_MODE)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--progress-every", type=int, default=1_000_000)
    args = ap.parse_args()
    t0 = time.time(); per_vid = {}; total = matched = 0
    pattern = compile_pattern(args.phenotype_mode)
    with gzip.open(args.variant_summary, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if args.max_rows and total >= args.max_rows: break
            total += 1
            if args.progress_every and total % args.progress_every == 0:
                print(f"[progress] {total:,} rows; {len(per_vid):,} variants", file=sys.stderr, flush=True)
            if not pattern.search(row.get("PhenotypeList", "")): continue
            matched += 1
            vid = row.get("VariationID", "").strip()
            if not vid: continue
            update_entry(per_vid.setdefault(vid, new_entry(row)), row, {args.phenotype_mode})
    export_tables(per_vid, Path(args.outdir), sha256_path(Path(args.variant_summary)), total, matched, time.time() - t0)
    print(args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
