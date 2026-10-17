#!/usr/bin/env python3
"""Analyze four record selections from fixed ClinVar, HPO and ClinGen inputs."""
import argparse
import csv
import gzip
import json
import re
import subprocess
import sys
import time
import tempfile
from collections import Counter
from functools import lru_cache
from datetime import datetime, timezone
from pathlib import Path

from build_evidence_map import export_tables, new_entry, update_entry, parse_gene_annotations, sha256_path, write_tsv
from phenotype_filters import compile_pattern
from analysis_inputs import ROOT, RAW_SHA


@lru_cache(maxsize=32768)
def normalize_name(text):
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


@lru_cache(maxsize=4096)
def canonical_ids(text):
    return frozenset(("ORPHA" if namespace.lower() in {"orphanet", "orpha"} else namespace.upper()) + ":" + number
            for namespace, number in re.findall(r"(OMIM|Orphanet|ORPHA|HP|DECIPHER):([0-9]+)", text, re.I))


SPILL_FIELDS = ["VariationID", "#AlleleID", "GeneSymbol", "GeneID", "ClinicalSignificance", "ReviewStatus", "PhenotypeList",
                "PhenotypeIDS", "RCVaccession", "Assembly", "NumberSubmitters", "SCVsForAggregateGermlineClassification", "inclusion_reasons"]
PARTITIONS = 32


def consolidate_partitions(parts, outdir, total, matched, elapsed):
    """Merge disjoint VariationID partitions with bounded memory."""
    outdir.mkdir(parents=True, exist_ok=True)
    gene_counts = {}; metrics = Counter(); statuses = Counter(); unknown = Counter()
    for part in parts:
        for row in read_tsv(part / "gene_evidence_map.tsv"):
            c = gene_counts.setdefault(row["gene_symbol"], Counter())
            for field, value in row.items():
                if field != "gene_symbol": c[field] += int(value)
        for row in read_tsv(part / "qc_summary.tsv"):
            key = row["metric"]
            if key not in {"variant_summary_sha256", "total_rows_processed", "hearing_loss_rows_record_level", "elapsed_seconds", "unique_genes_hearing_loss", "plp_review_stars_ge_2_unique_genes"}:
                metrics[key] += int(row["value"])
        for row in read_tsv(part / "review_status_top50.tsv"): statuses[row["review_status_normalized"]] += int(row["n"])
        for row in read_tsv(part / "unknown_review_statuses.tsv"): unknown[row["review_status"]] += int(row["n_variants"])
    for name in ["variant_aggregate.tsv", "variant_high_confidence.tsv", "variant_gene_links.tsv", "unresolved_gene_annotations.tsv"]:
        temp = outdir / (name + ".partial")
        with temp.open("w", newline="") as output:
            writer = None
            for part in parts:
                with (part / name).open() as source:
                    reader = csv.DictReader(source, delimiter="\t")
                    if writer is None:
                        writer = csv.DictWriter(output, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n"); writer.writeheader()
                    writer.writerows(reader)
        temp.replace(outdir / name)
    gene_rows = [{"gene_symbol": g, **dict(c)} for g, c in sorted(gene_counts.items())]
    write_tsv(outdir / "gene_evidence_map.tsv", gene_rows)
    final = {"variant_summary_sha256": RAW_SHA, "total_rows_processed": total, "hearing_loss_rows_record_level": matched,
             "unique_genes_hearing_loss": len(gene_counts), "plp_review_stars_ge_2_unique_genes": sum(c["n_plp_review_stars_ge_2"] > 0 for c in gene_counts.values()),
             "elapsed_seconds": round(elapsed, 2), **dict(metrics)}
    write_tsv(outdir / "qc_summary.tsv", [{"metric": k, "value": v} for k, v in final.items()], fieldnames=["metric", "value"])
    write_tsv(outdir / "review_status_top50.tsv", [{"review_status_normalized": k, "n": n} for k, n in statuses.most_common()], fieldnames=["review_status_normalized", "n"])
    write_tsv(outdir / "unknown_review_statuses.tsv", [{"review_status": k, "n_variants": n} for k, n in unknown.items()], fieldnames=["review_status", "n_variants"])
    return final


def read_tsv(path):
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))



def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modes", nargs="+", default=["primary", "expanded", "disease_supported", "curated_gene"], choices=["primary", "expanded", "disease_supported", "curated_gene"])
    ap.add_argument("--outdir", type=Path, default=ROOT / "results")
    ap.add_argument("--references", type=Path, default=ROOT / "data/references")
    args = ap.parse_args()
    raw = ROOT / "data/raw/clinvar/variant_summary.txt.gz"
    if sha256_path(raw) != RAW_SHA:
        raise RuntimeError("Pinned ClinVar checksum mismatch; refusing to run")
    args.outdir = args.outdir.resolve()
    reference_checksums = {}
    if any(mode in args.modes for mode in {"disease_supported", "curated_gene"}):
        manifest = json.loads((args.references / "references_manifest.json").read_text())
        for name, metadata in manifest["outputs"].items():
            checksum = sha256_path(args.references / name)
            if checksum != metadata["sha256"]: raise RuntimeError(f"Reference checksum mismatch: {name}")
            reference_checksums[name] = checksum
    diseases, names, genes, hp_terms = set(), set(), set(), set()
    if "disease_supported" in args.modes:
        for r in read_tsv(args.references / "hearing_diseases.tsv"):
            diseases.update(canonical_ids(r["disease_id"]))
            names.add(normalize_name(r["disease_name"]))
            hp_terms.update(canonical_ids(r["hpo_terms"]))
        terms_file = args.references / "hearing_hpo_terms.tsv"
        if terms_file.exists():
            hp_terms.update(r["hpo_id"] for r in read_tsv(terms_file))
    if "curated_gene" in args.modes:
        positive = {"Definitive", "Strong", "Moderate", "Limited", "Supportive"}
        genes = {r["gene_symbol"] for r in read_tsv(args.references / "hearing_genes.tsv") if r["gene_symbol"] and r["classification"] in positive}
        if not genes:
            raise RuntimeError("No curated comparator genes")
    primary, expanded = compile_pattern("hearing_loss_only"), compile_pattern("expanded_synonyms")
    records = {mode: set() for mode in args.modes}
    matched = dict.fromkeys(args.modes, 0)
    logdir = ROOT / "logs"; logdir.mkdir(parents=True, exist_ok=True)
    log = logdir / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_analysis.log")
    start = time.time(); total = 0
    code_sha = {p.name: sha256_path(p) for p in Path(__file__).parent.glob("*.py")}
    with tempfile.TemporaryDirectory(prefix="clinvar-analysis-") as work, log.open("w", encoding="utf-8") as lf, gzip.open(raw, "rt", encoding="utf-8") as f:
        workdir = Path(work); handles = {}; writers = {}
        for mode in args.modes:
            for partition in range(PARTITIONS):
                path = workdir / f"{mode}_{partition}.tsv"
                handle = path.open("w", newline="")
                handles[(mode, partition)] = handle
                writer = csv.DictWriter(handle, fieldnames=SPILL_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader(); writers[(mode, partition)] = writer
        def report(msg):
            print(msg, flush=True); lf.write(msg + "\n"); lf.flush()
        report(f"sha256={RAW_SHA}; modes={args.modes}; references={args.references}")
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            text = row["PhenotypeList"]
            is_primary = bool(primary.search(text)); is_expanded = bool(expanded.search(text))
            selected = {"primary": {"original_phrase"} if is_primary else set(), "expanded": {"expanded_text"} if is_expanded else set()}
            if "disease_supported" in records:
                ids = canonical_ids(row["PhenotypeIDS"])
                row_names = {normalize_name(x) for x in text.split("|")}
                reasons = {"expanded_text"} if is_expanded else set()
                if ids & diseases: reasons.add("HPO_disease_id")
                if row_names & names: reasons.add("HPO_disease_name_exact")
                if ids & hp_terms: reasons.add("HPO_hearing_term")
                selected["disease_supported"] = reasons
            if "curated_gene" in records:
                symbols, _, _ = parse_gene_annotations(row.get("GeneSymbol", ""), row.get("GeneID", ""))
                selected["curated_gene"] = {"curated_gene_annotation"} if symbols & genes else set()
            for mode, ids_in_mode in records.items():
                reasons = selected[mode]
                if not reasons: continue
                vid = row["VariationID"].strip()
                if not vid: continue
                matched[mode] += 1
                ids_in_mode.add(vid)
                writers[(mode, int(vid) % PARTITIONS)].writerow({k: row.get(k, "") for k in SPILL_FIELDS} | {"inclusion_reasons": "|".join(sorted(reasons))})
            if total % 1_000_000 == 0:
                report(f"rows={total}; variants={ {m: len(v) for m,v in records.items()} }; elapsed={time.time()-start:.1f}s")
        for handle in handles.values(): handle.close()
        summaries = []
        for mode in records:
            part_dirs = []
            for partition in range(PARTITIONS):
                per_vid = {}
                with (workdir / f"{mode}_{partition}.tsv").open() as source:
                    for row in csv.DictReader(source, delimiter="\t"):
                        vid = row["VariationID"]
                        if vid not in per_vid: per_vid[vid] = new_entry(row)
                        update_entry(per_vid[vid], row, set(row["inclusion_reasons"].split("|")))
                if not per_vid: continue
                partdir = workdir / f"{mode}_tables_{partition}"
                export_tables(per_vid, partdir, RAW_SHA, total, 0, 0)
                part_dirs.append(partdir)
                del per_vid
            if not part_dirs: raise RuntimeError(f"Empty analysis definition: {mode}")
            metrics = consolidate_partitions(part_dirs, args.outdir / mode, total, matched[mode], time.time() - start)
            summaries.append({"definition": mode, **metrics})
            table_script = Path(__file__).with_name("generate_utility_tables.py")
            for single in [False, True]:
                out = args.outdir / mode / "single_gene" if single else args.outdir / mode
                cmd = [sys.executable, str(table_script), "--variant-aggregate", str(args.outdir / mode / "variant_aggregate.tsv"), "--outdir", str(out)]
                if single: cmd.append("--single-gene-only")
                subprocess.run(cmd, check=True, stdout=lf, stderr=lf)
            report(f"completed {mode}: {len(records[mode])} unique variants")
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "definition_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    pairs = []
    for left in records:
        for right in records:
            a, b = set(records[left]), set(records[right])
            pairs.append({"left": left, "right": right, "n_left": len(a), "n_right": len(b), "n_intersection": len(a & b), "n_left_only": len(a - b), "n_right_only": len(b - a)})
    write_tsv(args.outdir / "definition_overlaps.tsv", pairs)
    metadata = {"started_utc": datetime.fromtimestamp(start, timezone.utc).isoformat(), "source_sha256": RAW_SHA,
                "total_rows": total, "definitions": args.modes, "log": str(log.relative_to(ROOT)), "python": sys.version,
                "references": str(args.references.relative_to(ROOT)), "reference_sha256": reference_checksums, "code_sha256": code_sha, "aggregation_partitions": PARTITIONS, "elapsed_seconds": time.time() - start}
    (args.outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
