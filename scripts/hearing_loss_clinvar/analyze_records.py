#!/usr/bin/env python3
"""Stream paired definitions, archive evidence and gene-annotation sensitivity.

Only VariationID sets and gene-sized tables persist in memory; variant rows and
additional-variant tables are never materialized as rich dictionaries/lists.
"""
import csv
import hashlib
import json
from collections import Counter
from contextlib import ExitStack

from build_evidence_map import aggregate_bucket, clin_sig_bucket, write_tsv
from analysis_inputs import ROOT

OUT = ROOT / "results"
REF = ROOT / "data/references"
MODES = ("primary", "expanded", "disease_supported", "curated_gene")
ADDITION_FIELDS = ["definition", "variation_id", "gene_symbols", "bucket_variant_level", "review_tier_proxy_max", "condition_names", "condition_ids", "inclusion_reasons"]
STRUCTURE_FIELDS = ["definition", "variation_id", "source_row_count", "source_classification_count", "source_bucket_count", "n_germline_scvs", "number_submitters_min", "number_submitters_max", "n_rcvs"]
VALIDITY_FIELDS = ["definition", "gene_symbol", "n_annotated_variants", "curation_match", "curation_sources", "classifications", "disease_ids", "submitters", "n_reference_curation_rows"]
SENSITIVITY_FIELDS = ["definition", "gene_symbol", "n_all_annotation_variants", "n_single_gene_variants", "n_all_conflicting", "n_single_gene_conflicting", "n_all_uncertainty", "n_single_gene_uncertainty", "uncertainty_proportion", "uncertainty_count_rank", "uncertainty_proportion_rank_n_ge_10"]
POLICIES = ("fixed_precedence", "require_same_bucket_else_other", "highest_review_source_rows")
HASH_MODULUS = 1 << 256


def rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


class TsvSink:
    """Publish a streamed table atomically only after the analysis succeeds."""
    def __init__(self, path, fields):
        self.path, self.temporary, self.fields = path, path.with_suffix(path.suffix + ".tmp"), fields

    def __enter__(self):
        self.handle = self.temporary.open("w", encoding="utf-8", newline="")
        self.writer = csv.DictWriter(self.handle, fieldnames=self.fields, delimiter="\t", lineterminator="\n")
        self.writer.writeheader()
        return self.writer

    def __exit__(self, kind, value, traceback):
        self.handle.close()
        if kind is None:
            self.temporary.replace(self.path)
        else:
            self.temporary.unlink(missing_ok=True)


def pair_hash(vid, gene):
    # Compare every variant-gene membership without storing millions of pairs.
    # Link uniqueness is checked separately, so duplicates cannot be hidden.
    return int.from_bytes(hashlib.sha256(f"{vid}\0{gene}".encode()).digest(), "big")


def summarize_stream(mode, primary_ids, additions, expert_additions, structure):
    ids, genes, review_status_genes = set(), set(), set()
    counts, changed, memberships_by_gene = Counter(), Counter(), Counter()
    differing = mixed = membership_count = membership_hash = 0
    for row in rows(OUT / mode / "variant_aggregate.tsv"):
        vid = row["variation_id"]
        if vid in ids:
            raise RuntimeError(f"Duplicate VariationID in final table: {mode}/{vid}")
        ids.add(vid)
        symbols = {gene for gene in row["gene_symbols"].split(";") if gene}
        if len(symbols) != int(row["gene_symbol_set_size"]):
            raise RuntimeError(f"Gene-set-size mismatch: {mode}/{vid}")
        genes.update(symbols)
        membership_count += len(symbols)
        for gene in symbols:
            memberships_by_gene[gene] += 1
            membership_hash = (membership_hash + pair_hash(vid, gene)) % HASH_MODULUS
        bucket, tier = row["bucket_variant_level"], int(row["review_tier_proxy_max"])
        counts["multi"] += len(symbols) > 1
        counts["single"] += len(symbols) == 1
        counts["no_gene"] += not symbols
        counts["unavailable"] += tier < 0
        counts["with_germline_scvs"] += int(row["n_germline_scvs"]) > 0
        counts["sum_germline_scv_counts"] += int(row["n_germline_scvs"])
        counts["expert"] += tier == 3
        counts["guideline"] += tier == 4
        counts["vus"] += bucket == "vus"
        counts["conflicting"] += bucket == "conflicting"
        if int(row["plp_review_stars_ge_2_same_row"]):
            counts["plp_review_stars_ge_2"] += 1
            review_status_genes.update(symbols)
        if mode != "primary" and vid not in primary_ids:
            addition = {"definition": mode, **{field: row[field] for field in ADDITION_FIELDS[1:]}}
            additions.writerow(addition)
            if mode == "disease_supported" and tier == 3:
                expert_additions.writerow(addition)
        buckets = set(row["source_bucket_set"].split("|"))
        differing += int(row["source_classification_count"]) > 1
        mixed += len(buckets) > 1
        highest_tier, highest_buckets = -1, set()
        for pair in row["source_classification_review_pairs"].split("|"):
            classification, stars = pair.rsplit("@@", 1)
            stars = int(stars)
            if stars > highest_tier:
                highest_tier, highest_buckets = stars, {clin_sig_bucket(classification)}
            elif stars == highest_tier:
                highest_buckets.add(clin_sig_bucket(classification))
        alternatives = {"fixed_precedence": aggregate_bucket(buckets),
                        "require_same_bucket_else_other": next(iter(buckets)) if len(buckets) == 1 else "other",
                        "highest_review_source_rows": aggregate_bucket(highest_buckets)}
        for policy in POLICIES:
            changed[policy] += alternatives[policy] != bucket
        structure.writerow({"definition": mode, **{field: row[field] for field in STRUCTURE_FIELDS[1:]}})
    n = len(ids)
    summary = {"definition": mode, "n_variation_ids": n, "n_gene_symbols": len(genes),
               "n_multi_gene_variants": counts["multi"], "n_plp_review_stars_ge_2": counts["plp_review_stars_ge_2"],
               "plp_review_stars_ge_2_proportion": counts["plp_review_stars_ge_2"] / n if n else 0,
               "n_plp_review_stars_ge_2_gene_symbols": len(review_status_genes), "n_expert_panel": counts["expert"],
               "n_practice_guideline": counts["guideline"], "n_vus": counts["vus"], "n_conflicting": counts["conflicting"],
               "n_single_gene_variants": counts["single"], "n_no_gene_variants": counts["no_gene"],
               "n_unavailable_review_status": counts["unavailable"],
               "n_variants_with_participating_germline_scvs": counts["with_germline_scvs"],
               "sum_variant_participating_germline_scv_counts": counts["sum_germline_scv_counts"]}
    policies = [{"definition": mode, "policy": policy, "n_variants": n,
                 "n_differing_source_classifications": differing, "n_mixed_source_buckets": mixed,
                 "n_changed_variant_buckets": changed[policy], "scope": "preaggregated_summary_rows_not_SCV_votes"}
                for policy in POLICIES]
    return ids, summary, policies, (membership_count, membership_hash, memberships_by_gene)


def validate_memberships(mode, ids, expected, all_genes):
    expected_count, expected_hash, expected_by_gene = expected
    seen_variants, current_genes, current_vid = set(), set(), None
    link_count = link_hash = 0
    observed_by_gene = Counter()
    for link in rows(OUT / mode / "variant_gene_links.tsv"):
        vid, gene = link["variation_id"], link["gene_symbol"]
        if vid not in ids:
            raise RuntimeError(f"Unknown VariationID in gene links: {mode}/{vid}")
        # Export groups a variant's links contiguously, including partition
        # concatenation. Enforce that contract instead of storing every pair.
        if vid != current_vid:
            if vid in seen_variants:
                raise RuntimeError(f"Non-contiguous or repeated gene-link group: {mode}/{vid}")
            seen_variants.add(vid)
            current_vid, current_genes = vid, set()
        if gene in current_genes:
            raise RuntimeError(f"Duplicate variant-gene link: {mode}/{vid}/{gene}")
        current_genes.add(gene)
        link_count += 1
        observed_by_gene[gene] += 1
        link_hash = (link_hash + pair_hash(vid, gene)) % HASH_MODULUS
    gene_counts = {gene: int(row["n_unique_variants"]) for gene, row in all_genes.items()}
    if (link_count != expected_count or link_hash != expected_hash or observed_by_gene != expected_by_gene
            or gene_counts != dict(expected_by_gene)):
        raise RuntimeError(f"Variant-gene membership/count/checksum mismatch: {mode}")
    return {"definition": mode, "unique_variants": len(ids), "gene_memberships": expected_count,
            "gene_sum": sum(gene_counts.values()), "membership_check": "passed",
            "membership_sha256_sum_mod_2_256": f"{expected_hash:064x}"}


def gene_rows(mode, refs, all_genes, validity, sensitivity, anchors):
    single_genes = {row["gene_symbol"]: row for row in rows(OUT / mode / "single_gene/utility_gene_all.tsv")}
    count_order = sorted(all_genes, key=lambda gene: (-int(all_genes[gene]["uncertainty_burden"]), gene))
    prop_order = sorted((gene for gene in all_genes if int(all_genes[gene]["n_unique_variants"]) >= 10),
                        key=lambda gene: (-int(all_genes[gene]["uncertainty_burden"]) / int(all_genes[gene]["n_unique_variants"]), gene))
    count_ranks = {gene: i for i, gene in enumerate(count_order, 1)}
    prop_ranks = {gene: i for i, gene in enumerate(prop_order, 1)}
    for gene, row in all_genes.items():
        evidence = refs.get(gene, [])
        validity.writerow({"definition": mode, "gene_symbol": gene, "n_annotated_variants": row["n_unique_variants"],
                           "curation_match": "matched" if evidence else "not_matched_in_selected_reference",
                           "curation_sources": "|".join(sorted({item.get("source", "") for item in evidence})),
                           "classifications": "|".join(sorted({item["classification"] for item in evidence})),
                           "disease_ids": "|".join(sorted({item["disease_id"] for item in evidence})),
                           "submitters": "|".join(sorted({item["submitter"] for item in evidence})),
                           "n_reference_curation_rows": len(evidence)})
        single = single_genes.get(gene, {})
        result = {"definition": mode, "gene_symbol": gene,
                  "n_all_annotation_variants": row["n_unique_variants"], "n_single_gene_variants": single.get("n_unique_variants", 0),
                  "n_all_conflicting": row["n_conflicting"], "n_single_gene_conflicting": single.get("n_conflicting", 0),
                  "n_all_uncertainty": row["uncertainty_burden"], "n_single_gene_uncertainty": single.get("uncertainty_burden", 0),
                  "uncertainty_proportion": int(row["uncertainty_burden"]) / int(row["n_unique_variants"]),
                  "uncertainty_count_rank": count_ranks[gene], "uncertainty_proportion_rank_n_ge_10": prop_ranks.get(gene, "")}
        sensitivity.writerow(result)
        if mode == "primary" and gene in {"DIAPH1", "MET"}:
            anchors.writerow(result)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    refs = {}
    for row in rows(REF / "hearing_genes.tsv"):
        refs.setdefault(row["gene_symbol"], []).append(row)
    id_sets, coverage, policies, validations = {}, [], [], []
    with ExitStack() as stack:
        additions = stack.enter_context(TsvSink(OUT / "additional_variants_vs_primary.tsv", ADDITION_FIELDS))
        experts = stack.enter_context(TsvSink(OUT / "additional_disease_expert_variants.tsv", ADDITION_FIELDS))
        structure = stack.enter_context(TsvSink(OUT / "source_evidence_structure.tsv", STRUCTURE_FIELDS))
        validity = stack.enter_context(TsvSink(OUT / "gene_disease_validity_comparison.tsv", VALIDITY_FIELDS))
        sensitivity = stack.enter_context(TsvSink(OUT / "gene_annotation_sensitivity.tsv", SENSITIVITY_FIELDS))
        anchors = stack.enter_context(TsvSink(OUT / "figure_gene_anchors.tsv", SENSITIVITY_FIELDS))
        for mode in MODES:
            ids, summary, mode_policies, expected = summarize_stream(mode, id_sets.get("primary", set()), additions, experts, structure)
            id_sets[mode] = ids
            coverage.append(summary)
            policies.extend(mode_policies)
            all_genes = {row["gene_symbol"]: row for row in rows(OUT / mode / "utility_gene_all.tsv")}
            validations.append(validate_memberships(mode, ids, expected, all_genes))
            gene_rows(mode, refs, all_genes, validity, sensitivity, anchors)
    overlaps = []
    for left in MODES:
        for right in MODES:
            a, b = id_sets[left], id_sets[right]
            intersection = sum(vid in b for vid in a)
            overlaps.append({"left": left, "right": right, "n_left": len(a), "n_right": len(b),
                             "n_intersection": intersection, "n_left_only": len(a) - intersection,
                             "n_right_only": len(b) - intersection})
    write_tsv(OUT / "coverage_summary.tsv", coverage)
    write_tsv(OUT / "definition_overlaps.tsv", overlaps)
    write_tsv(OUT / "aggregation_policy_sensitivity.tsv", policies)
    write_tsv(OUT / "numerical_verification.tsv", validations)
    facts = {row["definition"]: row for row in coverage}
    (OUT / "verified_findings.json").write_text(json.dumps(facts, indent=2) + "\n")
    print(json.dumps(facts, indent=2))


if __name__ == "__main__":
    main()
