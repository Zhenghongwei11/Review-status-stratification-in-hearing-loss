#!/usr/bin/env python3
"""Prepare versioned HPO disease and ClinGen gene-disease references.

Downloads are opt-in;
normal runs verify the recorded checksums before deriving tables.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[2]
HPO_RELEASE = "v2026-09-01"
HPO_BASE = f"https://github.com/obophenotype/human-phenotype-ontology/releases/download/{HPO_RELEASE}"
SOURCES = {
    "hp.obo": (f"{HPO_BASE}/hp.obo", HPO_RELEASE, "HPO custom license; attribution, version disclosure, integrity conditions"),
    "phenotype.hpoa": (f"{HPO_BASE}/phenotype.hpoa", HPO_RELEASE, "HPO custom license; attribution, version disclosure, integrity conditions"),
    "gencc_submissions.csv": ("https://thegencc.org/download/action/submissions-export-csv?format=new", "live export; retrieval date and checksum identify snapshot", "CC0 1.0 Universal"),
    "hpo_license.md": (f"https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/{HPO_RELEASE}/LICENSE.md", HPO_RELEASE, "HPO license pointer"),
    "hpo_legacy_license.html": ("https://human-phenotype-ontology.github.io/license.html", "retrieval-time official legacy license page", "HPO custom license"),
    "gencc_terms.html": ("https://thegencc.org/terms", "retrieval-time terms page", "CC0 1.0 Universal"),
    "clingen_gene_validity.csv": ("https://search.clinicalgenome.org/kb/gene-validity/download", "live export; retrieval date and checksum identify snapshot", "CC0 1.0 Universal"),
    "clingen_terms.html": ("https://clinicalgenome.org/about/terms-of-use/", "retrieval-time terms page", "CC0 1.0 Universal"),
}
ROOT_HPO = "HP:0000365"
DISEASE_TEXT = re.compile(r"\bhearing\b|\bdeaf(?:ness)?\b|\b(?:usher|pendred|waardenburg)\b", re.I)
POSITIVE_CLASSIFICATIONS = {"Definitive", "Strong", "Moderate", "Limited", "Supportive"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, path: Path) -> dict:
    """curl handles local proxy settings; atomically retain only complete downloads."""
    partial = path.with_suffix(path.suffix + ".part")
    headers = path.with_suffix(path.suffix + ".headers")
    previous_headers = headers.read_text() if headers.exists() else ""
    proxy = os.environ.get("REVISION_REFERENCE_PROXY")
    command = ["curl", "--http1.1", "-fL", "--retry", "2", "--connect-timeout", "30",
               "--max-time", "1800", *(["--proxy", proxy] if proxy else []), "-sS"]
    head = subprocess.run([*command, "--head", "--dump-header", str(headers),
                           "--write-out", "%{url_effective}", url, "-o", os.devnull],
                          check=True, capture_output=True, text=True)
    response_headers = headers.read_text()
    lengths = re.findall(r"^content-length:\s*(\d+)\s*$", response_headers, re.I | re.M)
    validators = re.findall(r"^(?:etag|last-modified):\s*(.+?)\s*$", response_headers, re.I | re.M)
    previous_validators = re.findall(r"^(?:etag|last-modified):\s*(.+?)\s*$", previous_headers, re.I | re.M)
    size = int(lengths[-1]) if lengths else 0
    # Range-capable official static exports can be downloaded in independently
    # checked pieces. This avoids losing a large file to a slow proxy timeout.
    if path.name != "gencc_submissions.csv" and size > 2 * 1024 * 1024 and "accept-ranges: bytes" in response_headers.lower():
        chunk_size = 2 * 1024 * 1024
        ranges = [(start, min(start + chunk_size - 1, size - 1)) for start in range(0, size, chunk_size)]
        def piece(bounds: tuple[int, int]) -> Path:
            start, end = bounds
            chunk = path.with_suffix(path.suffix + f".part.{start}")
            if validators and validators == previous_validators and chunk.exists() and chunk.stat().st_size == end - start + 1:
                return chunk
            suffix = ["--range", f"{start}-{end}",
                      *(["--header", f"If-Range: {validators[-1]}"] if validators else []),
                      "--write-out", "%{http_code}", head.stdout.strip(), "-o", str(chunk)]
            # Public release assets may be directly reachable even when the
            # normal browser/GitHub proxy is slow. Keep a bounded fallback.
            if urlsplit(head.stdout.strip()).hostname == "release-assets.githubusercontent.com":
                result = subprocess.run(["curl", "--http1.1", "--noproxy", "*", "-fL", "-sS",
                                         "--connect-timeout", "5", "--max-time", "30", *suffix],
                                        capture_output=True, text=True)
                if result.returncode:
                    result = subprocess.run([*command, *suffix], check=True, capture_output=True, text=True)
            else:
                result = subprocess.run([*command, *suffix], check=True, capture_output=True, text=True)
            if result.stdout.strip() != "206" or chunk.stat().st_size != end - start + 1:
                raise RuntimeError(f"Range request did not return exact source bytes: {path.name} {start}-{end}")
            return chunk
        with ThreadPoolExecutor(max_workers=4) as pool:
            pieces = list(pool.map(piece, ranges))
        with partial.open("wb") as handle:
            for chunk in pieces:
                handle.write(chunk.read_bytes())
                chunk.unlink()
        resolved_url = head.stdout.strip()
    else:
        result = subprocess.run([*command, "--dump-header", str(headers), "--write-out", "%{url_effective}",
                                 url, "-o", str(partial)], check=True, capture_output=True, text=True)
        resolved_url = result.stdout.strip()
    partial.replace(path)
    info = {"requested_url": url, "resolved_url": resolved_url,
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "sha256": sha256(path), "bytes": path.stat().st_size}
    path.with_suffix(path.suffix + ".download.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


def source_info(path: Path, *, refresh: bool = False) -> dict:
    url, version, license_name = SOURCES[path.name]
    meta = path.with_suffix(path.suffix + ".download.json")
    if refresh or not path.exists():
        info = download(url, path)
    elif meta.exists():
        info = json.loads(meta.read_text())
        if sha256(path) != info["sha256"]:
            raise RuntimeError(f"Cached source checksum changed: {path}")
    else:
        # Supports independently downloaded official source files; report the
        # file timestamp rather than pretend HTTP metadata were captured.
        info = {"requested_url": url, "resolved_url": None,
                "retrieved_at_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
                "retrieval_metadata": "pre-existing local download; timestamp from file mtime; no response headers captured",
                "sha256": sha256(path), "bytes": path.stat().st_size}
        meta.write_text(json.dumps(info, indent=2) + "\n")
    if path.stat().st_size == 0:
        raise RuntimeError(f"Empty source file is unavailable, not a valid cache: {path}; rerun --download")
    return {**info, "version": version, "license": license_name, "local_file": str(path.relative_to(ROOT))}


def parse_obo(path: Path) -> tuple[dict, str]:
    terms, current, version = {}, None, ""
    for line in path.read_text().splitlines():
        if line.startswith("data-version: "):
            version = line.split(": ", 1)[1]
        if line == "[Term]":
            current = {"parents": set(), "alt_ids": set(), "obsolete": False}
        elif line.startswith("["):
            current = None
        elif current is not None:
            if line.startswith("id: "):
                current["id"] = line[4:]
                terms[current["id"]] = current
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("is_a: "):
                current["parents"].add(line[6:].split()[0])
            elif line.startswith("alt_id: "):
                current["alt_ids"].add(line[8:])
            elif line == "is_obsolete: true":
                current["obsolete"] = True
    if ROOT_HPO not in terms or len(terms) < 10000:
        raise RuntimeError("Incomplete HPO ontology, or hearing-impairment root absent")
    return terms, version


def descendants(terms: dict) -> set[str]:
    children = defaultdict(set)
    for term_id, item in terms.items():
        if not item["obsolete"]:
            for parent in item["parents"]:
                children[parent].add(term_id)
    result, pending = {ROOT_HPO}, [ROOT_HPO]
    while pending:
        for child in children[pending.pop()]:
            if child not in result:
                result.add(child)
                pending.append(child)
    return result


def curie(value: str) -> str:
    value = value.strip()
    return re.sub(r"^(?:MIM|OMIM):", "OMIM:",
                  re.sub(r"^(?:ORPHA|Orphanet|ORPHANET):", "ORPHA:", value))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(source_dir: Path, outdir: Path, download_sources: bool, hpo_only: bool = False, gene_source: str = "clingen") -> dict:
    source_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    source_names = [name for name in SOURCES if not name.startswith(("gencc", "clingen")) or (not hpo_only and name.startswith(gene_source))]
    with ThreadPoolExecutor(max_workers=3) as pool:
        infos = list(pool.map(lambda name: source_info(source_dir / name, refresh=download_sources), source_names))
    sources = dict(zip(source_names, infos))
    terms, ontology_version = parse_obo(source_dir / "hp.obo")
    hearing_terms = descendants(terms)
    alt_ids = {alt: term for term, item in terms.items() for alt in item["alt_ids"]}
    comments, annotation_lines = [], []
    with (source_dir / "phenotype.hpoa").open() as handle:
        for line in handle:
            (comments if line.startswith("#") else annotation_lines).append(line)
    disease_records, annotation_count, negative_hearing_count = {}, 0, 0
    for row in csv.DictReader(annotation_lines, delimiter="\t"):
        annotation_count += 1
        hpo = alt_ids.get(row["hpo_id"], row["hpo_id"])
        if hpo not in hearing_terms:
            continue
        if "NOT" in row["qualifier"].split("|"):
            negative_hearing_count += 1
            continue
        if row.get("aspect", "P") != "P":
            continue
        disease_id = curie(row["database_id"])
        record = disease_records.setdefault(disease_id, {"disease_name": row["disease_name"], "hpo_terms": set()})
        record["hpo_terms"].add(hpo)
    if annotation_count < 100000 or len(disease_records) < 100:
        raise RuntimeError("Incomplete HPO annotation file; refuse plausible-looking partial reference")
    diseases = [{"disease_id": d, "disease_name": r["disease_name"],
                 "hpo_terms": "|".join(sorted(r["hpo_terms"])), "source": "HPO phenotype.hpoa"}
                for d, r in sorted(disease_records.items())]
    write_tsv(outdir / "hearing_diseases.tsv", diseases, ["disease_id", "disease_name", "hpo_terms", "source"])
    write_tsv(outdir / "hearing_hpo_terms.tsv", [{"hpo_id": t, "name": terms[t]["name"], "is_root": int(t == ROOT_HPO)} for t in sorted(hearing_terms)], ["hpo_id", "name", "is_root"])
    genes, match_counts, total_gencc = [], Counter(), 0
    if not hpo_only and gene_source == "gencc":
        with (source_dir / "gencc_submissions.csv").open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"gene_curie", "gene_symbol", "disease_curie", "disease_title", "classification_title", "submitter_title"}
            if not required.issubset(reader.fieldnames or []):
                raise RuntimeError("Unexpected GenCC new-format CSV schema")
            for row in reader:
                total_gencc += 1
                ids = {curie(row.get(key, "")) for key in ("disease_curie", "disease_original_curie", "submitted_as_disease_id")} - {""}
                id_hits = ids.intersection(disease_records)
                names = [row.get(key, "") for key in ("disease_title", "disease_original_title", "submitted_as_disease_name")]
                name_hits = [name for name in names if DISEASE_TEXT.search(name)]
                if not id_hits and not name_hits:
                    continue
                basis = "HPO_disease_id" if id_hits else "explicit_disease_name"
                match_counts[basis] += 1
                genes.append({"gene_symbol": row["gene_symbol"], "gene_id": row["gene_curie"],
                              "disease_id": curie(row["disease_curie"]), "disease_name": row["disease_title"],
                              "classification": row["classification_title"], "submitter": row["submitter_title"],
                              "match_basis": basis, "matched_disease_ids": "|".join(sorted(id_hits)),
                              "matched_disease_names": "|".join(sorted(set(name_hits))),
                              "source_submission_id": row.get("sgc_id", ""),
                              "source_submission_version": row.get("version_number", ""),
                              "original_disease_id": curie(row.get("disease_original_curie", "")),
                              "source": "GenCC"})
    if not hpo_only and gene_source == "gencc" and total_gencc < 10000:
        raise RuntimeError("Incomplete GenCC export; refuse partial reference")
    total_clingen = 0
    clingen_creation_header = None
    if not hpo_only and gene_source == "clingen":
        with (source_dir / "clingen_gene_validity.csv").open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            for row in reader:
                if row and row[0].startswith("FILE CREATED:"):
                    clingen_creation_header = row[0]
                if row and row[0] == "GENE SYMBOL":
                    fields = row
                    break
            else:
                raise RuntimeError("ClinGen CSV header absent")
            for values in reader:
                if not values or values[0].startswith("+"):
                    continue
                if len(values) != len(fields):
                    raise RuntimeError("Truncated or unexpected ClinGen CSV row")
                row = dict(zip(fields, values))
                total_clingen += 1
                disease_id = curie(row["DISEASE ID (MONDO)"])
                id_hit = disease_id in disease_records
                name_hit = bool(DISEASE_TEXT.search(row["DISEASE LABEL"]))
                panel_hit = "hearing loss" in row["GCEP"].lower()
                if not (id_hit or name_hit or panel_hit):
                    continue
                basis = "HPO_disease_id" if id_hit else "hearing_loss_expert_panel" if panel_hit else "explicit_disease_name"
                match_counts[basis] += 1
                genes.append({"gene_symbol": row["GENE SYMBOL"], "gene_id": row["GENE ID (HGNC)"],
                              "disease_id": disease_id, "disease_name": row["DISEASE LABEL"],
                              "classification": row["CLASSIFICATION"], "submitter": row["GCEP"],
                              "match_basis": basis, "matched_disease_ids": disease_id if id_hit else "",
                              "matched_disease_names": row["DISEASE LABEL"] if name_hit else "",
                              "source_submission_id": row["ONLINE REPORT"],
                              "source_submission_version": row["CLASSIFICATION DATE"],
                              "original_disease_id": "", "source": "ClinGen"})
        if total_clingen < 1000:
            raise RuntimeError("Incomplete ClinGen CSV; refuse partial reference")
    genes.sort(key=lambda row: (row["gene_symbol"], row["disease_id"], row["submitter"], row["source_submission_id"]))
    write_tsv(outdir / "hearing_diseases.tsv", diseases, ["disease_id", "disease_name", "hpo_terms", "source"])
    if not hpo_only:
        write_tsv(outdir / "hearing_genes.tsv", genes, ["gene_symbol", "gene_id", "disease_id", "disease_name", "classification", "submitter", "match_basis", "matched_disease_ids", "matched_disease_names", "source_submission_id", "source_submission_version", "original_disease_id", "source"])
    write_tsv(outdir / "hearing_hpo_terms.tsv", [{"hpo_id": t, "name": terms[t]["name"], "is_root": int(t == ROOT_HPO)} for t in sorted(hearing_terms)], ["hpo_id", "name", "is_root"])
    output_names = ["hearing_diseases.tsv", "hearing_hpo_terms.tsv"] + ([] if hpo_only else ["hearing_genes.tsv"])
    manifest = {"schema_version": 1, "reference_status": "HPO_complete_gene_curation_pending" if hpo_only else "complete_with_ClinGen_GenCC_unavailable" if gene_source == "clingen" else "complete", "gene_curation_source": None if hpo_only else gene_source, "clinvar_snapshot_accessed": "2026-05-23",
                "temporal_limitation": "HPO September 2026 and later gene-disease curations are applied to fixed May 2026 ClinVar records; these are not May 2026 annotations or an exhaustive gold standard.",
                "sources": sources, "hpo_ontology_data_version": ontology_version,
                "hpo_annotation_metadata": [line.rstrip() for line in comments],
                "selection": {"hpo_root": ROOT_HPO, "descendants_include_root": True,
                              "graph_edges": "is_a", "exclude_obsolete_terms": True,
                              "exclude_negative_NOT_annotations": True, "annotation_aspect": "P",
                              "gencc_id_match": "disease_curie OR disease_original_curie OR submitted_as_disease_id in HPO-selected disease IDs",
                              "gencc_name_pattern": DISEASE_TEXT.pattern,
                              "gencc_scope": "Only matched gene-disease submissions; other diseases for the same gene are not automatically included.",
                              "clingen_selection": "HPO disease CURIE match OR explicit hearing/deaf/Usher/Pendred/Waardenburg disease name OR Hearing Loss Gene Curation Expert Panel; only selected relationships retained",
                              "gene_id_namespace": "HGNC; not NCBI GeneID",
                              "positive_comparator_classifications": sorted(POSITIVE_CLASSIFICATIONS)},
                "counts": {"ontology_terms": len(terms), "hearing_hpo_terms_including_root": len(hearing_terms),
                           "hpo_annotation_rows": annotation_count, "negative_hearing_annotation_rows_excluded": negative_hearing_count,
                           "hearing_disease_ids": len(diseases), "gencc_input_submission_rows": total_gencc,
                           "hearing_gencc_submission_rows": len(genes), "hearing_gencc_unique_gene_symbols": len({g["gene_symbol"] for g in genes}),
                           "positive_evidence_gene_symbols": len({g["gene_symbol"] for g in genes if g["classification"] in POSITIVE_CLASSIFICATIONS}),
                           "gencc_match_basis_counts": dict(match_counts), "classification_counts": dict(Counter(g["classification"] for g in genes))},
                "outputs": {name: {"sha256": sha256(outdir / name), "bytes": (outdir / name).stat().st_size} for name in output_names}}
    if hpo_only:
        for field in list(manifest["counts"]):
            if "gencc" in field or field in {"classification_counts", "positive_evidence_gene_symbols"}:
                manifest["counts"][field] = None
    elif gene_source == "clingen":
        for field in list(manifest["counts"]):
            if "gencc" in field:
                manifest["counts"][field] = None
        manifest["counts"].update({"clingen_input_curation_rows": total_clingen,
                                  "hearing_clingen_curation_rows": len(genes),
                                  "hearing_clingen_unique_gene_symbols": len({g["gene_symbol"] for g in genes}),
                                  "clingen_match_basis_counts": dict(match_counts)})
        manifest["clingen_creation_header"] = clingen_creation_header
        manifest["unavailable_sources"] = [{"source": "GenCC", "requested_url": SOURCES["gencc_submissions.csv"][0],
                                            "retrieval_date": "2026-10-13", "status": "HTTP 429 Too Many Requests",
                                            "access_limit": "official download page reports 20 requests/day; no further requests made after rate limitation was identified",
                                            "partial_downloads_used": False,
                                            "replacement": "complete official ClinGen gene-disease validity CSV; not equivalent to the full multi-submitter GenCC database"}]
    (outdir / "references_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Refresh official sources; omit for checksum-checked offline derivation")
    parser.add_argument("--hpo-only", action="store_true", help="Derive complete HPO references while explicitly leaving GenCC pending")
    parser.add_argument("--gene-source", choices=("clingen", "gencc"), default="clingen", help="ClinGen is the available official alternative; GenCC requires a complete export and observing its download limit")
    parser.add_argument("--source-dir", type=Path, default=ROOT / "data/sources")
    parser.add_argument("--outdir", type=Path, default=ROOT / "data/references")
    args = parser.parse_args()
    manifest = prepare(args.source_dir.resolve(), args.outdir.resolve(), args.download, args.hpo_only, args.gene_source)
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
