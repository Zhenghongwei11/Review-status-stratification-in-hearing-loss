#!/usr/bin/env python3
"""Retrieve separately dated VCV examples and retain actual SCV classifications."""
import argparse
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

from build_evidence_map import sha256_path, write_tsv, review_tier
from analysis_inputs import ROOT

URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
IDS = ["41608", "45216", "2353", "17035"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    args = ap.parse_args()
    source = ROOT / "data/sources/clinvar_examples_vcv.xml"
    url = URL + "?" + urllib.parse.urlencode({"db": "clinvar", "rettype": "vcv", "is_variationid": "", "id": ",".join(IDS)})
    if args.download:
        source.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=90) as response:
            source.write_bytes(response.read())
    root = ET.parse(source).getroot()
    variants, submissions = [], []
    for v in root.findall("VariationArchive"):
        vid = v.get("VariationID")
        variants.append({"variation_id": vid, "vcv_accession": v.get("Accession"), "vcv_version": v.get("Version"),
            "name": v.get("VariationName"), "last_updated": v.get("DateLastUpdated"),
            "number_of_submissions": v.get("NumberOfSubmissions"), "number_of_submitters": v.get("NumberOfSubmitters"),
            "review_status": v.findtext("ClassifiedRecord/Classifications/GermlineClassification/ReviewStatus", ""),
            "classification": v.findtext("ClassifiedRecord/Classifications/GermlineClassification/Description", "")})
        for assertion in v.findall("ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion"):
            accession = assertion.find("ClinVarAccession")
            if accession is None:
                continue
            classification = assertion.find("Classification")
            if classification is None:
                continue
            status = classification.findtext("ReviewStatus", "")
            traits = assertion.findall("TraitSet/Trait")
            disease_names = sorted({e.text for t in traits for e in t.findall("Name/ElementValue") if e.text})
            disease_ids = sorted({e.get("DB", "") + ":" + e.get("ID", "") for t in traits for e in t.findall("XRef")})
            submissions.append({"variation_id": vid, "scv_accession": accession.get("Accession"), "scv_version": accession.get("Version"),
                "submitter": accession.get("SubmitterName"), "organization_category": accession.get("OrganizationCategory"),
                "classification": classification.findtext("GermlineClassification", ""), "review_status": status,
                "review_stars": review_tier(status), "date_last_evaluated": classification.get("DateLastEvaluated", ""),
                "date_updated": accession.get("DateUpdated", ""), "date_submitted": assertion.get("SubmissionDate", ""),
                "contributes_to_aggregate": assertion.get("ContributesToAggregateClassification", ""),
                "condition_names": "|".join(disease_names), "condition_ids": "|".join(disease_ids)})
    if {v["variation_id"] for v in variants} != set(IDS):
        raise RuntimeError("Requested examples were not all returned")
    out = ROOT / "results"
    write_tsv(out / "retrieved_vcv_examples.tsv", variants)
    write_tsv(out / "retrieved_scv_examples.tsv", submissions)
    manifest = {"url": url, "retrieved_utc": json.loads((ROOT / "data/references/scv_examples_manifest.json").read_text())["retrieved_utc"], "sha256": sha256_path(source),
        "source_path": str(source.relative_to(ROOT)), "variation_ids": IDS, "n_scvs": len(submissions),
        "snapshot_boundary": "Current VCV versions retrieved in October 2026; illustrative only. They do not reconstruct May 2026 SCV classifications. No population-level SCV analysis is inferred from four purposive examples.",
        "selection": "Selected conflict records in MET and DIAPH1; USH2A expert-panel syndromic omission; GJB2 zero-star conflict exception."}
    (ROOT / "data/references").mkdir(parents=True, exist_ok=True)
    (ROOT / "data/references/scv_examples_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Parsed {len(variants)} VCVs and {len(submissions)} SCVs")


if __name__ == "__main__":
    main()
