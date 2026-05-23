#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import io
import re
from pathlib import Path
import sys

import requests


DEFAULT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
HEARING_LOSS_RE = re.compile(r"\bhearing\s+loss\b", re.IGNORECASE)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Stream ClinVar variant_summary and write a small sample containing hearing-loss phenotype rows."
    )
    ap.add_argument("--url", default=DEFAULT_URL)
    ap.add_argument("--out", default="data/raw/clinvar/variant_summary.sample_hearing_loss.tsv.gz")
    ap.add_argument("--rows", type=int, default=20000, help="Number of hearing-loss rows to keep (plus header).")
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".partial")

    kept = 0
    with requests.get(args.url, stream=True, timeout=180) as r:
        r.raise_for_status()
        gz = gzip.GzipFile(fileobj=r.raw)
        text_stream = io.TextIOWrapper(gz, encoding="utf-8", errors="replace", newline="")

        header = text_stream.readline()
        if not header:
            raise RuntimeError("Empty stream (no header).")
        cols = header.strip("\r\n").split("\t")
        try:
            idx_pheno = cols.index("PhenotypeList")
        except ValueError as e:
            raise RuntimeError("PhenotypeList column not found in header.") from e

        scanned = 0
        with gzip.open(tmp, "wt", encoding="utf-8", newline="") as out_f:
            out_f.write(header)
            for line in text_stream:
                scanned += 1
                if scanned % 500000 == 0:
                    print(f"[progress] scanned={scanned} kept={kept}", file=sys.stderr)
                parts = line.rstrip("\\n").split("\\t")
                if idx_pheno >= len(parts):
                    continue
                pheno = parts[idx_pheno]
                if not pheno or HEARING_LOSS_RE.search(pheno) is None:
                    continue
                out_f.write(line)
                kept += 1
                if kept >= args.rows:
                    break

    if kept == 0:
        raise RuntimeError("No hearing-loss rows found; increase scan or check regex.")

    tmp.replace(out)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
