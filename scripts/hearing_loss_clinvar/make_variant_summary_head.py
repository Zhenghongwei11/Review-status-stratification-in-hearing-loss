#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import io
from pathlib import Path

import requests


DEFAULT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"


def main() -> int:
    ap = argparse.ArgumentParser(description="Create a small head subset of ClinVar variant_summary for smoke tests.")
    ap.add_argument("--url", default=DEFAULT_URL)
    ap.add_argument("--out", default="data/raw/clinvar/variant_summary.head.tsv.gz")
    ap.add_argument("--lines", type=int, default=200000)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Stream download, decompress, keep first N lines, then re-compress.
    with requests.get(args.url, stream=True, timeout=120) as r:
        r.raise_for_status()
        gz = gzip.GzipFile(fileobj=r.raw)
        text_stream = io.TextIOWrapper(gz, encoding="utf-8", errors="replace", newline="")
        with gzip.open(out, "wt", encoding="utf-8", newline="") as out_f:
            for i, line in enumerate(text_stream, start=1):
                out_f.write(line)
                if i >= args.lines:
                    break

    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

