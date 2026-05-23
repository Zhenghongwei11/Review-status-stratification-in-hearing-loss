#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import sys
from pathlib import Path

import requests


DEFAULT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_with_resume(url: str, dest: Path, *, timeout_s: int = 120) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")

    headers = {}
    existing = tmp.stat().st_size if tmp.exists() else 0
    if existing:
        headers["Range"] = f"bytes={existing}-"

    with requests.get(url, stream=True, headers=headers, timeout=timeout_s) as r:
        if r.status_code not in (200, 206):
            raise RuntimeError(f"Download failed: HTTP {r.status_code} for {url}")
        mode = "ab" if existing and r.status_code == 206 else "wb"
        with tmp.open(mode) as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

    tmp.replace(dest)


def main() -> int:
    ap = argparse.ArgumentParser(description="Download ClinVar tab-delimited files with checksum.")
    ap.add_argument("--url", default=DEFAULT_URL, help="ClinVar variant_summary.txt.gz URL")
    ap.add_argument(
        "--out",
        default="data/raw/clinvar/variant_summary.txt.gz",
        help="Output path for the downloaded file",
    )
    ap.add_argument(
        "--write-sha256",
        default="data/raw/clinvar/variant_summary.txt.gz.sha256",
        help="Write SHA256 checksum to this path",
    )
    args = ap.parse_args()

    out = Path(args.out)
    print(f"[download] {args.url} -> {out}")
    download_with_resume(args.url, out)

    digest = sha256_file(out)
    sha_path = Path(args.write_sha256)
    sha_path.parent.mkdir(parents=True, exist_ok=True)
    sha_path.write_text(digest + "\n", encoding="utf-8")
    print(f"[sha256] {digest}")

    stamp = dt.datetime.now().replace(microsecond=0).isoformat()
    meta_path = out.with_suffix(out.suffix + ".meta.txt")
    meta_path.write_text(f"downloaded_at_local={stamp}\nurl={args.url}\nsha256={digest}\n", encoding="utf-8")
    print(f"[meta] {meta_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

