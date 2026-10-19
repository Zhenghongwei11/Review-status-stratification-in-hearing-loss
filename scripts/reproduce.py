#!/usr/bin/env python3
"""Reproduce the complete ClinVar analysis, five figures and reported tables."""
import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT_URL = 'https://github.com/Zhenghongwei11/Review-status-stratification-in-hearing-loss/releases/download/v2.0.0/Pinned_Inputs.zip'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def inputs(archive):
    with (ROOT / 'docs/DATA_MANIFEST.tsv').open() as source:
        expected = {row['path']: row for row in csv.DictReader(source, delimiter='\t')}
    missing = [name for name in expected if not (ROOT / name).is_file()]
    if missing:
        if archive is None:
            archive = ROOT / 'inputs-download.zip'
            temporary = archive.with_suffix('.zip.partial')
            print('Downloading the fixed input archive:', INPUT_URL, flush=True)
            with urllib.request.urlopen(INPUT_URL, timeout=120) as response, temporary.open('wb') as output:
                shutil.copyfileobj(response, output)
            temporary.replace(archive)
        with zipfile.ZipFile(archive) as bundle:
            if set(bundle.namelist()) != set(expected):
                raise RuntimeError('Input archive does not match the fixed input manifest')
            for name in missing:
                destination = ROOT / name
                destination.parent.mkdir(parents=True, exist_ok=True)
                temporary = destination.with_suffix(destination.suffix + '.partial')
                with bundle.open(name) as source, temporary.open('wb') as output:
                    shutil.copyfileobj(source, output)
                if sha256(temporary) != expected[name]['sha256']:
                    temporary.unlink()
                    raise RuntimeError('Input checksum mismatch: ' + name)
                temporary.replace(destination)
    for name, record in expected.items():
        if sha256(ROOT / name) != record['sha256']:
            raise RuntimeError('Input checksum mismatch: ' + name)
    print('All fixed input checksums verified.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs-archive', type=Path, help='Local Pinned_Inputs.zip; otherwise missing inputs are downloaded from the versioned release')
    args = parser.parse_args()
    inputs(args.inputs_archive)
    programs = [
        ('prepare_references.py', '--gene-source', 'clingen'),
        ('run_analysis.py',),
        ('analyze_records.py',),
        ('summarize_hpo_selection.py',),
        ('compare_gene_evidence.py',),
        ('parse_submission_examples.py',),
        ('plot_figures.py',),
        ('export_tables.py',),
    ]
    for name, *options in programs:
        print('Running:', name, flush=True)
        subprocess.run([sys.executable, str(ROOT / 'scripts/hearing_loss_clinvar' / name), *options], check=True, cwd=ROOT)
    print('Complete: plots/fig1–5, tables/Table1.tsv, tables/Table2.tsv, tables/Table3.tsv and tables/Supplementary_Tables.xlsx', flush=True)


if __name__ == '__main__':
    main()
