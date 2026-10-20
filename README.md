# Record coverage and disease context in hearing-loss-related ClinVar records

Analysis code and results for comparisons of text-based, HPO-linked and ClinGen gene-based ClinVar record selection. The analysis uses fixed public inputs and summarizes germline aggregate classifications, official review stars and annotated gene memberships.

## Reproduce the analysis

Python 3.11 or later is required; the recorded analysis environment used Python 3.13. Install the scientific dependencies and run the analysis entrypoint:

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 scripts/reproduce.py
```

The entrypoint downloads the fixed `Pinned_Inputs.zip` release asset when source files are absent and verifies every input against `docs/DATA_MANIFEST.tsv`. With a local copy:

```sh
python3 scripts/reproduce.py --inputs-archive /path/to/Pinned_Inputs.zip
```

It derives the phenotype and gene references from their original source files, computes all four record selections, summarizes selection overlap and gene profiles, parses the four saved VCV examples, and exports all five figures and reported tables. It requires several GB of free disk space and can take tens of minutes on a laptop. The rolling current ClinVar FTP file is not interchangeable with the fixed input.

## Outputs

- `plots/fig1.png` through `plots/fig5.png`, plus PDF/TIFF and plotting data.
- `tables/Table1.tsv`: record coverage under the four selection definitions.
- `tables/Table2.tsv`: disease contexts and classification evidence among HPO-linked additions.
- `tables/Table3.tsv`: review-status composition under the phrase rule.
- `tables/Supplementary_Tables.xlsx`: S1–S15, a Readme and a field Codebook.
- `results/`: selection summaries, gene profiles, comparisons and individual-submission examples.

To export tables or figures from existing result files:

```sh
python3 scripts/hearing_loss_clinvar/export_tables.py --main-workbook
python3 scripts/hearing_loss_clinvar/plot_figures.py
```
The study-design figure distinguishes three nested phenotype selections from the separate ClinGen gene comparator. It can be plotted individually with `python3 scripts/hearing_loss_clinvar/plot_figures.py --figures 1`.

The table exporter also creates `Main_Tables.xlsx` when `--main-workbook` is specified. `results/disease_context_display.tsv` and its provenance JSON identify the selected records, source hashes and nonadditive disease contributions behind Table 2. `docs/OUTPUT_PROVENANCE.tsv` maps every main figure and main/supporting table to its generating code and input results. Table 2 combines a single-record example, the complete four-star group and three overlapping disease-identifier contributions; its rows are not mutually exclusive categories.

## Inputs and interpretation

The fixed ClinVar input was accessed on 2026-05-23. HPO release v2026-09-01, a separately dated ClinGen gene-validity export and four saved VCV examples are specified in the input manifest. These references do not reconstruct historical classifications at every submission level. Record selection uses condition and gene annotations; a retrieved aggregate classification is not necessarily specific to hearing loss. Variant–gene membership is many-to-many.

Large raw inputs and complete record-level outputs are distributed as release assets. `Reproducibility.zip` contains the scientific code, all derived result files, tables and figures. The repository tree contains the code, summary results, tables and figures; large intermediate and record-level result files can be regenerated or obtained from that archive. Checksums accompany the release.

## Citation and license

Software citation metadata are provided in `CITATION.cff`. Code is available under the MIT license. Source data retain their respective terms; see `docs/DATA_MANIFEST.tsv` and `docs/DATA_LICENSES.md`. We acknowledge ClinVar and ClinGen contributors and the Human Phenotype Ontology Consortium.

The curated-gene comparator counts every explicitly named gene on each selected record. Its 1,484 annotated symbols comprise 164 qualifying reference genes and 1,320 coannotated symbols outside that reference. Figure 4 distinguishes highlighted groups by both color and marker shape; Figure 5 ranks annotations in phrase-selected records, whose disease labels can differ or overlap.
