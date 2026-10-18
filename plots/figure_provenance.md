# Figure data and provenance

The five figures summarize fixed-input ClinVar record selections, review status and gene profiles. Each `plots/figN_source_data.tsv` provides the exact values plotted. `figure_provenance.json` distinguishes the recorded computational hashes from the public code layout and released metadata.

Figures use white backgrounds, an accessible discrete palette, TrueType PDF fonts and 300 dpi PNG/TIFF exports. Counts use distinct VariationIDs. A record can contribute to several named genes; gene totals are nonadditive. The P/LP subset requires both that classification and at least two official review stars on the same selected source row. Gene plots require at least ten annotated records.

Figure 4 labels DIAPH1, MET, MYO7A and SLC26A4; axes show counts and within-gene proportions. Figure 5 ranks conflict counts and proportions separately, with ties ordered alphabetically. The comparisons are descriptive; no sampling error bars or significance tests are implied.

## Study-design figure

Figure 1 shows the three nested phenotype selections and the separate 164-gene ClinGen comparator. Counts remain 18,384, 25,661, 473,704 and 162,683 distinct VariationIDs. The gene-reference size and nesting are taken from `results/clingen_grade_sensitivity.tsv` and `results/definition_overlaps.tsv`. Summaries are calculated within each selection set. Figure-specific generation hashes are recorded in `figure_provenance.json`.
