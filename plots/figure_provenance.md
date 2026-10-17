# Figure provenance

Renderer: `scripts/hearing_loss_clinvar/plot_revision.py`. The renderer resolves paths from its own location and never writes to submitted figure paths.

Figures use white backgrounds, an accessible discrete palette, embedded TrueType PDF fonts, and 300 dpi PNG exports. Category labels, classification hatches, labelled points and marker shapes supplement colour for grayscale interpretation. The dataset is a complete selected snapshot; no sampling error bars or significance tests are implied.

Each `plots/hearing_loss_clinvar_revision1/figN_source_data.tsv` contains the exact values displayed. `figure_provenance.json` records checksums of every input table read by this run.

- Global variant counts use unique VariationIDs within each definition.
- Gene plots include genes with at least 10 annotated unique variants. A variant annotated to several genes contributes once to each; totals across genes are not additive.
- The review-status-defined P/LP subset combines the P/LP classification bucket with at least 2 official ClinVar review stars. It is a descriptive subset and does not represent clinical reportability or diagnostic yield.
- Figure 4 labels DIAPH1, MET and the two genes with the largest counts in that review-status-defined P/LP subset. All displayed proportions are recomputed from integer numerators and denominators, not rounded utility-table rates.
- Figure 5 reports independent rankings of the 12 highest conflict counts and the 12 highest conflict proportions among genes with at least 10 annotated VariationIDs; the panels may include different genes. Single-gene sensitivity results are reported in S5 Table.
- Figure 2 compares different selection universes; curated-gene annotation includes records for other disorders. HPO and gene curation versions must be reported separately from the fixed ClinVar input date.

## Input checksums

- `data/revision1/references/references_manifest.json`: `dc4511ff0c2b8a08c0a72d94ca834accaa339c73e2a772d6c5f03def2f72ea32`
- `results/hearing_loss_clinvar_revision1/coverage_summary.tsv`: `ffbfbc80c88765bb77c78b86a3a46da37737227be35004483cca29b446e6d5d1`
- `results/hearing_loss_clinvar_revision1/curated_gene/variant_aggregate.tsv`: `56ed5f57ab5f762049eb5494700d411004468b3c6bb14f568760c755a9766703`
- `results/hearing_loss_clinvar_revision1/disease_supported/variant_aggregate.tsv`: `f1d8e6c74de1047e3fade9df2734ba9ee21175456cf8087aa9447c6ad971b5ea`
- `results/hearing_loss_clinvar_revision1/expanded/variant_aggregate.tsv`: `e48e0c8628e520a94c5ff5c7660c3ac6ee10abefeea4c9d8ffa4deb6a7435ab7`
- `results/hearing_loss_clinvar_revision1/primary/utility_gene_all.tsv`: `c559233a41ca454561649677630fdaedb268beb866078aea1e2c3ef101097bf7`
- `results/hearing_loss_clinvar_revision1/primary/variant_aggregate.tsv`: `3f02c6c0377760926fefbb12f8d0cc6f14116e54201da4a209b3c383bb0446fe`
- `results/hearing_loss_clinvar_revision1/run_metadata.json`: `beb558fc2b4e53f9fc184b8feea143afa5f7504e89f0abf61a6b6bd0d84bb3ba`

## External reference versions

- hp.obo: version `v2026-09-01`; retrieved `2026-10-13T08:56:38.123011+00:00`; SHA256 `93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22`; source https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/hp.obo
- phenotype.hpoa: version `v2026-09-01`; retrieved `2026-10-13T08:56:38.190350+00:00`; SHA256 `e89aa39c8f97bf5a52c5d160f9250a603681d630ade8ec2679d84ac5aece1f72`; source https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype.hpoa
- hpo_license.md: version `v2026-09-01`; retrieved `2026-10-13T08:56:41.931309+00:00`; SHA256 `534f0ba30b74d13590247ddb5e24d749f5733b41cf19a730010b55edf4266336`; source https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/v2026-09-01/LICENSE.md
- hpo_legacy_license.html: version `retrieval-time official legacy license page`; retrieved `2026-10-13T08:56:51.163641+00:00`; SHA256 `d63799ff8d381891114b22e78c8164f198a9642f9ed0ca8b53bc73238c318b39`; source https://human-phenotype-ontology.github.io/license.html
- clingen_gene_validity.csv: version `live export; retrieval date and checksum identify snapshot`; retrieved `2026-10-13T08:59:52.798247+00:00`; SHA256 `e76ad9ff65af9515d0f909b351f440808d584c67842dc914026d3a8d0b20db10`; source https://search.clinicalgenome.org/kb/gene-validity/download
- clingen_terms.html: version `retrieval-time terms page`; retrieved `2026-10-13T09:02:12.955026+00:00`; SHA256 `bbeed8bdb6ba5029ee803b17b2e9a46cc92c9787c2cc5420a9c93477c50dc6ec`; source https://clinicalgenome.org/about/terms-of-use/
