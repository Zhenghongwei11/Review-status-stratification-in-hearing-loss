# Figure provenance

Renderer: `scripts/hearing_loss_clinvar/plot_figures.py`. The renderer resolves input and output paths from its own location.

Figures use white backgrounds, an accessible discrete palette, embedded TrueType PDF fonts, and 300 dpi PNG exports. Category labels, classification hatches, labelled points and marker shapes supplement colour for grayscale interpretation. The dataset is a complete selected snapshot; no sampling error bars or significance tests are implied.

Each `plots/figN_source_data.tsv` contains the exact values displayed. `figure_provenance.json` records checksums of every input table read by this run.

- Global variant counts use unique VariationIDs within each definition.
- Gene plots include genes with at least 10 annotated unique variants. A variant annotated to several genes contributes once to each; totals across genes are not additive.
- The review-status-defined P/LP subset combines the P/LP classification bucket with at least 2 official ClinVar review stars. It is a descriptive subset and does not represent clinical reportability or diagnostic yield.
- Figure 4 labels DIAPH1, MET and the two genes with the largest counts in that review-status-defined P/LP subset. All displayed proportions are recomputed from integer numerators and denominators, not rounded utility-table rates.
- Figure 5 reports independent rankings of the 12 highest conflict counts and the 12 highest conflict proportions among genes with at least 10 annotated VariationIDs; the panels may include different genes. Single-gene sensitivity results are reported in S5 Table.
- Figure 2 compares different selection universes; curated-gene annotation includes records for other disorders. HPO and gene curation versions must be reported separately from the fixed ClinVar input date.

## Input checksums

- `data/references/references_manifest.json`: `31c8e7cec3b59f137d843c62c1face903fa6c09991ff568f1b7c52d3bf46e093`
- `results/primary/utility_gene_all.tsv`: `c559233a41ca454561649677630fdaedb268beb866078aea1e2c3ef101097bf7`
- `results/run_metadata.json`: `f9b34e2387d71034afef6b296a3f51bbb5f1ff0be88f3cbb1a5ec233171ba459`

## External reference versions

- hp.obo: version `v2026-09-01`; retrieved `2026-10-13T08:56:38.105503+00:00`; SHA256 `93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22`; source https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/hp.obo
- phenotype.hpoa: version `v2026-09-01`; retrieved `2026-10-13T08:56:38.171859+00:00`; SHA256 `e89aa39c8f97bf5a52c5d160f9250a603681d630ade8ec2679d84ac5aece1f72`; source https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype.hpoa
- hpo_license.md: version `v2026-09-01`; retrieved `2026-10-13T08:56:41.923564+00:00`; SHA256 `534f0ba30b74d13590247ddb5e24d749f5733b41cf19a730010b55edf4266336`; source https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/v2026-09-01/LICENSE.md
- hpo_legacy_license.html: version `retrieval-time official legacy license page`; retrieved `2026-10-13T08:56:51.158477+00:00`; SHA256 `d63799ff8d381891114b22e78c8164f198a9642f9ed0ca8b53bc73238c318b39`; source https://human-phenotype-ontology.github.io/license.html
- clingen_gene_validity.csv: version `live export; retrieval date and checksum identify snapshot`; retrieved `2026-10-13T08:59:52.798247+00:00`; SHA256 `e76ad9ff65af9515d0f909b351f440808d584c67842dc914026d3a8d0b20db10`; source https://search.clinicalgenome.org/kb/gene-validity/download
- clingen_terms.html: version `retrieval-time terms page`; retrieved `2026-10-13T09:02:12.955026+00:00`; SHA256 `bbeed8bdb6ba5029ee803b17b2e9a46cc92c9787c2cc5420a9c93477c50dc6ec`; source https://clinicalgenome.org/about/terms-of-use/
