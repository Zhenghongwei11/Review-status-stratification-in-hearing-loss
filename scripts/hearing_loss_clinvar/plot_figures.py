#!/usr/bin/env python3
"""Render five figures from fixed-input ClinVar analysis results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import FuncFormatter, MaxNLocator, PercentFormatter
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
COLORS = ["#0072B2", "#E69F00", "#CC79A7", "#009E73", "#888888", "#444444"]
HATCHES = ["", "///", "...", "\\\\", "xx", "++"]
BUCKETS = ["pathogenic_or_likely_pathogenic", "vus", "conflicting", "benign_or_likely_benign", "other", "unknown"]
LABELS = ["P/LP", "VUS", "Conflicting", "B/LB", "Other", "Unknown"]
DEFINITIONS = [("primary", "Hearing-loss phrase"), ("expanded", "Expanded text"),
               ("disease_supported", "Text + HPO-linked disease"), ("curated_gene", "Curated-gene comparator")]


def read_tsv(path):
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def file_sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No source rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def integer(row, key):
    return int(row.get(key) or 0)


def style():
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "legend.fontsize": 8, "pdf.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "axes.spines.top": False, "axes.spines.right": False})


def panel(ax, letter, title):
    ax.set_title(title, loc="left", pad=12)
    ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=12)


def count_format(ax, axis="y"):
    target = ax.yaxis if axis == "y" else ax.xaxis
    target.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))


class Renderer:
    def __init__(self, results, outdir, docs):
        self.results, self.outdir, self.docs = results, outdir, docs
        outdir.mkdir(parents=True, exist_ok=True); docs.mkdir(parents=True, exist_ok=True)
        self.inputs = {}; self.captions = {}; self.output_files = []

    def load(self, relative):
        return list(self.iterload(relative))

    def iterload(self, relative):
        path = self.results / relative
        self.inputs[str(path)] = file_sha(path)
        with path.open(encoding="utf-8") as stream:
            yield from csv.DictReader(stream, delimiter="\t")

    def save(self, fig, number, rows, caption):
        stem = self.outdir / f"fig{number}"
        for suffix in ["png", "pdf"]:
            fig.savefig(stem.with_suffix("." + suffix), dpi=300, facecolor="white", bbox_inches="tight")
            if suffix == "png":
                # Export opaque RGB; journals may reject a PNG alpha channel even
                # when every alpha value is 255. This changes no plotted values.
                with Image.open(stem.with_suffix(".png")) as png:
                    png.convert("RGB").save(stem.with_suffix(".png"), dpi=(300,300))
            self.output_files.append(str(stem.with_suffix("." + suffix)))
        with Image.open(stem.with_suffix(".png")) as png:
            png.convert("RGB").save(self.outdir / f"Fig{number}.tif", compression="tiff_lzw", dpi=(300, 300))
        self.output_files.append(str(self.outdir / f"Fig{number}.tif"))
        write_tsv(self.outdir / f"fig{number}_source_data.tsv", rows)
        plt.close(fig)
        self.captions[number] = caption

    def plot_study_design_and_record_selection(self):
        coverage = {row["definition"]: row for row in self.load("coverage_summary.tsv")}
        definitions = ["primary", "expanded", "disease_supported", "curated_gene"]
        if any(name not in coverage for name in definitions):
            raise ValueError("Figure 1 requires all four selection definitions in coverage_summary.tsv")
        counts = {name: integer(coverage[name], "n_variation_ids") for name in definitions}
        if any(count <= 0 for count in counts.values()):
            raise ValueError("Figure 1 requires a positive VariationID count for each definition")
        fig, ax = plt.subplots(figsize=(7.4, 4.4)); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        boxes = [
            (0.22, 0.81, 0.56, 0.16, "ClinVar variation records", "Accessed 2026-05-23"),
            (0.02, 0.47, 0.45, 0.27, "Phenotype-based selection", "Hearing-loss phrase: n = " + f"{counts['primary']:,}" + " distinct VariationIDs" +
             "\nExpanded text: n = " + f"{counts['expanded']:,}" + " distinct VariationIDs" +
             "\nText + HPO diseases: n = " + f"{counts['disease_supported']:,}" + " distinct VariationIDs"),
            (0.53, 0.47, 0.45, 0.27, "Gene-based selection", "Curated-gene comparator:" +
             "\nn = " + f"{counts['curated_gene']:,}" + " distinct VariationIDs" +
             "\nRecords across annotated conditions"),
            (0.12, 0.23, 0.76, 0.17, "One aggregate record per VariationID", "Aggregate classifications and review status\nGene memberships and disease annotations"),
            (0.12, 0.015, 0.76, 0.13, "Descriptive comparisons", "Coverage / review status / gene membership counts and proportions")]
        for x, y, w, h, title, detail in boxes:
            ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.009",linewidth=0.9,edgecolor="#555555",facecolor="#F6F8FA"))
            ax.text(x+w/2, y+h*0.73, title, ha="center", va="center", fontweight="bold", fontsize=9)
            ax.text(x+w/2, y+h*0.31, detail, ha="center", va="center", fontsize=8, linespacing=1.4)
        for xy, xytext in [((.245,.753),(.40,.80)),((.755,.753),(.60,.80)),((.35,.412),(.245,.46)),((.65,.412),(.755,.46)),((.5,.158),(.5,.22))]:
            ax.annotate("",xy=xy,xytext=xytext,arrowprops={"arrowstyle":"->","color":"#444444","lw":1.1})
        labels = {"primary": "Hearing-loss phrase", "expanded": "Expanded text", "disease_supported": "Text plus HPO-linked diseases", "curated_gene": "Curated-gene annotation comparator"}
        rows = [{"definition": name, "selection_rule": labels[name], "unique_variation_ids": counts[name]} for name in definitions]
        caption = ("Study design and record selection. Three phenotype-based rules and a gene-based comparator were applied to the same ClinVar export. "
            "Each branch reports its number of distinct VariationIDs; annotated gene memberships are summarized separately. "
            "The curated-gene comparator includes records across annotated conditions.")
        self.save(fig,1,rows,caption)

    def plot_review_status_distribution(self):
        variants = self.load("primary/variant_aggregate.tsv")
        rows = []; stars = list(range(5))
        for star in stars:
            values = [v for v in variants if integer(v,"review_tier_proxy_max") == star]
            for bucket, label in zip(BUCKETS,LABELS):
                rows.append({"review_stars":star,"classification_bucket":bucket,"classification_label":label,
                    "n_variants":sum(v["bucket_variant_level"]==bucket for v in values),"n_in_star_group":len(values)})
        unknown = sum(integer(v,"review_tier_proxy_max") < 0 for v in variants)
        if unknown:
            raise ValueError("Unknown review statuses need an explicit Figure 2 category before publication")
        fig, axes = plt.subplots(1,2,figsize=(7.5,3.7),gridspec_kw={"width_ratios":[1.35,1]})
        bottom = np.zeros(5)
        for bucket, label, color, hatch in zip(BUCKETS,LABELS,COLORS,HATCHES):
            vals = np.array([r["n_variants"] for r in rows if r["classification_bucket"]==bucket])
            if vals.sum():
                axes[0].bar(stars,vals,bottom=bottom,label=label,color=color,width=.65,edgecolor="white",linewidth=.4,hatch=hatch)
            bottom += vals
        for x, total in zip(stars,bottom):
            axes[0].text(x,total+max(bottom)*.025,f"{int(total):,}",ha="center",fontsize=8)
        axes[0].set_ylim(0,max(bottom)*1.18); axes[0].set_xticks(stars,["0","1","2","3","4"])
        axes[0].set_xlabel("Official ClinVar review stars"); axes[0].set_ylabel("Unique variants"); count_format(axes[0])
        panel(axes[0],"A",f"Hearing-loss phrase (n = {len(variants):,})")
        expert = [r for r in rows if r["review_stars"] == 3 and r["n_variants"]]
        y = np.arange(len(expert))
        axes[1].barh(y,[r["n_variants"] for r in expert],color=[COLORS[BUCKETS.index(r["classification_bucket"])] for r in expert],height=.65)
        axes[1].set_yticks(y,[r["classification_label"] for r in expert]); axes[1].invert_yaxis()
        for yi,r in zip(y,expert): axes[1].text(r["n_variants"]+1.5,yi,str(r["n_variants"]),va="center",fontsize=8)
        axes[1].set_xlim(0,max(r["n_variants"] for r in expert)*1.22); axes[1].set_xlabel("Unique variants")
        axes[1].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=4))
        n_expert = sum(r["n_variants"] for r in expert)
        panel(axes[1],"B",f"Expert panel: 3 stars (n = {n_expert:,})")
        handles, labels=axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc="lower center",ncol=5,frameon=False,bbox_to_anchor=(.5,-.01))
        fig.subplots_adjust(left=.09,right=.97,top=.87,bottom=.22,wspace=.40)
        self.save(fig,3,rows,f"Review-status distribution under the phrase rule. (A) Classification categories by official review stars for {len(variants):,} VariationIDs. (B) Composition of the {n_expert:,} expert-panel records. The Other category represents one drug-response record (n = 1). No VariationIDs had 4-star review status (practice guidelines; n = 0). P/LP, pathogenic/likely pathogenic; VUS, variant of uncertain significance; B/LB, benign/likely benign.")

    def plot_gene_level_classification_profiles(self):
        genes = self.load("primary/utility_gene_all.tsv")
        selected = [g for g in genes if integer(g,"n_unique_variants") >= 10]
        names = {"DIAPH1", "MET"}
        leaders = sorted(selected,key=lambda g:(-integer(g,"n_plp_review_stars_ge_2"),g["gene_symbol"]))[:2]
        names.update(g["gene_symbol"] for g in leaders)
        leader_positions = [(.67,.58),(.57,.86)]
        count_positions = {"DIAPH1":(.03,.91),"MET":(.03,.57)}
        proportion_positions = {"DIAPH1":(.05,.64),"MET":(.05,.96)}
        for g,pos in zip(leaders,leader_positions):
            count_positions[g["gene_symbol"]]=pos
            proportion_positions[g["gene_symbol"]]=pos
        fig, axes=plt.subplots(1,2,figsize=(7.5,4.1))
        rows=[]
        for g in selected:
            n=integer(g,"n_unique_variants"); high=integer(g,"n_plp_review_stars_ge_2"); unc=integer(g,"uncertainty_burden")
            conf=integer(g,"n_conflicting")
            row={"gene_symbol":g["gene_symbol"],"n_unique_variants":n,"n_plp_review_stars_ge_2":high,"n_vus_plus_conflicting":unc,
                 "n_conflicting":conf,"plp_review_stars_ge_2_proportion":high/n,"vus_plus_conflicting_proportion":unc/n,"conflicting_proportion":conf/n,"annotated":g["gene_symbol"] in names}
            rows.append(row)
        for ax,proportion in zip(axes,[False,True]):
            xv=[r["plp_review_stars_ge_2_proportion"] if proportion else r["n_plp_review_stars_ge_2"] for r in rows]
            yv=[r["vus_plus_conflicting_proportion"] if proportion else r["n_vus_plus_conflicting"] for r in rows]
            ax.scatter(xv,yv,s=16,color="#A6ADB4",alpha=.72,linewidths=.2,edgecolors="white",zorder=2)
            if proportion:
                ax.set_xlim(-.02,1.04); ax.set_ylim(-.02,1.10)
                ax.xaxis.set_major_formatter(PercentFormatter(1)); ax.yaxis.set_major_formatter(PercentFormatter(1))
                ax.set_xlabel("P/LP with ≥2 review stars / gene variants"); ax.set_ylabel("(VUS + conflicting) / gene variants")
            else:
                ax.set_xscale("symlog",linthresh=1); ax.set_yscale("symlog",linthresh=1)
                ax.set_xlim(-.1,max(xv)*1.5); ax.set_ylim(-.1,max(yv)*1.9)
                ax.set_xticks([0,1,10,100]); ax.set_yticks([0,1,10,100,1000]); count_format(ax); count_format(ax,"x")
                ax.set_xlabel("P/LP variants with ≥2 review stars"); ax.set_ylabel("VUS + conflicting variants")
            for r in rows:
                if not r["annotated"]: continue
                target=(r["plp_review_stars_ge_2_proportion"],r["vus_plus_conflicting_proportion"]) if proportion else (r["n_plp_review_stars_ge_2"],r["n_vus_plus_conflicting"])
                color="#D55E00" if r["gene_symbol"] in {"DIAPH1","MET"} else "#0072B2"
                ax.scatter(*target,s=30,color=color,edgecolor="white",linewidth=.5,zorder=3)
                positions = proportion_positions if proportion else count_positions
                ax.annotate(r["gene_symbol"],xy=target,xytext=positions[r["gene_symbol"]],textcoords="axes fraction",
                    fontsize=8,fontstyle="italic",color=color,ha="left",va="center",
                    bbox={"boxstyle":"round,pad=.2","facecolor":"white","edgecolor":"none","alpha":.96},
                    arrowprops={"arrowstyle":"->","lw":.7,"color":color,"connectionstyle":"arc3,rad=.03"},zorder=4)
        panel(axes[0],"A","Gene membership counts")
        panel(axes[1],"B","Within-gene proportions")
        fig.subplots_adjust(left=.09,right=.98,bottom=.18,top=.88,wspace=.40)
        self.save(fig,4,rows,f"Gene-level classification profiles under the phrase rule. Points represent {len(selected)} genes with at least ten distinct VariationIDs. Orange highlights DIAPH1 and MET; blue highlights MYO7A and SLC26A4; grey marks the other {len(selected) - 4} genes. (A) Counts of P/LP records with at least two review stars versus VUS-plus-conflicting counts; values from 0 to 1 are shown linearly and values above 1 logarithmically on symmetric-log axes. (B) Corresponding within-gene proportions. Each VariationID contributes once to each annotated gene, so gene counts are nonadditive. P/LP, pathogenic/likely pathogenic; VUS, variant of uncertain significance.")

    def plot_conflict_gene_rankings(self):
        genes = self.load("primary/utility_gene_all.tsv")
        eligible = [g for g in genes if integer(g,"n_unique_variants") >= 10]
        by_count = sorted(eligible,key=lambda g:(-integer(g,"n_conflicting"),g["gene_symbol"]))[:12]
        by_proportion = sorted(eligible,key=lambda g:(-integer(g,"n_conflicting")/integer(g,"n_unique_variants"),g["gene_symbol"]))[:12]
        rows=[]
        ranked = {"A":by_count,"B":by_proportion}
        for panel_label, ranked_genes in ranked.items():
            for rank,g in enumerate(ranked_genes,1):
                n=integer(g,"n_unique_variants"); count=integer(g,"n_conflicting")
                rows.append({"panel":panel_label,"rank":rank,"gene":g["gene_symbol"],
                    "count":count,"n":n,"proportion":count/n})

        fig,axes=plt.subplots(1,2,figsize=(7.5,4.8),gridspec_kw={"width_ratios":[1,1.15]})
        for ax,panel_label in zip(axes,["A","B"]):
            panel_rows=[r for r in rows if r["panel"]==panel_label]
            y=np.arange(len(panel_rows))
            ax.set_yticks(y,[r["gene"] for r in panel_rows])
            ax.invert_yaxis()
            for label in ax.get_yticklabels():label.set_fontstyle("italic")
            if panel_label=="A":
                values=[r["count"] for r in panel_rows]
                ax.barh(y,values,height=.68,color="#CC79A7",edgecolor="white",linewidth=.4)
                xmax=max(values) if values else 1
                ax.set_xlim(0,xmax*1.18)
                for yi,r in zip(y,panel_rows):
                    ax.text(r["count"]+xmax*.025,yi,f'{r["count"]:,}',va="center",fontsize=8)
                ax.set_xlabel("Conflicting VariationIDs")
                count_format(ax,"x")
                panel(ax,"A","Highest conflict counts")
            else:
                values=[r["proportion"] for r in panel_rows]
                xmax=max(values) if values else 1
                ax.barh(y,values,height=.68,color="#CC79A7",edgecolor="white",linewidth=.4)
                ax.set_xlim(0,min(1.0,xmax*1.18 if xmax else 1.0))
                ax.xaxis.set_major_formatter(PercentFormatter(1))
                ax.set_xlabel("Conflicting / all gene VariationIDs")
                for yi,r in zip(y,panel_rows):
                    ax.text(1.03,yi,f'n={r["n"]:,}',transform=ax.get_yaxis_transform(),
                        ha="left",va="center",fontsize=8,color="#444444",clip_on=False)
                panel(ax,"B","Highest conflict proportions")
        fig.subplots_adjust(left=.11,right=.87,bottom=.17,top=.89,wspace=.58)
        self.save(fig,5,rows,"Gene rankings by conflict count and proportion. Genes had at least ten annotated VariationIDs. (A) Twelve largest conflict counts. (B) Twelve largest conflict proportions, with gene-level denominators shown as n. Panels are ranked independently; tied ranks are ordered alphabetically by gene symbol. Single-gene sensitivity results are in S5 Table.")

    def plot_selection_coverage_and_classification_composition(self):
        sets={}; rows=[]
        for mode,label in DEFINITIONS:
            variant_ids, genes, buckets, n_expert, n = set(), set(), Counter(), 0, 0
            for v in self.iterload(f"{mode}/variant_aggregate.tsv"):
                n += 1
                variant_ids.add(v["variation_id"])
                genes.update(g for g in v["gene_symbols"].split(";") if g)
                buckets[v["bucket_variant_level"]] += 1
                n_expert += integer(v,"review_tier_proxy_max") == 3
            if n != len(variant_ids):
                raise ValueError(f"Duplicated VariationIDs in {mode}")
            sets[mode]=variant_ids
            for bucket,blabel in zip(BUCKETS,LABELS):
                c=buckets[bucket]
                rows.append({"definition":mode,"definition_label":label,"n_unique_variants":n,"n_annotated_genes":len(genes),
                    "classification_bucket":bucket,"n_in_bucket":c,"bucket_proportion":c/n if n else 0,
                    "n_expert_panel":n_expert})
        fig,axes=plt.subplots(1,2,figsize=(7.5,4.4),gridspec_kw={"width_ratios":[1.1,1]})
        y=np.arange(4); totals=[len(sets[m]) for m,_ in DEFINITIONS]
        axes[0].barh(y,totals,color=["#D9D9D9","#BDBDBD","#969696","#636363"],height=.64)
        for yi,total in zip(y,totals):axes[0].text(total+max(totals)*.02,yi,f"{total:,}",va="center",fontsize=8)
        axes[0].set_yticks(y,[l for _,l in DEFINITIONS]); axes[0].invert_yaxis(); axes[0].set_xlim(0,max(totals)*1.24)
        axes[0].set_xlabel("Distinct VariationIDs");count_format(axes[0],"x");axes[0].xaxis.set_major_locator(MaxNLocator(nbins=3))
        left=np.zeros(4)
        for bucket,label,color,hatch in zip(BUCKETS,LABELS,COLORS,HATCHES):
            vals=np.array([next(r["bucket_proportion"] for r in rows if r["definition"]==m and r["classification_bucket"]==bucket) for m,_ in DEFINITIONS])
            if vals.sum(): axes[1].barh(y,vals,left=left,height=.64,color=color,label=label,edgecolor="white",linewidth=.4,hatch=hatch)
            left+=vals
        axes[1].set_yticks(y,[]);axes[1].invert_yaxis();axes[1].set_xlim(0,1);axes[1].xaxis.set_major_formatter(PercentFormatter(1))
        axes[1].set_xlabel("Within-definition classification proportion")
        panel(axes[0],"A","Record-set size");panel(axes[1],"B","Classification composition")
        handles,labels=axes[1].get_legend_handles_labels();fig.legend(handles,labels,ncol=5,loc="lower center",frameon=False,bbox_to_anchor=(.5,.015))
        fig.subplots_adjust(left=.29,right=.97,bottom=.20,top=.86,wspace=.33)
        self.save(fig,2,rows,"Record-set size and classification composition across selection rules. (A) Distinct VariationIDs selected by the phrase, expanded-text and HPO-linked rules and the gene comparator. (B) Classification proportions among all selected records in each set, regardless of review status. P/LP, pathogenic/likely pathogenic; VUS, variant of uncertain significance; B/LB, benign/likely benign.")

    def document(self):
        caption_path=self.docs/"figure_captions.md"
        prior=caption_path.read_text() if caption_path.exists() else ""
        # A partial figure run updates only those sections while retaining completed captions.
        for number in range(1,6):
            if number not in self.captions and f"## Figure {number}\n" in prior:
                self.captions[number]=prior.split(f"## Figure {number}\n",1)[1].split("\n## Figure ",1)[0].strip()
        caption_path.write_text("# Figure captions\n\n"+"\n\n".join(f"## Figure {n}\n\n{c}" for n,c in sorted(self.captions.items()))+"\n")
        metadata_path=self.results/"run_metadata.json"
        analysis_metadata=json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
        if metadata_path.exists():self.inputs[str(metadata_path)]=file_sha(metadata_path)
        references_path=ROOT/"data/references/references_manifest.json"
        references=json.loads(references_path.read_text()) if references_path.exists() else {}
        if references_path.exists():self.inputs[str(references_path)]=file_sha(references_path)
        # Public stable source URLs suffice; expiring signed redirect URLs are not copied.
        resources={name:{key:source.get(key) for key in ["requested_url","version","retrieved_at_utc","sha256"]}
            for name,source in references.get("sources",{}).items()}
        provenance={"renderer":str(Path(__file__).resolve().relative_to(ROOT)),
            "renderer_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"inputs_sha256":{str(Path(path).relative_to(ROOT)) if Path(path).is_relative_to(ROOT) else path: sha for path, sha in self.inputs.items()},
            "analysis_metadata":analysis_metadata,
            "reference_resources":resources,"reference_selection":references.get("selection",{}),
            "hpo_annotation_metadata":references.get("hpo_annotation_metadata",[]),
            "temporal_limitation":references.get("temporal_limitation",""),
            "outputs":[str(Path(path).relative_to(ROOT)) if Path(path).is_relative_to(ROOT) else path for path in self.output_files],"png_dpi":300,"png_color_mode":"RGB","pdf_fonttype":42,
            "classification_buckets":dict(zip(BUCKETS,LABELS)),"classification_hatches":dict(zip(BUCKETS,HATCHES)),"gene_minimum_variants":10,
            "plp_review_stars_ge_2":"P/LP and official ClinVar review stars >= 2 on the same selected source row","gene_memberships":"Each unique variant counted once per annotated gene; gene counts are not additive",
            "statistics":"Descriptive complete-snapshot counts and proportions; no inference, error bars, or significance tests"}
        (self.outdir/"figure_provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
        lines=["# Figure provenance", "", "Renderer: `scripts/hearing_loss_clinvar/plot_figures.py`. The renderer resolves input and output paths from its own location.", "",
            "Figures use white backgrounds, an accessible discrete palette, embedded TrueType PDF fonts, and 300 dpi PNG exports. Category labels, classification hatches, labelled points and marker shapes supplement colour for grayscale interpretation. The dataset is a complete selected snapshot; no sampling error bars or significance tests are implied.", "",
            "Each `plots/figN_source_data.tsv` contains the exact values displayed. `figure_provenance.json` records checksums of every input table read by this run.", "",
            "- Global variant counts use unique VariationIDs within each definition.",
            "- Gene plots include genes with at least 10 annotated unique variants. A variant annotated to several genes contributes once to each; totals across genes are not additive.",
            "- The review-status-defined P/LP subset combines the P/LP classification bucket with at least 2 official ClinVar review stars. It is a descriptive subset and does not represent clinical reportability or diagnostic yield.",
            "- Figure 4 labels DIAPH1, MET and the two genes with the largest counts in that review-status-defined P/LP subset. All displayed proportions are recomputed from integer numerators and denominators, not rounded utility-table rates.",
            "- Figure 5 reports independent rankings of the 12 highest conflict counts and the 12 highest conflict proportions among genes with at least 10 annotated VariationIDs; the panels may include different genes. Single-gene sensitivity results are reported in S5 Table.",
            "- Figure 2 compares different selection universes; curated-gene annotation includes records for other disorders. HPO and gene curation versions must be reported separately from the fixed ClinVar input date.", "", "## Input checksums", ""]
        lines.extend(f"- `{Path(path).relative_to(ROOT) if Path(path).is_relative_to(ROOT) else path}`: `{sha}`" for path,sha in sorted(self.inputs.items()))
        if resources:
            lines.extend(["", "## External reference versions", ""])
            lines.extend(f"- {name}: version `{info['version']}`; retrieved `{info['retrieved_at_utc']}`; SHA256 `{info['sha256']}`; source {info['requested_url']}" for name,info in resources.items())
        (self.docs/"figure_provenance.md").write_text("\n".join(lines)+"\n")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results",type=Path,default=ROOT/"results")
    ap.add_argument("--outdir",type=Path,default=ROOT/"plots")
    ap.add_argument("--docs",type=Path,default=ROOT/"plots")
    ap.add_argument("--figures",type=int,nargs="+",choices=range(1,6),default=list(range(1,6)))
    args=ap.parse_args()
    style()
    renderer=Renderer(args.results.resolve(),args.outdir.resolve(),args.docs.resolve())
    renderers = {
        1: renderer.plot_study_design_and_record_selection,
        2: renderer.plot_selection_coverage_and_classification_composition,
        3: renderer.plot_review_status_distribution,
        4: renderer.plot_gene_level_classification_profiles,
        5: renderer.plot_conflict_gene_rankings,
    }
    for number in args.figures: renderers[number]()
    renderer.document(); print(f"Rendered figures {args.figures} into {renderer.outdir}")


if __name__=="__main__":main()
