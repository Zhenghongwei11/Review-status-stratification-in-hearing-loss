#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        return list(r)


def to_int(x: str) -> int:
    try:
        return int(x)
    except Exception:
        return 0


def to_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate publication-style figures for the hearing-loss ClinVar route.")
    ap.add_argument("--outdir", default="plots/hearing_loss_clinvar")
    ap.add_argument("--utility-tier", default="results/hearing_loss_clinvar/utility_overall_by_tier.tsv")
    ap.add_argument("--utility-gene", default="results/hearing_loss_clinvar/utility_gene_all.tsv")
    ap.add_argument("--annotate-n", type=int, default=12)
    ap.add_argument("--sensitivity-summary", default="results/hearing_loss_clinvar_sensitivity/summary.tsv")
    args = ap.parse_args()

    # Local import so the repo can be used without plotting deps until needed.
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from plot_style import apply_plos_one_style

    apply_plos_one_style()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------
    # Figure 1: analysis overview schematic (no raw data)
    # -----------------------------------------
    # Research-only overview (avoid "workflow/triage/operational" phrasing).
    fig, ax = plt.subplots(figsize=(7.0, 2.45))
    ax.axis("off")

    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Compact, journal-friendly boxes (readable in 1–2 column widths).
    boxes = [
        ("Data source", "ClinVar variant summary export"),
        ("Phenotype definition", "Hearing-loss text definition\n+ sensitivity run"),
        ("Variant aggregation", "Unique VariationID\nReview-status tier proxy (0–4)\nClinical-significance categories"),
        ("Summary metrics", "Uncertainty burden by tier\nGene-level conflict patterns"),
    ]

    # Box geometry in axes coordinates.
    y0 = 0.20
    h = 0.62
    w = 0.215
    xs = [0.03, 0.275, 0.52, 0.765]

    face = "#eaf2fb"  # very light blue
    edge = "#2f3a4a"

    for i, (x, (hdr, body)) in enumerate(zip(xs, boxes), start=1):
        box = FancyBboxPatch(
            (x, y0),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=1.0,
            edgecolor=edge,
            facecolor=face,
            transform=ax.transAxes,
        )
        ax.add_patch(box)

        # Header
        ax.text(
            x + w / 2,
            y0 + h * 0.78,
            hdr,
            ha="center",
            va="center",
            fontsize=9.2,
            fontweight="bold",
            color="#111111",
            transform=ax.transAxes,
        )
        # Body
        ax.text(
            x + w / 2,
            y0 + h * 0.40,
            body,
            ha="center",
            va="center",
            fontsize=8.1,
            color="#111111",
            linespacing=1.15,
            transform=ax.transAxes,
        )

        # Step number badge (subtle, top-left).
        ax.text(
            x + 0.018,
            y0 + h - 0.028,
            f"{i}",
            ha="left",
            va="top",
            fontsize=8.0,
            fontweight="bold",
            color=edge,
            transform=ax.transAxes,
        )

    # Arrows between boxes.
    for x in xs[:-1]:
        start = (x + w, y0 + h / 2)
        end = (x + w + (xs[1] - xs[0] - w) * 0.72, y0 + h / 2)
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=10,
                linewidth=1.0,
                color=edge,
                transform=ax.transAxes,
            )
        )

    fig.tight_layout(pad=0.20)
    fig.savefig(outdir / "fig1_workflow_schematic.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig1_workflow_schematic.pdf", bbox_inches="tight")
    plt.close(fig)

    # -------------------------
    # Figure 2: tier yield bars
    # -------------------------
    tier_rows = read_tsv(Path(args.utility_tier))
    # sort tiers descending numeric
    tier_rows.sort(key=lambda r: to_int(r.get("review_tier_proxy", r.get("review_tier_proxy_max", "0"))), reverse=True)

    tiers = [r.get("review_tier_proxy", r.get("review_tier_proxy_max", "")) for r in tier_rows]
    n_plp = [to_int(r.get("n::pathogenic_or_likely_pathogenic", "0")) for r in tier_rows]
    n_vus = [to_int(r.get("n::vus", "0")) for r in tier_rows]
    n_conf = [to_int(r.get("n::conflicting", "0")) for r in tier_rows]
    n_benign = [to_int(r.get("n::benign_or_likely_benign", "0")) for r in tier_rows]
    n_other = [to_int(r.get("n::other", "0")) + to_int(r.get("n::unknown", "0")) for r in tier_rows]

    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    x = list(range(len(tiers)))
    bottom = [0] * len(tiers)

    def stack(values, label, color):
        nonlocal bottom
        ax.bar(x, values, bottom=bottom, label=label, color=color, edgecolor="white", linewidth=0.4)
        bottom = [b + v for b, v in zip(bottom, values)]

    stack(n_plp, "P/LP", "#1b9e77")
    stack(n_vus, "VUS", "#7570b3")
    stack(n_conf, "Conflicting", "#d95f02")
    stack(n_benign, "B/LB", "#66a61e")
    stack(n_other, "Other/Unknown", "#999999")

    ax.set_xticks(x)
    ax.set_xticklabels([f"Tier {t}" for t in tiers])
    ax.set_ylabel("Unique variants (VariationID)")
    ax.set_title("ClinVar hearing-loss variants by review tier: yield vs uncertainty")
    ax.legend(frameon=False, ncol=2, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "fig2_tier_yield_stackedbar.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig2_tier_yield_stackedbar.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---------------------------------------------
    # Figure 3: yield vs uncertainty scatter
    # ---------------------------------------------
    gene_rows = read_tsv(Path(args.utility_gene))
    for r in gene_rows:
        r["n_unique_variants"] = to_int(r.get("n_unique_variants", "0"))
        r["uncertainty_burden"] = to_int(r.get("uncertainty_burden", "0"))
        r["highconf_rate"] = to_float(r.get("highconf_rate", "0"))
        r["conflict_rate"] = to_float(r.get("conflict_rate", "0"))

    # Keep genes with at least a minimal evidence footprint for plot readability.
    gene_rows_plot = [r for r in gene_rows if r["n_unique_variants"] >= 10]

    xs = [r["highconf_rate"] for r in gene_rows_plot]
    ys = [r["uncertainty_burden"] for r in gene_rows_plot]
    sizes = [max(12, min(220, r["n_unique_variants"] ** 0.7)) for r in gene_rows_plot]
    colors = [r["conflict_rate"] for r in gene_rows_plot]

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="viridis", alpha=0.85, edgecolor="white", linewidth=0.25)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Conflict Rate")

    ax.set_xlabel("Higher-review P/LP rate (tier≥2)")
    ax.set_ylabel("Uncertainty Burden (unique variants)")
    ax.set_title("Gene-level yield vs uncertainty (ClinVar hearing-loss evidence map)")
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate extremes:
    # - highest P/LP yield under higher-review threshold
    # - top hotspots (uncertainty burden)
    gene_rows_plot_sorted_action = sorted(
        gene_rows_plot, key=lambda r: (-to_int(r.get("n_highconf_plp_tier_ge_2", "0")), -r["n_unique_variants"], r["gene_symbol"])
    )
    gene_rows_plot_sorted_hot = sorted(gene_rows_plot, key=lambda r: (-r["uncertainty_burden"], -r["n_unique_variants"], r["gene_symbol"]))

    to_annotate = {r["gene_symbol"] for r in gene_rows_plot_sorted_action[: args.annotate_n // 2]}
    to_annotate |= {r["gene_symbol"] for r in gene_rows_plot_sorted_hot[: args.annotate_n // 2]}

    # Simple, deterministic annotation offsets by quadrant to reduce collisions.
    for r in gene_rows_plot:
        if r["gene_symbol"] not in to_annotate:
            continue
        dx = 5 if r["highconf_rate"] < 0.25 else -5
        dy = 4 if r["uncertainty_burden"] < 150 else -4
        ha = "left" if dx > 0 else "right"
        va = "bottom" if dy > 0 else "top"
        ax.annotate(
            r["gene_symbol"],
            (r["highconf_rate"], r["uncertainty_burden"]),
            textcoords="offset points",
            xytext=(dx, dy),
            fontsize=8,
            ha=ha,
            va=va,
        )

    fig.tight_layout()
    fig.savefig(outdir / "fig3_actionability_vs_uncertainty.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig3_actionability_vs_uncertainty.pdf", bbox_inches="tight")
    plt.close(fig)

    # -------------------------------------------------
    # Figure 4: conflict hotspot lollipop (top 25 genes)
    # -------------------------------------------------
    hot_path = Path("results/hearing_loss_clinvar/utility_gene_conflict_hotspots_top50.tsv")
    if hot_path.exists():
        hot = read_tsv(hot_path)
        # Keep top 25 for readability
        hot = hot[:25]
        genes = [r["gene_symbol"] for r in hot][::-1]
        n_conf = [to_int(r.get("n_conflicting", "0")) for r in hot][::-1]
        n_unc = [to_int(r.get("uncertainty_burden", "0")) for r in hot][::-1]

        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        y = list(range(len(genes)))
        ax.hlines(y, [0] * len(y), n_conf, color="#d95f02", lw=1.2, alpha=0.9, label="Conflicting variants")
        ax.plot(n_conf, y, "o", color="#d95f02", markersize=5)
        # Overlay uncertainty burden as light points (context)
        ax.plot(n_unc, y, "o", color="#7570b3", markersize=4, alpha=0.65, label="Uncertainty Burden (VUS + conflicting)")

        ax.set_yticks(y)
        ax.set_yticklabels(genes)
        ax.set_xlabel("Count (unique variants)")
        ax.set_title("Conflict hotspot genes in ClinVar hearing-loss records (top 25)")
        ax.legend(frameon=False, fontsize=9, loc="lower right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(outdir / "fig4_conflict_hotspots_lollipop.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / "fig4_conflict_hotspots_lollipop.pdf", bbox_inches="tight")
        plt.close(fig)

    # --------------------------------------------
    # Figure 5 / Supplementary: phenotype sensitivity
    # --------------------------------------------
    sens_path = Path(args.sensitivity_summary)
    if sens_path.exists():
        sens = read_tsv(sens_path)
        # Expect two rows: hearing_loss_only and expanded_synonyms
        # Build a compact multi-panel "slope" summary to avoid duplicating bar-chart style used elsewhere.
        def get(label: str, key: str) -> int:
            for r in sens:
                if r.get("run_label") == label:
                    return to_int(r.get(key, "0"))
            return 0

        labels = ["hearing_loss_only", "expanded_synonyms"]
        primary = labels[0]
        sensitivity = labels[1]

        uniq_primary = get(primary, "hearing_loss_variation_ids_unique")
        uniq_sens = get(sensitivity, "hearing_loss_variation_ids_unique")
        high_primary = get(primary, "highconf_plp_tier_ge_2_variation_ids_unique")
        high_sens = get(sensitivity, "highconf_plp_tier_ge_2_variation_ids_unique")
        genes_primary = get(primary, "unique_genes_hearing_loss")
        genes_sens = get(sensitivity, "unique_genes_hearing_loss")

        metrics = [
            ("Unique HL variants", uniq_primary, uniq_sens, "Count (unique VariationID)"),
            ("Higher-review P/LP", high_primary, high_sens, "Count (unique VariationID)"),
            ("Unique HL genes", genes_primary, genes_sens, "Count (unique genes)"),
        ]

        fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.6))
        for ax, (title, v1, v2, ylab) in zip(axes, metrics):
            ax.plot([0, 1], [v1, v2], color="#444444", lw=1.2, zorder=1)
            ax.scatter([0], [v1], color="#4c72b0", s=28, zorder=2)
            ax.scatter([1], [v2], color="#dd8452", s=28, zorder=2)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Primary", "Sensitivity"], rotation=0)
            ax.set_title(title, fontsize=10)
            ax.set_ylabel(ylab)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            # Annotate values
            ax.annotate(f"{v1:,}", (0, v1), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
            ax.annotate(f"{v2:,}", (1, v2), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)

        fig.suptitle("Sensitivity to phenotype text definition (ClinVar PhenotypeList)", y=1.02, fontsize=11)
        fig.tight_layout()
        # Main-figure filename
        fig.savefig(outdir / "fig5_phenotype_sensitivity.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / "fig5_phenotype_sensitivity.pdf", bbox_inches="tight")
        # Keep a stable supplementary filename as well.
        fig.savefig(outdir / "figS1_sensitivity_summary.png", dpi=300, bbox_inches="tight")
        fig.savefig(outdir / "figS1_sensitivity_summary.pdf", bbox_inches="tight")
        plt.close(fig)

    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
