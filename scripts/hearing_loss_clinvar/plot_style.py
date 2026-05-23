from __future__ import annotations


def apply_plos_one_style() -> None:
    """
    Apply a consistent, journal-friendly matplotlib style tuned for PLOS ONE:
    - Sans-serif fonts (Arial/Helvetica fallback)
    - Legible sizes when downscaled
    - Clean axes, light grids off by default
    - Vector PDF text as editable (fonttype 42)
    """

    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.titlepad": 10,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.fontsize": 9,
            "legend.frameon": False,
            "lines.linewidth": 1.2,
            "lines.markersize": 5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
