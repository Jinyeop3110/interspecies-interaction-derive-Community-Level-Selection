#!/usr/bin/env python3
"""Figure 6 panels: coalescence in natural sample-derived communities.

Generates:
  6b  similarity map per medium, with the parental-label-swapped counterpart
  6c  fraction of Dominance outcomes across the nutrient gradient

    python pipeline/04_figure_generation/fig6/generate_figure.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import (  # noqa: E402
    MEDIUM_COLORS, MM, SWAPPED_ALPHA, SWAPPED_COLOR,
    draw_similarity_boundaries, marginal_pdi_histogram, save,
    stacked_outcome_series, use_paper_style,
)

from coalescence import io, outcomes  # noqa: E402


def panel_b(table):
    """Similarity map for the natural sample-derived communities."""
    fig, axes = plt.subplots(2, 3, figsize=(180 * MM, 83 * MM),
                             gridspec_kw={"height_ratios": [60, 22.6]},
                             facecolor="w", edgecolor="k")
    fig.subplots_adjust(wspace=0.3, hspace=0.45)

    for column, medium in enumerate(io.MEDIA_ORDER):
        ax = axes[0, column]
        label = io.MEDIUM_LABELS[medium]
        group = table[table.Medium == medium]

        ax.scatter(group.x1, group.x2, s=25, color=MEDIUM_COLORS[label],
                   alpha=0.7, linewidths=0)
        ax.scatter(group.x2, group.x1, s=25, color=SWAPPED_COLOR,
                   alpha=SWAPPED_ALPHA, linewidths=0)

        draw_similarity_boundaries(ax)
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_title(label, fontsize=7, style="italic",
                     color=MEDIUM_COLORS[label])

        histogram = axes[1, column]
        marginal_pdi_histogram(histogram, group.PDI.to_numpy(dtype=float),
                               1 - group.PDI.to_numpy(dtype=float),
                               MEDIUM_COLORS[label])
        if column == 1:
            histogram.set_xlabel("Parental Dominance Index")

    axes[0, 0].set_ylabel("Similarity(C,B)")
    axes[0, 1].set_xlabel("Similarity(C,A)")
    axes[1, 0].set_ylabel("Density")
    return save(fig, "fig6b_similarity_map")


def panel_c(fractions):
    """Outcome fractions across media.

    The published panel is the full three-outcome stacked chart, in the same
    style as Fig. 4d -- not a Dominance-only series.
    """
    fig, ax = plt.subplots(figsize=(60 * MM, 60 * MM), facecolor="w", edgecolor="k")
    stacked_outcome_series(ax, fractions, annotate=False)
    ax.set_xlabel("Media")
    ax.set_ylabel("Coalescence outcome fraction")
    return save(fig, "fig6c_outcome_fractions")


def main():
    use_paper_style()
    table = outcomes.outcome_table("natural")
    fractions = outcomes.outcome_fractions(table)

    print("Figure 6 - natural sample-derived communities")
    print(fractions.to_string(float_format=lambda x: f"{x:.3f}"))

    panel_b(table)
    panel_c(fractions)


if __name__ == "__main__":
    main()
