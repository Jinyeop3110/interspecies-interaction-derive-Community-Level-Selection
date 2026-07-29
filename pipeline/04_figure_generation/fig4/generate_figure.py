#!/usr/bin/env python3
"""Figure 4 panels: nutrient concentration and coalescence outcomes.

Generates:
  4c  similarity map per medium, with the parental-label-swapped counterpart
  4d  outcome fractions across the nutrient gradient

Run from anywhere:  python figures/fig4/make_panels.py
"""

import numpy as np

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import (  # noqa: E402
    MEDIUM_COLORS, MM, POOL_MARKERS, SWAPPED_ALPHA, SWAPPED_COLOR,
    draw_similarity_boundaries, save, stacked_outcome_series, use_paper_style,
)

from coalescence import io, outcomes  # noqa: E402


def panel_c(table):
    """Similarity map, one axis per medium.

    Each event appears twice: once under the original parental labelling
    (coloured) and once with the two parental labels swapped (grey).  The map
    has no intrinsic parent ordering, so plotting both makes the symmetry of
    the space explicit rather than implying a preferred axis.
    """
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 60 * MM), sharex=True, sharey=True)

    for ax, medium in zip(axes, io.MEDIA_ORDER):
        label = io.MEDIUM_LABELS[medium]
        group = table[table.Medium == medium]

        for pool, marker in POOL_MARKERS.items():
            subset = group[group.pool == pool]
            if subset.empty:
                continue
            ax.scatter(subset.x1, subset.x2, s=25, marker=marker,
                       color=MEDIUM_COLORS[label], alpha=0.7, linewidths=0,
                       label=f"Pool {pool}" if medium == io.MEDIA_ORDER[0] else None)
        for pool, marker in POOL_MARKERS.items():
            subset = group[group.pool == pool]
            if subset.empty:
                continue
            ax.scatter(subset.x2, subset.x1, s=25, marker=marker,
                       color=SWAPPED_COLOR, alpha=SWAPPED_ALPHA, linewidths=0)

        draw_similarity_boundaries(ax)
        ax.set_title(label, fontsize=7, style="italic",
                     color=MEDIUM_COLORS[label])
        ax.set_xlabel("Similarity(C,A)")

    axes[0].set_ylabel("Similarity(C,B)")
    axes[0].legend(loc="upper right", fontsize=6, frameon=False)
    return save(fig, "fig4c_similarity_map")


def panel_d(fractions):
    """Outcome fractions across the nutrient gradient.

    A stacked step-fill with dashed separators between media and the
    percentages written into each band, as in the published panel, rather than
    grouped bars with an external legend.
    """
    fig, ax = plt.subplots(figsize=(60 * MM, 60 * MM), facecolor="w", edgecolor="k")
    stacked_outcome_series(ax, fractions)
    ax.set_xlabel("Media")
    ax.set_ylabel("Coalescence outcome fraction")
    return save(fig, "fig4d_outcome_fractions")


def main():
    use_paper_style()
    table = outcomes.outcome_table("synthetic")
    fractions = outcomes.outcome_fractions(table)

    print("Figure 4 - synthetic communities")
    print(fractions.to_string(float_format=lambda x: f"{x:.3f}"))

    panel_c(table)
    panel_d(fractions)


if __name__ == "__main__":
    main()
