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
    MEDIUM_COLORS, MM, OUTCOME_COLORS, SWAPPED_ALPHA, SWAPPED_COLOR,
    draw_similarity_boundaries, save, use_paper_style,
)  # noqa: E402

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

        ax.scatter(group.x1, group.x2, s=25, color=MEDIUM_COLORS[label],
                   alpha=0.7, linewidths=0)
        ax.scatter(group.x2, group.x1, s=25, color=SWAPPED_COLOR,
                   alpha=SWAPPED_ALPHA, linewidths=0)

        draw_similarity_boundaries(ax)
        ax.set_title(f"{label} (n = {len(group)})", fontsize=7)
        ax.set_xlabel("Similarity(C,A)")

    axes[0].set_ylabel("Similarity(C,B)")
    return save(fig, "fig4c_similarity_map")


def panel_d(fractions):
    """Stacked outcome fractions along the nutrient gradient."""
    fig, ax = plt.subplots(figsize=(55 * MM, 45 * MM))

    labels = list(fractions.index)
    bottom = np.zeros(len(labels))
    for outcome in ("Dominance", "Mixture", "Restructuring"):
        values = fractions[outcome].to_numpy(dtype=float)
        ax.bar(labels, values, bottom=bottom, color=OUTCOME_COLORS[outcome],
               alpha=0.85, width=0.65, label=outcome, edgecolor="none")
        bottom += values

    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction")
    ax.legend(fontsize=6, loc="upper right", bbox_to_anchor=(1.45, 1.0))
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
