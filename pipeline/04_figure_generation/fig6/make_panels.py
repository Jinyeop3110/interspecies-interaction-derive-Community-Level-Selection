#!/usr/bin/env python3
"""Figure 6 panels: coalescence in natural sample-derived communities.

Generates:
  6b  similarity map per medium, with the parental-label-swapped counterpart
  6c  fraction of Dominance outcomes across the nutrient gradient

Run from anywhere:  python figures/fig6/make_panels.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import MEDIUM_COLORS, MM, OUTCOME_COLORS, save, use_paper_style  # noqa: E402

from coalescence import io, outcomes  # noqa: E402


def panel_b(table):
    """Similarity map for the natural sample-derived communities."""
    fig, axes = plt.subplots(1, 3, figsize=(120 * MM, 42 * MM), sharex=True, sharey=True)

    for ax, medium in zip(axes, io.MEDIA_ORDER):
        label = io.MEDIUM_LABELS[medium]
        group = table[table.Medium == medium]

        ax.scatter(group.x2, group.x1, s=8, color="#9a9a9a", alpha=0.55, linewidths=0)
        ax.scatter(group.x1, group.x2, s=8, color=MEDIUM_COLORS[label], linewidths=0)

        theta = np.linspace(0, np.pi / 2, 200)
        radius = 1 / np.sqrt(2)
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), lw=0.7, color="k", ls="--")
        ax.plot([0, 1], [0, 1], lw=0.7, color="k", alpha=0.4)

        ax.set_title(f"{label} (n = {len(group)})", fontsize=7)
        ax.set_xlabel("similarity to parent A")
        ax.set_aspect("equal")
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)

    axes[0].set_ylabel("similarity to parent B")
    return save(fig, "fig6b_similarity_map")


def panel_c(fractions):
    """Dominance fraction across media."""
    fig, ax = plt.subplots(figsize=(50 * MM, 45 * MM))

    labels = list(fractions.index)
    values = fractions["Dominance"].to_numpy(dtype=float)
    ax.bar(labels, values, color=OUTCOME_COLORS["Dominance"], alpha=0.85,
           width=0.6, edgecolor="none")

    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction Dominance")
    return save(fig, "fig6c_dominance_fraction")


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
