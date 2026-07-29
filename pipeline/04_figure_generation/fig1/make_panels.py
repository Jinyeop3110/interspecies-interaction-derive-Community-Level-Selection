#!/usr/bin/env python3
"""Figure 1 panels: coalescence of synthetic communities frequently yields Dominance.

Generates:
  1e  coalescence outcomes in similarity space, Base medium, and the
      accompanying outcome-fraction bar

Panel 1d is a time course and is NOT generated here; see README.md.

Run from anywhere:  python figures/fig1/make_panels.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import MEDIUM_COLORS, MM, OUTCOME_COLORS, save, use_paper_style  # noqa: E402

from coalescence import io, outcomes  # noqa: E402

#: Figure 1 characterizes the standard condition before the nutrient gradient
#: is introduced in Figure 4.
MEDIUM = "M"


def panel_e(table):
    """Similarity map for the Base medium.

    Grey points are the same events with the two parental labels swapped; the
    map has no intrinsic parent ordering.  See figures/fig4/README.md.
    """
    fig, ax = plt.subplots(figsize=(52 * MM, 52 * MM))

    ax.scatter(table.x2, table.x1, s=9, color="#9a9a9a", alpha=0.55, linewidths=0)
    ax.scatter(table.x1, table.x2, s=9, color=MEDIUM_COLORS["Base"], linewidths=0)

    theta = np.linspace(0, np.pi / 2, 200)
    radius = 1 / np.sqrt(2)
    ax.plot(radius * np.cos(theta), radius * np.sin(theta), lw=0.7, color="k", ls="--")
    ax.plot([0, 1], [0, 1], lw=0.7, color="k", alpha=0.4)

    ax.set_xlabel("similarity to parent A")
    ax.set_ylabel("similarity to parent B")
    ax.set_title(f"Base (n = {len(table)})", fontsize=7)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    return save(fig, "fig1e_similarity_map")


def panel_e_fractions(fractions):
    """Outcome-fraction bar shown alongside the similarity map in panel e."""
    fig, ax = plt.subplots(figsize=(40 * MM, 45 * MM))

    names = ["Dominance", "Mixture", "Restructuring"]
    values = [fractions.loc["Base", name] for name in names]
    ax.bar(names, values, color=[OUTCOME_COLORS[n] for n in names],
           alpha=0.85, width=0.6, edgecolor="none")

    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of events")
    ax.tick_params(axis="x", rotation=45)
    return save(fig, "fig1e_outcome_fractions")


def main():
    use_paper_style()
    table = outcomes.outcome_table("synthetic")
    base_only = table[table.Medium == MEDIUM]
    fractions = outcomes.outcome_fractions(table)

    print(f"Figure 1 - Base medium (n = {len(base_only)})")
    print(fractions.loc[["Base"]].to_string(float_format=lambda x: f"{x:.3f}"))

    panel_e(base_only)
    panel_e_fractions(fractions)


if __name__ == "__main__":
    main()
