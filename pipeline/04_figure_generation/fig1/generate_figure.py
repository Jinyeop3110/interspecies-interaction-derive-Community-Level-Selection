#!/usr/bin/env python3
"""Figure 1 panels: coalescence of synthetic communities frequently yields Dominance.

Generates:
  1e  outcomes in similarity space, Base medium, plus the stacked outcome bar
      drawn beside it

Panel 1d is a time course across seven growth-dilution cycles and is NOT
generated here; the archived tables carry only the final timepoint. See
README.md.

    python pipeline/04_figure_generation/fig1/generate_figure.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import (  # noqa: E402
    MEDIUM_COLORS, MM, POOL_ALPHAS, POOL_MARKERS, SWAPPED_ALPHA, SWAPPED_COLOR,
    draw_similarity_boundaries, save, stacked_outcome_bar, use_paper_style,
)

from coalescence import outcomes  # noqa: E402

#: Figure 1 characterizes the standard condition before the nutrient gradient
#: is introduced in Figure 4.
MEDIUM = "M"


def panel_e_map(table):
    """Similarity map for the Base medium.

    Marker shape encodes the initial species-pool size, as the caption
    specifies: circles for 6 species, squares for 12, triangles for 24.

    Every event is plotted twice, once under the original parental labelling
    and once with the labels swapped. The map has no intrinsic ordering of the
    two parents, so showing both makes that symmetry explicit. Coloured points
    are drawn first and the grey twins over them, as in the originals.
    """
    fig, ax = plt.subplots(figsize=(60 * MM, 60 * MM), facecolor="w", edgecolor="k")

    # Alpha deepens with pool size, as in the published panel.
    for (pool, marker), alpha in zip(POOL_MARKERS.items(), POOL_ALPHAS):
        group = table[table.pool == pool]
        if group.empty:
            continue
        ax.scatter(group.x1, group.x2, s=25, marker=marker,
                   color=MEDIUM_COLORS["Base"], alpha=alpha, linewidths=0,
                   label=f"Pool {pool}", zorder=2)
    for pool, marker in POOL_MARKERS.items():
        group = table[table.pool == pool]
        if group.empty:
            continue
        ax.scatter(group.x2, group.x1, s=25, marker=marker,
                   color=SWAPPED_COLOR, alpha=SWAPPED_ALPHA, linewidths=0,
                   zorder=1)

    draw_similarity_boundaries(ax)
    ax.set_xlabel("Similarity(C,A)")
    ax.set_ylabel("Similarity(C,B)")
    ax.set_title("Coalescence Outcomes", fontsize=8)
    ax.legend(loc="upper right", fontsize=6, frameon=True,
              framealpha=1.0, edgecolor="0.7")
    return save(fig, "fig1e_similarity_map")


def panel_e_bar(table):
    """Stacked outcome bar for the Base medium, annotated with counts."""
    fig, ax = plt.subplots(figsize=(22 * MM, 60 * MM), facecolor="w", edgecolor="k")
    stacked_outcome_bar(ax, table, annotate="count")
    return save(fig, "fig1e_outcome_fractions")


def main():
    use_paper_style()
    table = outcomes.outcome_table("synthetic")
    base_only = table[table.Medium == MEDIUM]

    counts = base_only.outcome.value_counts()
    print(f"Figure 1 - Base medium, n = {len(base_only)}")
    for code, name in ((0, "Dominance"), (2, "Restructuring"), (1, "Mixture")):
        n = int(counts.get(code, 0))
        print(f"  {name:14s} n = {n:3d}  ({n / len(base_only):.0%})")
    print("  published:       Dominance 54 (65%), Restructuring 26 (31%), "
          "Mixture 3 (4%)")
    print(f"  pool sizes:      {base_only.pool.value_counts().sort_index().to_dict()}")

    panel_e_map(base_only)
    panel_e_bar(base_only)


if __name__ == "__main__":
    main()
