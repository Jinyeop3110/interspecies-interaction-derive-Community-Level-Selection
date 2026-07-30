#!/usr/bin/env python3
"""Figure 3 panels: interaction strength controls the outcome transition.

Generates:
  3a  similarity maps with PDI histograms at three representative mu
  3b  outcome fractions across the mu sweep

Everything here is simulated, so this script runs the model. The sweep is the
expensive part; expect several minutes.

    python pipeline/04_figure_generation/fig3/generate_figure.py --replicates 200
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "03_simulation"))

import matplotlib.pyplot as plt  # noqa: E402
from _common import (  # noqa: E402
    DENSITY_RAMP, MM, OUTCOME_COLORS, density_map,
    draw_similarity_boundaries, marginal_pdi_histogram, save, use_paper_style,
)  # noqa: E402

from coalescence.decomposition import parental_dominance_index  # noqa: E402
from coalescence_simulation import outcome_fractions, sweep  # noqa: E402

#: The three interaction strengths shown as maps in panel a.
REPRESENTATIVE_MU = [0.3, 0.6, 0.8]

#: Sweep range for panel b. The endpoint is included, unlike a bare np.arange.
MU_MIN, MU_MAX, MU_STEP = 0.0, 1.2, 0.1


def panel_a(table):
    """Similarity maps at three interaction strengths, with PDI histograms."""
    # Three 60 mm maps with the histogram row beneath, proportioned as the
    # originals (maps 60 mm square, histograms about 22.6 mm tall).
    fig, axes = plt.subplots(2, 3, figsize=(180 * MM, 83 * MM),
                             gridspec_kw={"height_ratios": [60, 22.6]},
                             facecolor="w", edgecolor="k")
    fig.subplots_adjust(wspace=0.3, hspace=0.45)

    for column, mu in enumerate(REPRESENTATIVE_MU):
        events = table[np.isclose(table.mu, mu)]
        ax = axes[0, column]

        for outcome, name in enumerate(("Dominance", "Mixture", "Restructuring")):
            selection = events[events.outcome == outcome]
            ax.scatter(selection.x1, selection.x2, s=3, linewidths=0,
                       color=OUTCOME_COLORS[name], alpha=0.6)

        draw_similarity_boundaries(ax)
        ax.set_title(f"$\\mu$ = {mu}", fontsize=7, style="italic")
        if column == 0:
            ax.set_ylabel("Similarity(C,B)")

        # PDI histogram: unimodal at weak interactions, bimodal at strong.
        histogram = axes[1, column]
        pdi = parental_dominance_index(events.x1.to_numpy(), events.x2.to_numpy())
        marginal_pdi_histogram(histogram, pdi, 1 - pdi, DENSITY_RAMP[0])
        histogram.set_xlabel("Parental Dominance Index")
        histogram.set_xlim(0, 1)
        if column == 0:
            histogram.set_ylabel("Density")

    return save(fig, "fig3a_similarity_maps")


def panel_b(table):
    """Outcome fractions across the interaction-strength sweep.

    Stacked step-fill with the three representative interaction strengths
    marked, and the class letters written into the bands, as published.
    """
    fractions = outcome_fractions(table)
    print(fractions.to_string(float_format=lambda x: f"{x:.3f}"))

    fig, ax = plt.subplots(figsize=(60 * MM, 60 * MM), facecolor="w", edgecolor="k")
    mu_values = fractions.index.to_numpy(dtype=float)

    # step="post" holds each value until the next x, so append a trailing
    # point or the final mu band is never drawn.
    step = mu_values[-1] - mu_values[-2]
    x = np.append(mu_values, mu_values[-1] + step)
    bottom = np.zeros(len(x))
    midpoint = len(mu_values) // 2
    for name in ("Dominance", "Mixture", "Restructuring"):
        values = np.append(fractions[name].to_numpy(dtype=float), 0.0)
        ax.fill_between(x, bottom, bottom + values, step="post",
                        color=OUTCOME_COLORS[name], alpha=0.85, linewidth=0)
        if values[midpoint] > 0.08:
            ax.text(mu_values[midpoint], bottom[midpoint] + values[midpoint] / 2,
                    name[0], ha="center", va="center", fontsize=8,
                    style="italic")
        bottom = bottom + values

    for mu in REPRESENTATIVE_MU:
        ax.axvline(mu, color="black", linestyle="--", alpha=0.5, linewidth=1)
        ax.text(mu, 1.02, f"{mu}", ha="center", va="bottom", fontsize=6)

    ax.set_xlabel(r"Interaction Strength $\mu$")
    ax.set_ylabel("Fraction")
    ax.set_xlim(mu_values.min(), mu_values.max() + step)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"])
    return save(fig, "fig3b_outcome_fractions")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--communities", type=int, default=2)
    parser.add_argument("--species", type=int, default=12)
    args = parser.parse_args()

    use_paper_style()

    # Include the endpoint: np.arange(0, 1.2, 0.1) would stop at 1.1.
    mu_values = np.round(np.arange(MU_MIN, MU_MAX + MU_STEP / 2, MU_STEP), 10)
    mu_values = np.unique(np.concatenate([mu_values, REPRESENTATIVE_MU]))

    print(f"Figure 3 - sweeping {len(mu_values)} values of mu, "
          f"{args.replicates} replicates each")
    table = sweep(mu_values, args.replicates, args.communities, args.species)

    panel_a(table)
    panel_b(table)


if __name__ == "__main__":
    main()
