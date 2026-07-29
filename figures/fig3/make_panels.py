#!/usr/bin/env python3
"""Figure 3 panels: interaction strength controls the outcome transition.

Generates:
  3a  similarity maps with PDI histograms at three representative mu
  3b  outcome fractions across the mu sweep

Everything here is simulated, so this script runs the model. The sweep is the
expensive part; expect several minutes.

    python figures/fig3/make_panels.py --replicates 200
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipeline" / "03_simulation"))

import matplotlib.pyplot as plt  # noqa: E402
from _common import MM, OUTCOME_COLORS, save, use_paper_style  # noqa: E402

from coalescence.decomposition import parental_dominance_index  # noqa: E402
from coalescence_simulation import outcome_fractions, sweep  # noqa: E402

#: The three interaction strengths shown as maps in panel a.
REPRESENTATIVE_MU = [0.3, 0.6, 0.8]

#: Sweep range for panel b. The endpoint is included, unlike a bare np.arange.
MU_MIN, MU_MAX, MU_STEP = 0.0, 1.2, 0.1


def panel_a(table):
    """Similarity maps at three interaction strengths, with PDI histograms."""
    fig, axes = plt.subplots(2, 3, figsize=(120 * MM, 62 * MM),
                             gridspec_kw={"height_ratios": [3, 1]})

    for column, mu in enumerate(REPRESENTATIVE_MU):
        events = table[np.isclose(table.mu, mu)]
        ax = axes[0, column]

        for outcome, name in enumerate(("Dominance", "Mixture", "Restructuring")):
            selection = events[events.outcome == outcome]
            ax.scatter(selection.x1, selection.x2, s=3, linewidths=0,
                       color=OUTCOME_COLORS[name], alpha=0.6)

        theta = np.linspace(0, np.pi / 2, 200)
        radius = 1 / np.sqrt(2)
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), lw=0.7, color="k", ls="--")
        ax.plot([0, 1], [0, 1], lw=0.7, color="k", alpha=0.4)
        ax.set_title(f"$\\mu$ = {mu}  (n = {len(events)})", fontsize=7)
        ax.set_aspect("equal")
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)
        if column == 0:
            ax.set_ylabel("similarity to parent B")

        # PDI histogram: unimodal at weak interactions, bimodal at strong.
        histogram = axes[1, column]
        pdi = parental_dominance_index(events.x1.to_numpy(), events.x2.to_numpy())
        histogram.hist(pdi[np.isfinite(pdi)], bins=20, range=(0, 1),
                       color="#7a7a7a", linewidth=0)
        histogram.set_xlabel("PDI")
        histogram.set_xlim(0, 1)
        if column == 0:
            histogram.set_ylabel("count")

    return save(fig, "fig3a_similarity_maps")


def panel_b(table):
    """Outcome fractions across the interaction-strength sweep."""
    fractions = outcome_fractions(table)
    print(fractions.to_string(float_format=lambda x: f"{x:.3f}"))

    fig, ax = plt.subplots(figsize=(62 * MM, 45 * MM))
    mu_values = fractions.index.to_numpy(dtype=float)

    bottom = np.zeros(len(mu_values))
    for name in ("Dominance", "Mixture", "Restructuring"):
        values = fractions[name].to_numpy(dtype=float)
        ax.fill_between(mu_values, bottom, bottom + values, step="post",
                        color=OUTCOME_COLORS[name], alpha=0.7, label=name,
                        linewidth=0)
        bottom += values

    ax.set_xlabel(r"interaction strength $\mu$")
    ax.set_ylabel("fraction of events")
    ax.set_xlim(mu_values.min(), mu_values.max())
    ax.set_ylim(0, 1)
    ax.legend(fontsize=5, loc="center left", bbox_to_anchor=(1.02, 0.5))
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
