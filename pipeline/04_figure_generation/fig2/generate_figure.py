#!/usr/bin/env python3
"""Figure 2 panels: generalized Lotka-Volterra model of community coalescence.

Generates:
  2b  simulated outcomes in similarity space at mu = 0.6
  2c  mean interaction coefficient before and after assembly
  2d  pairwise selection correlation, within vs across parental communities

Everything here is simulated, so this script runs the model rather than reading
archived data. Expect a few minutes at the default replicate count.

    python pipeline/04_figure_generation/fig2/generate_figure.py --replicates 300
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "03_simulation"))

import matplotlib.pyplot as plt  # noqa: E402
from _common import (  # noqa: E402
    MM, OUTCOME_COLORS, draw_similarity_boundaries, save, use_paper_style,
)  # noqa: E402

from coalescence.decomposition import classify, decompose, retention_and_asymmetry  # noqa: E402
from coalescence.selection_correlation import (  # noqa: E402
    event_concordance, permutation_test,
)
from coalescence_simulation import (  # noqa: E402
    mean_interaction_before_after, simulate_one_replicate,
)

#: Representative interaction strength used for Figure 2.
MU = 0.6


def run(replicates, n_communities=2, species_per_community=12):
    """Simulate once and collect everything the three panels need."""
    events, before, after = [], [], []

    for replicate in range(replicates):
        pairs, interactions, library = simulate_one_replicate(
            MU, replicate, n_communities, species_per_community,
            return_interactions=True,
        )
        events.extend(pairs)

        # Recover each parental community's equilibrium once. Every pair shares
        # the same parents, so take the first appearance of each.
        parents = [None] * n_communities
        index = 0
        for a in range(n_communities):
            for b in range(a + 1, n_communities):
                parents[a] = parents[a] if parents[a] is not None else pairs[index][0]
                parents[b] = parents[b] if parents[b] is not None else pairs[index][1]
                index += 1

        for seeded, survived in mean_interaction_before_after(
            interactions, library, parents
        ):
            before.append(seeded)
            after.append(survived)

    return events, np.array(before), np.array(after)


def panel_b(events):
    """Simulated outcomes in the two-parent similarity map."""
    coordinates = []
    for parent_a, parent_b, coalesced in events:
        x1, x2, _ = decompose(parent_a, parent_b, coalesced)
        r, asymmetry = retention_and_asymmetry(x1, x2)
        coordinates.append((x1, x2, classify(r, asymmetry)))
    coordinates = np.array(coordinates)

    fig, ax = plt.subplots(figsize=(60 * MM, 60 * MM), facecolor="w", edgecolor="k")
    for outcome, name in enumerate(("Dominance", "Mixture", "Restructuring")):
        selection = coordinates[coordinates[:, 2] == outcome]
        ax.scatter(selection[:, 0], selection[:, 1], s=5, linewidths=0,
                   color=OUTCOME_COLORS[name], label=name, alpha=0.7)

    draw_similarity_boundaries(ax)

    ax.set_xlabel("Similarity(C,A)")
    ax.set_ylabel("Similarity(C,B)")
    ax.set_title(f"$\\mu$ = {MU}", fontsize=7, style="italic")
    ax.legend(fontsize=5, loc="lower left")
    return save(fig, "fig2b_similarity_map")


def panel_c(before, after):
    """Mean interaction coefficient before and after assembly.

    Competitive exclusion removes strongly competing species during assembly,
    so the surviving set interacts more weakly than the community as seeded.

    Drawn as a jittered strip plot with a filled square for the mean and s.e.m.,
    matching the published panel, rather than connecting lines.
    """
    t_statistic, p_value = stats.ttest_rel(before, after)
    print(f"  paired t = {t_statistic:.2f}, df = {len(before) - 1}, p = {p_value:.2e}")
    print(f"  mean alpha before = {before.mean():.4f}, after = {after.mean():.4f}")
    print("  published: t_599 = 29.26, p = 1.54e-117")

    rng = np.random.default_rng(42)
    # The original specifies this panel in inches, not mm.
    fig, ax = plt.subplots(figsize=(2.2, 2.2), facecolor="w", edgecolor="k")

    for position, values, colour in ((0, before, "#8B7AB8"), (1, after, "#F4A582")):
        jitter = rng.normal(0, 0.1, len(values))
        ax.scatter(position + jitter, values, s=15, color=colour, alpha=0.3,
                   edgecolors="none")
        ax.errorbar(position, values.mean(),
                    yerr=values.std(ddof=1) / np.sqrt(len(values)),
                    fmt="s", markersize=12, color=colour,
                    markeredgecolor="black", markeredgewidth=0.5,
                    ecolor="black", capsize=5, capthick=1.5, linewidth=1.5,
                    zorder=10)

    # Significance bracket over the two categories.
    top = max(before.max(), after.max())
    ax.plot([0, 0, 1, 1], [top * 1.02, top * 1.06, top * 1.06, top * 1.02],
            color="black", linewidth=1)
    ax.text(0.5, top * 1.07, "***", ha="center", va="bottom", fontsize=14,
            fontweight="bold")

    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Initial", "Post-assembly"])
    ax.set_yticks([0, 0.4, 0.8])
    ax.set_ylabel("Mean pairwise interaction")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return save(fig, "fig2c_assembly_effect")


def panel_d(events):
    """Pairwise selection correlation, within vs across parental communities.

    Per-event values as a jittered strip plot, filled square means with s.e.m.,
    and the shuffled-label null as a grey baseline -- the published idiom.
    """
    same, cross = [], []
    for parent_a, parent_b, coalesced in events:
        rho_same, rho_cross = event_concordance(parent_a, parent_b, coalesced)
        if np.isfinite(rho_same):
            same.append(rho_same)
        if np.isfinite(rho_cross):
            cross.append(rho_cross)
    same, cross = np.array(same), np.array(cross)

    observed, p_value, null = permutation_test(events, n_permutations=200)
    print(f"  rho_same = {same.mean():.3f}, rho_cross = {cross.mean():.3f}, "
          f"delta = {observed:.3f}, permutation p = {p_value:.3g}")

    rng = np.random.default_rng(0)
    # As with panel c, the original specifies inches.
    fig, ax = plt.subplots(figsize=(2.5, 2.5), facecolor="w", edgecolor="k")

    for position, values, colour in ((0, same, "#e74c3c"), (1, cross, "#3498db")):
        jitter = rng.normal(0, 0.08, len(values))
        ax.scatter(position + jitter, values, s=15, color=colour, alpha=0.3,
                   edgecolors="none")
        ax.errorbar(position, values.mean(),
                    yerr=values.std(ddof=1) / np.sqrt(len(values)),
                    fmt="s", markersize=6, color=colour,
                    markeredgecolor="black", markeredgewidth=0.5,
                    ecolor="black", capsize=5, capthick=1.5, linewidth=1.5,
                    zorder=10)

    ax.axhline(np.nanmean(null), color="#95a5a6", linestyle="-", linewidth=2,
               alpha=0.7, zorder=5)

    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.2, 0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Same\nParent", "Cross\nParents"])
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.set_ylabel("Pairwise selection correlation")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return save(fig, "fig2d_selection_correlation")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--replicates", type=int, default=300)
    args = parser.parse_args()

    use_paper_style()
    print(f"Figure 2 - simulating {args.replicates} replicates at mu = {MU}")
    events, before, after = run(args.replicates)

    panel_b(events)
    panel_c(before, after)
    panel_d(events)


if __name__ == "__main__":
    main()
