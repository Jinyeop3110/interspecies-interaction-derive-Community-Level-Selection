#!/usr/bin/env python3
"""Figure 2 panels: generalized Lotka-Volterra model of community coalescence.

Generates:
  2b  simulated outcomes in similarity space at mu = 0.6
  2c  mean interaction coefficient before and after assembly
  2d  pairwise selection correlation, within vs across parental communities

Everything here is simulated, so this script runs the model rather than reading
archived data. Expect a few minutes at the default replicate count.

    python figures/fig2/make_panels.py --replicates 300
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipeline" / "03_simulation"))

import matplotlib.pyplot as plt  # noqa: E402
from _common import MM, OUTCOME_COLORS, save, use_paper_style  # noqa: E402

from coalescence.decomposition import classify, decompose, retention_and_asymmetry  # noqa: E402
from coalescence.selection_correlation import selection_correlation  # noqa: E402
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

    fig, ax = plt.subplots(figsize=(52 * MM, 52 * MM))
    for outcome, name in enumerate(("Dominance", "Mixture", "Restructuring")):
        selection = coordinates[coordinates[:, 2] == outcome]
        ax.scatter(selection[:, 0], selection[:, 1], s=5, linewidths=0,
                   color=OUTCOME_COLORS[name], label=name, alpha=0.7)

    theta = np.linspace(0, np.pi / 2, 200)
    radius = 1 / np.sqrt(2)
    ax.plot(radius * np.cos(theta), radius * np.sin(theta), lw=0.7, color="k", ls="--")
    ax.plot([0, 1], [0, 1], lw=0.7, color="k", alpha=0.4)

    ax.set_xlabel("similarity to parent A")
    ax.set_ylabel("similarity to parent B")
    ax.set_title(f"$\\mu$ = {MU} (n = {len(coordinates)})", fontsize=7)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=5, loc="lower left")
    return save(fig, "fig2b_similarity_map")


def panel_c(before, after):
    """Mean interaction coefficient before and after assembly.

    Competitive exclusion removes strongly competing species during assembly,
    so the surviving set interacts more weakly than the community as seeded.
    """
    t_statistic, p_value = stats.ttest_rel(before, after)
    print(f"  paired t = {t_statistic:.2f}, df = {len(before) - 1}, p = {p_value:.2e}")
    print(f"  mean alpha before = {before.mean():.4f}, after = {after.mean():.4f}")

    fig, ax = plt.subplots(figsize=(40 * MM, 45 * MM))
    for b, a in zip(before, after):
        ax.plot([0, 1], [b, a], color="#9a9a9a", alpha=0.06, lw=0.5)
    ax.errorbar([0, 1], [before.mean(), after.mean()],
                yerr=[before.std(ddof=1) / np.sqrt(len(before)),
                      after.std(ddof=1) / np.sqrt(len(after))],
                color="#003f5c", marker="o", ms=3, lw=1.2, capsize=2)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["pre-\nassembly", "post-\nassembly"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel(r"mean $\alpha_{ij}$")
    return save(fig, "fig2c_assembly_effect")


def panel_d(events):
    """Pairwise selection correlation, within vs across parental communities."""
    rho_same, rho_cross, delta = selection_correlation(events)
    print(f"  rho_same = {rho_same:.3f}, rho_cross = {rho_cross:.3f}, delta = {delta:.3f}")

    fig, ax = plt.subplots(figsize=(40 * MM, 45 * MM))
    ax.bar(["same\nparent", "different\nparent"], [rho_same, rho_cross],
           color=["#bc5090", "#003f5c"], alpha=0.85, width=0.6, edgecolor="none")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_ylabel("selection correlation $\\rho$")
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
