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
import pandas as pd
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "03_simulation"))

import matplotlib.pyplot as plt  # noqa: E402
from _common import (  # noqa: E402
    MM, OUTCOME_COLORS, density_map, draw_similarity_boundaries, save,
    stacked_outcome_bar, use_paper_style, warn_if_underpowered,
)  # noqa: E402

from coalescence import io  # noqa: E402
from coalescence.decomposition import (  # noqa: E402
    classify, decompose, retention_and_asymmetry,
)
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


def panel_b(events, table):
    """Simulated outcomes as a density field over the similarity map.

    The published panel shows a smoothed density rather than individual points,
    with a relative-density colourbar and the stacked outcome bar beside it.
    """
    coordinates = []
    for parent_a, parent_b, coalesced in events:
        x1, x2, _ = decompose(parent_a, parent_b, coalesced)
        if np.isfinite(x1) and np.isfinite(x2):
            coordinates.append((x1, x2))
    coordinates = np.array(coordinates)

    fig, ax = plt.subplots(figsize=(60 * MM, 60 * MM), facecolor="w", edgecolor="k")
    mappable = density_map(ax, coordinates[:, 0], coordinates[:, 1])
    draw_similarity_boundaries(ax)
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xlabel("Similarity(C,A)")
    ax.set_ylabel("Similarity(C,B)")
    ax.set_title(f"$\\mu$ = {MU}", fontsize=7, style="italic")

    colourbar = fig.colorbar(mappable, ax=ax, fraction=0.046, pad=0.08)
    colourbar.set_label("rel. density", style="italic", fontsize=7)
    colourbar.set_ticks([])
    save(fig, "fig2b_similarity_map")

    # The stacked outcome bar that accompanies the map.
    fig, ax = plt.subplots(figsize=(22 * MM, 60 * MM), facecolor="w", edgecolor="k")
    stacked_outcome_bar(ax, table, annotate="percent")
    return save(fig, "fig2b_outcome_fractions")


def panel_c(before, after):
    """Mean interaction coefficient before and after assembly.

    Competitive exclusion removes strongly competing species during assembly,
    so the surviving set interacts more weakly than the community as seeded.

    Drawn as a jittered strip plot with an open square for the mean and s.e.m.,
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
                    fmt="s", markersize=12, markerfacecolor="white",
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


#: The experimental sub-panel of 2d shows the Base medium, the standard
#: condition the rest of Figures 1 and 2 characterize. Its published values,
#: same +0.20 and cross -0.05, identify it.
EXPERIMENT_MEDIUM = "M"


def _strip_pair(ax, same, cross, ylim, null_level=None, title=None, seed=0):
    """One Same-vs-Cross sub-panel: jittered points, square means, null line."""
    rng = np.random.default_rng(seed)

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

    if null_level is not None:
        ax.axhline(null_level, color="#95a5a6", linestyle="-", linewidth=2,
                   alpha=0.7, zorder=5)

    # Significance bracket over the two categories.
    top = ylim[1] - 0.12 * (ylim[1] - ylim[0])
    ax.plot([0, 0, 1, 1], [top, top + 0.02, top + 0.02, top], color="black",
            linewidth=0.8)
    ax.text(0.5, top + 0.03, "***", ha="center", va="bottom", fontsize=9)

    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Same\nParent", "Cross\nParents"])
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    if title:
        ax.set_title(title, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _event_values(events):
    """Per-event same and cross concordance, over the events that support both."""
    same, cross = [], []
    for parent_a, parent_b, coalesced in events:
        result = event_concordance(parent_a, parent_b, coalesced)
        if result is None:
            continue
        rho_same, rho_cross = result
        if np.isfinite(rho_same) and np.isfinite(rho_cross):
            same.append(rho_same)
            cross.append(rho_cross)
    return np.array(same), np.array(cross)


def experimental_events(medium=EXPERIMENT_MEDIUM):
    """Coalescence events for one experimental medium, exclusions applied."""
    sequences = io.load_sequences("synthetic")
    columns = io.composition_matrix(sequences)
    events = []
    for _, event in io.load_coalescence_events("synthetic").iterrows():
        if event.Medium != medium:
            continue
        ids = {str(event.SampleIDX), str(event.SampleIDX_Sub1),
               str(event.SampleIDX_Sub2)}
        if not ids.isdisjoint(io.EXCLUDED_SAMPLES):
            continue
        triple = tuple(io.sample_vector(sequences, columns, str(s)) for s in
                       (event.SampleIDX_Sub1, event.SampleIDX_Sub2, event.SampleIDX))
        if any(v is None for v in triple):
            continue
        events.append(triple)
    return events


def panel_d(events):
    """Pairwise selection correlation, simulation beside experiment.

    Two sub-panels as published: the simulated events at mu = 0.6, and the
    experimental events in the Base medium. Per-event values as jittered
    points, filled square means with s.e.m., and the shuffled-label null as a
    grey baseline.
    """
    fig, axes = plt.subplots(1, 2, figsize=(3.6, 2.5), facecolor="w", edgecolor="k")

    same, cross = _event_values(events)
    _, _, null = permutation_test(events, n_permutations=200)
    print(f"  simulation: rho_same = {same.mean():.3f}, "
          f"rho_cross = {cross.mean():.3f}, delta = {same.mean() - cross.mean():.3f}")
    _strip_pair(axes[0], same, cross, (-0.8, 0.8), np.nanmean(null),
                title="Simulation")
    axes[0].set_ylabel("Pairwise selection correlation")

    experiment = experimental_events()
    exp_same, exp_cross = _event_values(experiment)
    _, exp_p, exp_null = permutation_test(experiment, n_permutations=200)
    print(f"  experiment ({io.MEDIUM_LABELS[EXPERIMENT_MEDIUM]}, n = {len(exp_same)}): "
          f"rho_same = {exp_same.mean():.3f}, rho_cross = {exp_cross.mean():.3f}, "
          f"delta = {exp_same.mean() - exp_cross.mean():.3f}, p = {exp_p:.3g}")
    print("  published: same +0.20, cross -0.05, delta 0.235")
    _strip_pair(axes[1], exp_same, exp_cross, (-0.4, 0.6), np.nanmean(exp_null),
                title="Experiment", seed=1)

    fig.subplots_adjust(wspace=0.45)
    return save(fig, "fig2d_selection_correlation")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--replicates", type=int, default=300)
    args = parser.parse_args()

    use_paper_style()
    print(f"Figure 2 - simulating {args.replicates} replicates at mu = {MU}")
    events, before, after = run(args.replicates)

    classified = []
    for parent_a, parent_b, coalesced in events:
        x1, x2, _ = decompose(parent_a, parent_b, coalesced)
        if not (np.isfinite(x1) and np.isfinite(x2)):
            continue
        r, asymmetry = retention_and_asymmetry(x1, x2)
        classified.append(classify(r, asymmetry))
    table = pd.DataFrame({"outcome": classified})
    warn_if_underpowered(len(table), "Figure 2")
    counts = table.outcome.value_counts(normalize=True)
    print(f"  n = {len(table)}  Dominance {counts.get(0, 0):.1%}  "
          f"Mixture {counts.get(1, 0):.1%}  Restructuring {counts.get(2, 0):.1%}")
    print("  published: Dominance 61%, Mixture 13%, Restructuring 26%")

    panel_b(events, table)
    panel_c(before, after)
    panel_d(events)


if __name__ == "__main__":
    main()
