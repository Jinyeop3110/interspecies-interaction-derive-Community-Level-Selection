#!/usr/bin/env python3
"""Figure 5 panels: predictability of Dominance direction.

Generates:
  5b  relative abundance of the dominant species across media

Panel 5c is not generated here; it needs the pairwise invasion assay results,
which are not part of the archived data. See README.md.

    python pipeline/04_figure_generation/fig5/generate_figure.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import MEDIUM_COLORS, MM, save, use_paper_style  # noqa: E402  # noqa: E402

from coalescence import io  # noqa: E402


def dominant_species_fraction(sequences, columns, sample_idx):
    """Relative abundance of the single most abundant species in one community."""
    if sample_idx not in sequences.index:
        return np.nan
    values = np.nan_to_num(
        sequences.loc[sample_idx, columns].to_numpy(dtype=float).ravel()[: len(columns)]
    )
    total = values.sum()
    return values.max() / total if total > 0 else np.nan


#: Parental community index range for the 12-species richness class. The
#: published panel is restricted to this class.
POOL_12_RANGE = (9, 18)


def dominant_fraction_by_medium(metadata, sequences):
    """Mean dominant-species fraction per medium.

    Restricted to the **12-species parental communities**, and averaged over
    individual cultures rather than over communities, matching
    ``Community_PermutateList("F", "S", medium, "S", 12, -1)`` followed by a
    flat mean and standard error in the original script. Both choices matter:
    pooling all three richness classes, or averaging replicates within a
    community first, shifts the Base value by about nine percentage points.
    """
    columns = io.composition_matrix(sequences)
    low, high = POOL_12_RANGE
    rows, values = {}, {}

    for medium in io.MEDIA_ORDER:
        selection = metadata[
            (metadata.Timepoint == "F")
            & (metadata.CommunityOrigin == "S")
            & (metadata.Medium == medium)
            & (metadata.CoalescenceType == "S")
        ].copy()
        selection = selection[~selection.SampleIDX.astype(str).isin(io.EXCLUDED_SAMPLES)]

        community_index = selection.CommunityIDX.astype(int)
        selection = selection[(community_index > low) & (community_index <= high)]

        fractions = np.array([
            dominant_species_fraction(sequences, columns, str(s))
            for s in selection.SampleIDX
        ])
        fractions = fractions[np.isfinite(fractions)]

        rows[io.MEDIUM_LABELS[medium]] = {
            "mean": fractions.mean(),
            "sem": fractions.std() / np.sqrt(len(fractions)),
            "n": len(fractions),
        }
        values[io.MEDIUM_LABELS[medium]] = fractions

    return pd.DataFrame(rows).T, values


def panel_b(summary, per_community):
    """Dominant-species abundance across the nutrient gradient.

    Pale bars at alpha 0.3 with the individual parental communities jittered
    over them, error bars in a darker shade of the medium colour, and the
    percentage written above each bar -- the idiom of the published panel.
    """
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(35 * MM, 40 * MM), facecolor="w", edgecolor="k")

    labels = list(summary.index)
    positions = np.arange(len(labels))

    for position, label in zip(positions, labels):
        colour = MEDIUM_COLORS[label]
        mean = float(summary.loc[label, "mean"])
        sem = float(summary.loc[label, "sem"])

        ax.bar(position, mean, width=0.8, color=colour, alpha=0.3, linewidth=0.5,
               edgecolor=colour)

        values = per_community[label]
        jitter = 0.1 * (rng.random(len(values)) - 0.5)
        ax.scatter(position + jitter, values, s=0.4, color=colour, alpha=0.3,
                   linewidths=0)

        ax.errorbar(position, mean, yerr=sem, fmt="none", ecolor=colour,
                    elinewidth=1.2, capsize=3, capthick=1.2, alpha=0.5)
        ax.text(position, mean + 0.02, f"{round(mean * 100)}%", ha="center",
                va="bottom", fontsize=6)

    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    for tick, label in zip(ax.get_xticklabels(), labels):
        tick.set_color(MEDIUM_COLORS[label])
        tick.set_style("italic")
    ax.set_ylabel("Dominant species relative abundance")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return save(fig, "fig5b_dominant_species_abundance")


def main():
    use_paper_style()
    metadata = io.load_metadata()
    sequences = io.load_sequences("synthetic")

    summary, per_community = dominant_fraction_by_medium(metadata, sequences)
    print("Figure 5b - dominant species relative abundance")
    print(summary.to_string(float_format=lambda x: f"{x:.3f}"))
    print("  published: 44 +- 2%, 51 +- 5%, 67 +- 4%")

    panel_b(summary, per_community)


if __name__ == "__main__":
    main()
