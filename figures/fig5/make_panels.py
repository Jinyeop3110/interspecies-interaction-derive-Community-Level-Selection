#!/usr/bin/env python3
"""Figure 5 panels: predictability of Dominance direction.

Generates:
  5b  relative abundance of the dominant species across media

Panel 5c is not generated here; it needs the pairwise invasion assay results,
which are not part of the archived data. See README.md.

Run from anywhere:  python figures/fig5/make_panels.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import MEDIUM_COLORS, MM, save, use_paper_style  # noqa: E402

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


def dominant_fraction_by_medium(metadata, sequences):
    """Mean dominant-species fraction per medium, over parental communities.

    The two replicates of each parental community are averaged before taking
    the mean across communities, so the unit of replication is the community
    (n = 30 per medium) rather than the individual culture.
    """
    columns = io.composition_matrix(sequences)
    rows = {}

    for medium in io.MEDIA_ORDER:
        selection = metadata[
            (metadata.Timepoint == "F")
            & (metadata.CommunityOrigin == "S")
            & (metadata.Medium == medium)
            & (metadata.CoalescenceType == "S")
        ].copy()
        selection["fraction"] = [
            dominant_species_fraction(sequences, columns, str(s))
            for s in selection.SampleIDX
        ]

        per_community = selection.groupby("CommunityIDX")["fraction"].mean().dropna()
        rows[io.MEDIUM_LABELS[medium]] = {
            "mean": per_community.mean(),
            "sem": per_community.std(ddof=1) / np.sqrt(len(per_community)),
            "n": len(per_community),
        }

    return pd.DataFrame(rows).T


def panel_b(summary):
    """Dominant-species abundance across the nutrient gradient."""
    fig, ax = plt.subplots(figsize=(50 * MM, 45 * MM))

    labels = list(summary.index)
    ax.bar(
        labels,
        summary["mean"].to_numpy(dtype=float),
        yerr=summary["sem"].to_numpy(dtype=float),
        color=[MEDIUM_COLORS[label] for label in labels],
        alpha=0.85, width=0.6, edgecolor="none",
        error_kw={"elinewidth": 0.8, "capthick": 0.8, "capsize": 2},
    )

    ax.set_ylim(0, 1)
    ax.set_ylabel("dominant species relative abundance")
    return save(fig, "fig5b_dominant_species_abundance")


def main():
    use_paper_style()
    metadata = io.load_metadata()
    sequences = io.load_sequences("synthetic")

    summary = dominant_fraction_by_medium(metadata, sequences)
    print("Figure 5b - dominant species relative abundance")
    print(summary.to_string(float_format=lambda x: f"{x:.3f}"))

    panel_b(summary)


if __name__ == "__main__":
    main()
