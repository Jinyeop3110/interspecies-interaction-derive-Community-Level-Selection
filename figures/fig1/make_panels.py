#!/usr/bin/env python3
"""Figure 1 panels: coalescence of synthetic communities frequently yields Dominance.

Generates:
  1d  composition of representative coalescence events
  1e  coalescence outcomes in similarity space, Base medium
  1f  outcome fractions, Base medium

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


#: Colours for the coalesced bar, by which parent each species came from.
ORIGIN_COLORS = {
    "A only": "#bc5090",
    "B only": "#003f5c",
    "both": "#8a8a8a",
    "neither": "#ffa600",
}

#: A species must exceed this relative abundance in a parent to count as
#: present there, matching the finest threshold used in the summary tables.
PRESENCE_THRESHOLD = 0.001


def _stacked_composition(ax, vector, colors, width=0.6, x=0.0):
    """Draw one community as a stacked bar, largest species at the bottom."""
    order = np.argsort(vector)[::-1]
    bottom = 0.0
    for index in order:
        height = vector[index]
        if height <= 0:
            continue
        ax.bar(x, height, bottom=bottom, width=width, color=colors[index],
               linewidth=0)
        bottom += height


def panel_d(events, sequences, columns):
    """Composition of representative coalescence events.

    For each event, three stacked bars: the two parental communities and the
    coalesced community.  Species in the coalesced bar are coloured by which
    parent they came from, which is what makes Dominance visible — one parent's
    species make up nearly all of the coalesced community.
    """
    from coalescence import io as _io

    for label, event in events.items():
        a = _io.sample_vector(sequences, columns, str(event.SampleIDX))
        a1 = _io.sample_vector(sequences, columns, str(event.SampleIDX_Sub1))
        a2 = _io.sample_vector(sequences, columns, str(event.SampleIDX_Sub2))
        if a is None or a1 is None or a2 is None:
            print(f"  skipped {label}: missing composition data")
            continue

        in_1 = a1 > PRESENCE_THRESHOLD
        in_2 = a2 > PRESENCE_THRESHOLD
        origin = np.where(in_1 & in_2, "both",
                          np.where(in_1, "A only",
                                   np.where(in_2, "B only", "neither")))
        coalesced_colors = [ORIGIN_COLORS[o] for o in origin]

        fig, ax = plt.subplots(figsize=(38 * MM, 42 * MM))
        _stacked_composition(ax, a1, [ORIGIN_COLORS["A only"]] * len(a1), x=0)
        _stacked_composition(ax, a2, [ORIGIN_COLORS["B only"]] * len(a2), x=1)
        _stacked_composition(ax, a, coalesced_colors, x=2)

        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["A", "B", "C"])
        ax.set_ylim(0, 1)
        ax.set_ylabel("relative abundance")
        ax.set_title(label, fontsize=7)
        save(fig, f"fig1d_composition_{label.lower()}")


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


def panel_f(fractions):
    """Outcome fractions in the Base medium."""
    fig, ax = plt.subplots(figsize=(40 * MM, 45 * MM))

    names = ["Dominance", "Mixture", "Restructuring"]
    values = [fractions.loc["Base", name] for name in names]
    ax.bar(names, values, color=[OUTCOME_COLORS[n] for n in names],
           alpha=0.85, width=0.6, edgecolor="none")

    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of events")
    ax.tick_params(axis="x", rotation=45)
    return save(fig, "fig1f_outcome_fractions")


def pick_representatives(table, events, sequences, columns):
    """Choose one Dominance and one Restructuring event to illustrate.

    Both are selected deterministically from the Base medium, and both are
    chosen for how clearly they show the defining feature of their class:

    * Dominance — the highest angular asymmetry, so one parent's species make
      up nearly all of the coalesced community.
    * Restructuring — the largest abundance of species present in neither
      parent, since low retention is what defines the class but novel species
      are what make it legible in a composition plot.

    Selecting on novelty matters because the experimental parental communities
    are drawn from a shared 54-isolate pool and therefore overlap heavily.  An
    event picked purely on lowest retention is usually dominated by species
    shared between both parents, which renders as one flat colour and shows
    the reader nothing.
    """
    from coalescence import io as _io

    base = table[table.Medium == MEDIUM]
    by_sample = events.copy()
    by_sample.index = by_sample.SampleIDX.astype(str)

    def novel_fraction(sample_idx):
        if str(sample_idx) not in by_sample.index:
            return -1.0
        event = by_sample.loc[str(sample_idx)]
        a = _io.sample_vector(sequences, columns, str(event.SampleIDX))
        a1 = _io.sample_vector(sequences, columns, str(event.SampleIDX_Sub1))
        a2 = _io.sample_vector(sequences, columns, str(event.SampleIDX_Sub2))
        if a is None or a1 is None or a2 is None or a.sum() <= 0:
            return -1.0
        novel = ~((a1 > PRESENCE_THRESHOLD) | (a2 > PRESENCE_THRESHOLD))
        return a[novel].sum() / a.sum()

    chosen = {}

    dominance = base[base.outcome == 0]
    if len(dominance):
        chosen["Dominance"] = dominance.loc[dominance.asymmetry.idxmax()]

    restructuring = base[base.outcome == 2].copy()
    if len(restructuring):
        restructuring["novel"] = [novel_fraction(s) for s in restructuring.SampleIDX]
        chosen["Restructuring"] = restructuring.loc[restructuring.novel.idxmax()]

    return {
        label: by_sample.loc[str(row.SampleIDX)]
        for label, row in chosen.items()
        if str(row.SampleIDX) in by_sample.index
    }


def main():
    use_paper_style()
    table = outcomes.outcome_table("synthetic")
    base_only = table[table.Medium == MEDIUM]
    fractions = outcomes.outcome_fractions(table)

    print(f"Figure 1 - Base medium (n = {len(base_only)})")
    print(fractions.loc[["Base"]].to_string(float_format=lambda x: f"{x:.3f}"))

    events = io.load_coalescence_events("synthetic")
    sequences = io.load_sequences("synthetic")
    columns = io.composition_matrix(sequences)
    representatives = pick_representatives(table, events, sequences, columns)
    print("representative events:",
          {k: v.SampleIDX for k, v in representatives.items()})

    panel_d(representatives, sequences, columns)
    panel_e(base_only)
    panel_f(fractions)


if __name__ == "__main__":
    main()
