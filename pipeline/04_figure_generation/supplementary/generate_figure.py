#!/usr/bin/env python3
"""Supplementary Fig. 14: parental OD controls and OD-pH association.

Generates the active Supplementary Fig. 14 PDF used by the LaTeX manuscript:
top row, signed parental OD difference versus PDI; bottom row,
parental-community endpoint OD versus endpoint pH.

Run from anywhere:
  python figures/supplementary/make_supp_fig14_od_ph.py
"""

import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
from _common import MEDIUM_COLORS, MM, use_paper_style  # noqa: E402
from coalescence import io, outcomes  # noqa: E402


#: Panels are written next to the script, like every other figure directory.
#: Set COALESCENCE_FIGURE_OUT to send them somewhere else.
OUT_BASE = Path(
    os.environ.get("COALESCENCE_FIGURE_OUT",
                   Path(__file__).resolve().parent / "panels")
) / "supp_fig14_od_vs_ph"


def format_p(value):
    if value < 1e-3:
        return f"{value:.1e}"
    return f"{value:.3f}"


def parent_measurements():
    metadata = io.load_metadata()
    parents = metadata[
        (metadata.Timepoint == "F")
        & (metadata.CommunityOrigin == "S")
        & (metadata.CoalescenceType == "S")
    ].copy()
    return parents.set_index("SampleIDX")


def analysis_table():
    table = outcomes.outcome_table("synthetic")
    event_parents = io.load_coalescence_events("synthetic")[
        ["SampleIDX", "SampleIDX_Sub1", "SampleIDX_Sub2"]
    ]
    table = table.merge(event_parents, on="SampleIDX", how="left")
    parents = parent_measurements()

    rows = []
    for _, event in table.iterrows():
        if event.SampleIDX_Sub1 not in parents.index or event.SampleIDX_Sub2 not in parents.index:
            continue
        parent_a = parents.loc[event.SampleIDX_Sub1]
        parent_b = parents.loc[event.SampleIDX_Sub2]
        rows.append(
            {
                "Medium": event.Medium,
                "PDI": event.PDI,
                "delta_OD": parent_a.fieldOD7 - parent_b.fieldOD7,
            }
        )
    return pd.DataFrame(rows)


def annotate_spearman(ax, x, y, loc="upper left"):
    valid = np.isfinite(x) & np.isfinite(y)
    rho, p_value = spearmanr(x[valid], y[valid])
    text = f"n = {valid.sum()}\n$\\rho$ = {rho:.2f}\n$p$ = {format_p(p_value)}"
    anchor = {
        "upper left": (0.05, 0.95, "left", "top"),
        "lower right": (0.95, 0.05, "right", "bottom"),
    }[loc]
    ax.text(
        anchor[0],
        anchor[1],
        text,
        transform=ax.transAxes,
        ha=anchor[2],
        va=anchor[3],
        fontsize=6.5,
    )


def symmetric_limits(values, pad=0.08):
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return (-1, 1)
    bound = np.nanmax(np.abs(finite)) * (1 + pad)
    if bound == 0:
        bound = 1
    return (-bound, bound)


def plot():
    use_paper_style()
    table = analysis_table()
    parents = parent_measurements().reset_index()

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(150 * MM, 84 * MM),
        sharey="row",
        constrained_layout=False,
    )

    for col, medium in enumerate(io.MEDIA_ORDER):
        label = io.MEDIUM_LABELS[medium]
        color = MEDIUM_COLORS[label]
        group = table[table.Medium == medium].copy()

        x_od = group.delta_OD.to_numpy(dtype=float)
        y_pdi = group.PDI.to_numpy(dtype=float)

        ax = axes[0, col]
        ax.scatter(-x_od, 1 - y_pdi, s=8, color="#bdbdbd", alpha=0.45, linewidths=0)
        ax.scatter(x_od, y_pdi, s=9, color=color, alpha=0.82, linewidths=0)
        ax.axvline(0, color="0.6", lw=0.6, ls="--")
        ax.axhline(0.5, color="0.6", lw=0.6, ls="--")
        ax.set_xlim(symmetric_limits(x_od))
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(label, fontsize=7, color=color)
        annotate_spearman(ax, x_od, y_pdi, loc="upper left")
        if col == 0:
            ax.set_ylabel("PDI")

        parent_group = parents[parents.Medium == medium].copy()
        parent_od = parent_group.fieldOD7.to_numpy(dtype=float)
        parent_ph = parent_group.fieldPH7.to_numpy(dtype=float)

        ax = axes[1, col]
        ax.scatter(parent_od, parent_ph, s=9, color=color, alpha=0.82, linewidths=0)
        annotate_spearman(ax, parent_od, parent_ph, loc="lower right")
        ax.set_xlabel(r"Endpoint OD$_{600}$")
        if col == 0:
            ax.set_ylabel("Endpoint pH")

    fig.subplots_adjust(left=0.08, right=0.99, top=0.91, bottom=0.14, wspace=0.28, hspace=0.36)
    for suffix in [".pdf"]:
        out = OUT_BASE.with_suffix(suffix)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, bbox_inches="tight")
        print(f"Wrote {out}")
    plt.close(fig)

    print("\nSpearman summaries")
    for medium in io.MEDIA_ORDER:
        label = io.MEDIUM_LABELS[medium]
        group = table[table.Medium == medium]
        valid = np.isfinite(group.delta_OD) & np.isfinite(group.PDI)
        rho, p_value = spearmanr(group.loc[valid, "delta_OD"], group.loc[valid, "PDI"])
        print(f"{label:5s} delta_OD vs PDI: n={valid.sum():2d}, rho={rho:.3f}, p={p_value:.3g}")

        parent_group = parents[parents.Medium == medium]
        valid = np.isfinite(parent_group.fieldOD7) & np.isfinite(parent_group.fieldPH7)
        rho, p_value = spearmanr(
            parent_group.loc[valid, "fieldOD7"], parent_group.loc[valid, "fieldPH7"]
        )
        print(f"{label:5s} parent OD vs pH: n={valid.sum():2d}, rho={rho:.3f}, p={p_value:.3g}")


if __name__ == "__main__":
    plot()
