#!/usr/bin/env python3
"""
analyze_additive_null.py
========================
Per-event additive null model analysis for Reviewer 3, Point #1.

For each coalescence event, we compare the OBSERVED outcome classification
(Dominance / Mixture / Restructuring) to a NULL expectation obtained by
assuming purely additive (50:50) mixing of the two parental communities:

    n_C_null = normalize(n_A + n_B)

Both observed n_C and null n_C_null are classified through the same
vector-decomposition pipeline used throughout the manuscript.

Output figures (saved as SVG, PDF, PNG):
  1. Sankey-style alluvial / paired comparison (null vs observed classification)
  2. Histogram of PDI_observed - PDI_null
  3. Contingency-table heatmap (observed x null)
  4. Per-medium breakdown panels
"""

import sys
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, BASE_DIR)  # so we can import common_setup, COLORMAP

from COLORMAP import (
    PHASE_DIAGRAM_COLORS,
    get_medium_color,
)

# ---------------------------------------------------------------------------
# Import canonical analysis functions from common_setup.py
# (sys.path already configured above so common_setup is importable)
# ---------------------------------------------------------------------------
from common_setup import (
    normalize,
    metric_VectorDecomposition_onlyPositive,
    calculate_assymetricity,
    characterize_case,
)


CLASS_NAMES = {0: "Dominance", 1: "Mixture", 2: "Restructuring"}
CLASS_COLORS = {
    0: PHASE_DIAGRAM_COLORS["Dominance"],
    1: PHASE_DIAGRAM_COLORS["Mixing"],
    2: PHASE_DIAGRAM_COLORS["Restructuring"],
}

MEDIUM_MAP = {"L": "LN", "M": "MN", "H": "HN"}

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
sns.set_style("ticks")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["figure.dpi"] = 200
mpl.rcParams["font.size"] = 8
mpl.rcParams["axes.linewidth"] = 0.5
mpl.rcParams["xtick.minor.width"] = 0.4
mpl.rcParams["xtick.major.width"] = 0.5
mpl.rcParams["ytick.minor.width"] = 0.4
mpl.rcParams["ytick.major.width"] = 0.5
plt.rcParams["text.usetex"] = False
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"

SAVE_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# Data comes from the archive directory; set COALESCENCE_DATA to relocate it.
from coalescence import io as _io

DATA = _io.data_dir()
coalescence_path = DATA / "processed_CoalescenceEvent_synthetic.xlsx"
sequences_path = DATA / "processed_Sequences_synthetic.xlsx"
raw_count_path = DATA / "M_OTUtableGreenGenes.csv"

print("Loading data...")
coal_df = pd.read_excel(coalescence_path)
seq_df = pd.read_excel(sequences_path)
raw_count_df = pd.read_csv(raw_count_path)

# Build a lookup: SampleIDX -> abundance vector (numpy array, float64)
seq_lookup = {}
for _, row in seq_df.iterrows():
    sid = row["SampleIDX"]
    vec = row.iloc[1:].values.astype(float)
    seq_lookup[sid] = vec

raw_ids = raw_count_df.iloc[:, 0].astype(str).str.replace("_F_filt.fastq.gz", "", regex=False)
raw_values = raw_count_df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0)
raw_lookup = {
    sid: raw_values.iloc[i].to_numpy(float)
    for i, sid in enumerate(raw_ids)
}

# Quality-control exclusions; single definition in coalescence.io.
exception_list = set(_io.EXCLUDED_SAMPLES)
_unused_exception_list = set(
    ["P4-02", "P4-03", "P4-23", "P4-24", "P7-97", "P8-12"]
    + ["P8-91"]
    + ["P5-73", "P5-69", "P5-64", "P5-61", "P5-59", "P5-56"]
    + ["P5-47", "P5-50"]
    + ["P5-39", "P5-87", "P5-54", "P6-02", "P6-47", "P6-74", "P6-57"]
)


# ---------------------------------------------------------------------------
# Per-event analysis
# ---------------------------------------------------------------------------
def classify_vector(n_C, n_A, n_B):
    """Classify a composition n_C relative to parents n_A, n_B.
    Returns (class_int, x_coord, y_coord, u_coeff, v_coeff, k_coeff).
    Returns None if the computation fails (zero vectors etc.)."""
    n_A = np.array(n_A, dtype=float)
    n_B = np.array(n_B, dtype=float)
    n_C = np.array(n_C, dtype=float)

    # Safety: skip zero vectors
    if np.sum(n_A) < 1e-12 or np.sum(n_B) < 1e-12 or np.sum(n_C) < 1e-12:
        return None

    try:
        u_coeff, v_coeff, k_coeff = metric_VectorDecomposition_onlyPositive(n_A, n_B, n_C)
    except (np.linalg.LinAlgError, FloatingPointError, ZeroDivisionError):
        return None

    # Guard NaN/Inf
    if not (np.isfinite(u_coeff) and np.isfinite(v_coeff) and np.isfinite(k_coeff)):
        return None

    x, y = calculate_assymetricity(u_coeff, v_coeff, k_coeff)
    if not (np.isfinite(x) and np.isfinite(y)):
        return None

    cls = characterize_case(x, y)
    return cls, x, y, u_coeff, v_coeff, k_coeff


print("Classifying events...")
results = []

for _, event in coal_df.iterrows():
    sid = event["SampleIDX"]
    if sid in exception_list:
        continue

    sid_sub1 = event["SampleIDX_Sub1"]
    sid_sub2 = event["SampleIDX_Sub2"]
    medium = event["Medium"]

    # Look up abundance vectors
    n_C = seq_lookup.get(sid)
    n_A = seq_lookup.get(sid_sub1)
    n_B = seq_lookup.get(sid_sub2)

    if n_C is None or n_A is None or n_B is None:
        continue

    # ------ Observed classification ------
    obs = classify_vector(n_C, n_A, n_B)
    if obs is None:
        continue
    obs_cls, obs_x, obs_y, obs_u, obs_v, obs_k = obs

    # ------ Additive null ------
    n_C_null_raw = np.array(n_A, dtype=float) + np.array(n_B, dtype=float)
    s = np.sum(n_C_null_raw)
    if s < 1e-12:
        continue
    n_C_null = n_C_null_raw / s  # normalize to relative abundance

    null = classify_vector(n_C_null, n_A, n_B)
    if null is None:
        continue
    null_cls, null_x, null_y, null_u, null_v, null_k = null

    results.append({
        "SampleIDX": sid,
        "Medium": medium,
        "obs_class": obs_cls,
        "null_class": null_cls,
        "obs_x": obs_x,
        "obs_y": obs_y,
        "null_x": null_x,
        "null_y": null_y,
        "obs_PDI": obs_y,   # PDI ~ y-coordinate (asymmetricity)
        "null_PDI": null_y,
        "obs_u": obs_u,
        "obs_v": obs_v,
        "obs_k": obs_k,
        "null_u": null_u,
        "null_v": null_v,
        "null_k": null_k,
    })

df = pd.DataFrame(results)
print(f"\nTotal events analysed: {len(df)}")

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("SUMMARY: Per-event additive null model")
print("=" * 60)

# Overall contingency
print("\nOverall contingency (rows=observed, cols=null):")
ct = pd.crosstab(
    df["obs_class"].map(CLASS_NAMES),
    df["null_class"].map(CLASS_NAMES),
    margins=True,
)
ct = ct.reindex(index=["Dominance", "Mixture", "Restructuring", "All"],
                columns=["Dominance", "Mixture", "Restructuring", "All"],
                fill_value=0)
print(ct.to_string())

print(f"\nEvents where observed=Dominance but null=Mixture: "
      f"{len(df[(df['obs_class']==0) & (df['null_class']==1)])}")
print(f"Events where observed=Dominance but null=Restructuring: "
      f"{len(df[(df['obs_class']==0) & (df['null_class']==2)])}")
print(f"Events where both observed & null = Dominance: "
      f"{len(df[(df['obs_class']==0) & (df['null_class']==0)])}")

# Agreement rate
agree = (df["obs_class"] == df["null_class"]).sum()
print(f"\nAgreement rate (obs==null): {agree}/{len(df)} = {agree/len(df):.1%}")

# Per-medium summary
for med_code in ["L", "M", "H"]:
    sub = df[df["Medium"] == med_code]
    if len(sub) == 0:
        continue
    med_label = MEDIUM_MAP[med_code]
    print(f"\n--- {med_label} (n={len(sub)}) ---")
    for cls_i, cls_name in CLASS_NAMES.items():
        obs_n = (sub["obs_class"] == cls_i).sum()
        null_n = (sub["null_class"] == cls_i).sum()
        print(f"  {cls_name:15s}  obs: {obs_n:3d} ({obs_n/len(sub):.0%})   null: {null_n:3d} ({null_n/len(sub):.0%})")

# PDI shift
delta_PDI = df["obs_PDI"] - df["null_PDI"]
print(f"\nPDI_obs - PDI_null:  mean={delta_PDI.mean():.4f}  "
      f"median={delta_PDI.median():.4f}")
tstat, pval = stats.ttest_rel(df["obs_PDI"], df["null_PDI"])
print(f"Paired t-test: t={tstat:.3f}, p={pval:.2e}")
wstat, wpval = stats.wilcoxon(df["obs_PDI"], df["null_PDI"])
print(f"Wilcoxon signed-rank: W={wstat:.1f}, p={wpval:.2e}")


# ---------------------------------------------------------------------------
# Helper: save figure in 3 formats
# ---------------------------------------------------------------------------
def save_fig(fig, basename):
    for ext in ["svg", "pdf", "png"]:
        path = os.path.join(SAVE_DIR, f"{basename}.{ext}")
        dpi_val = 300 if ext == "png" else None
        fig.savefig(path, bbox_inches="tight", dpi=dpi_val, transparent=(ext != "png"))
    print(f"  Saved {basename}.{{svg,pdf,png}}")


# ===================================================================
# FIGURE 1: Paired comparison — stacked bar (null vs observed) per event
# ===================================================================
print("\nGenerating Figure 1: Paired classification comparison...")

# We'll show a grouped bar chart: for each class, count obs vs null
fig1, ax1 = plt.subplots(figsize=(3.2, 2.4))

class_order = [0, 1, 2]
class_labels = [CLASS_NAMES[c] for c in class_order]
obs_counts = [int((df["obs_class"] == c).sum()) for c in class_order]
null_counts = [int((df["null_class"] == c).sum()) for c in class_order]

x_pos = np.arange(len(class_order))
width = 0.35

bars_obs = ax1.bar(x_pos - width / 2, obs_counts, width,
                   color=[CLASS_COLORS[c] for c in class_order],
                   edgecolor="k", linewidth=0.4, label="Observed")
bars_null = ax1.bar(x_pos + width / 2, null_counts, width,
                    color=[CLASS_COLORS[c] for c in class_order],
                    edgecolor="k", linewidth=0.4, alpha=0.4, label="Null (additive)")

ax1.set_xticks(x_pos)
ax1.set_xticklabels(class_labels, fontsize=8)
ax1.set_ylabel("Number of events", fontsize=8)
from matplotlib.patches import Patch
legend_handles = [
    Patch(facecolor="0.5", edgecolor="k", linewidth=0.4, label="Observed"),
    Patch(facecolor="0.5", edgecolor="k", linewidth=0.4, alpha=0.4,
          label="Null (additive)"),
]
ax1.legend(handles=legend_handles, fontsize=6, frameon=False)
ax1.set_title("Observed vs. additive null classification", fontsize=8)
sns.despine(ax=ax1)

# Add count labels
for bar_group in [bars_obs, bars_null]:
    for bar in bar_group:
        h = bar.get_height()
        if h > 0:
            ax1.text(bar.get_x() + bar.get_width() / 2, h + 1,
                     str(int(h)), ha="center", va="bottom", fontsize=6)

fig1.tight_layout()
save_fig(fig1, "fig1_paired_classification")
plt.close(fig1)


# ===================================================================
# FIGURE 2: Histogram of PDI_observed - PDI_null
# ===================================================================
print("Generating Figure 2: Delta-PDI histogram...")

fig2, ax2 = plt.subplots(figsize=(3.2, 2.4))

delta = df["obs_PDI"] - df["null_PDI"]
ax2.hist(delta, bins=30, color="#888888", edgecolor="k", linewidth=0.3, alpha=0.8)
ax2.axvline(0, color="k", linestyle="--", linewidth=0.8)
ax2.axvline(delta.mean(), color="#D32F2F", linestyle="-", linewidth=1.0,
            label=f"mean = {delta.mean():.3f}")
ax2.axvline(delta.median(), color="#1976D2", linestyle="-", linewidth=1.0,
            label=f"median = {delta.median():.3f}")

ax2.set_xlabel("PDI$_{obs}$ - PDI$_{null}$", fontsize=8)
ax2.set_ylabel("Count", fontsize=8)
ax2.set_title("Shift in Parent Dominance Index", fontsize=8)
ax2.legend(fontsize=6, frameon=False)
sns.despine(ax=ax2)
fig2.tight_layout()
save_fig(fig2, "fig2_delta_PDI_histogram")
plt.close(fig2)


# ===================================================================
# FIGURE 3: Contingency-table heatmap
# ===================================================================
print("Generating Figure 3: Contingency heatmap...")

fig3, ax3 = plt.subplots(figsize=(3.0, 2.6))

ct_raw = pd.crosstab(
    df["obs_class"].map(CLASS_NAMES),
    df["null_class"].map(CLASS_NAMES),
)
# Ensure order
ordered = ["Dominance", "Mixture", "Restructuring"]
ct_plot = ct_raw.reindex(index=ordered, columns=ordered, fill_value=0)

sns.heatmap(
    ct_plot,
    annot=True, fmt="d",
    cmap="Greys",
    linewidths=0.5,
    linecolor="white",
    ax=ax3,
    cbar_kws={"shrink": 0.6, "label": "Count"},
    annot_kws={"size": 9},
)
ax3.set_ylabel("Observed class", fontsize=8)
ax3.set_xlabel("Null class", fontsize=8)
ax3.set_title("Observed vs. null contingency", fontsize=8)
ax3.tick_params(axis="both", labelsize=8)

fig3.tight_layout()
save_fig(fig3, "fig3_contingency_heatmap")
plt.close(fig3)


# ===================================================================
# FIGURE 4: Per-medium breakdown panels
# ===================================================================
print("Generating Figure 4: Per-medium panels...")

medium_codes = ["L", "M", "H"]
medium_labels_short = {"L": "LN (Nutr$-$)", "M": "MN (Base)", "H": "HN (Nutr$+$)"}

fig4, axes4 = plt.subplots(1, 3, figsize=(7.5, 2.5), sharey=True)

for ax_idx, med_code in enumerate(medium_codes):
    ax = axes4[ax_idx]
    sub = df[df["Medium"] == med_code]

    obs_c = [int((sub["obs_class"] == c).sum()) for c in class_order]
    null_c = [int((sub["null_class"] == c).sum()) for c in class_order]

    bars_o = ax.bar(x_pos - width / 2, obs_c, width,
                    color=[CLASS_COLORS[c] for c in class_order],
                    edgecolor="k", linewidth=0.4)
    bars_n = ax.bar(x_pos + width / 2, null_c, width,
                    color=[CLASS_COLORS[c] for c in class_order],
                    edgecolor="k", linewidth=0.4, alpha=0.4)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(["Dom", "Mix", "Rest"], fontsize=6)
    med_label = medium_labels_short[med_code]
    med_color = get_medium_color(med_code)
    ax.set_title(med_label, fontsize=8, color=med_color, fontweight="bold")
    sns.despine(ax=ax)

    # Annotate counts
    for bar_group in [bars_o, bars_n]:
        for bar in bar_group:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                        str(int(h)), ha="center", va="bottom", fontsize=5)

axes4[0].set_ylabel("Number of events", fontsize=8)

leg_handles = [
    Patch(facecolor="#888888", edgecolor="k", linewidth=0.4, label="Observed"),
    Patch(facecolor="#888888", edgecolor="k", linewidth=0.4, alpha=0.4, label="Null"),
]
axes4[2].legend(handles=leg_handles, fontsize=6, frameon=False, loc="upper right")

# Panel labels
for i, ax in enumerate(axes4):
    ax.text(-0.15, 1.05, f"({chr(ord('a') + i)})", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="right")

fig4.tight_layout()
save_fig(fig4, "fig4_per_medium_panels")
plt.close(fig4)


# ===================================================================
# FIGURE 5: Per-medium contingency heatmaps
# ===================================================================
print("Generating Figure 5: Per-medium contingency heatmaps...")

fig5, axes5 = plt.subplots(1, 3, figsize=(8.0, 2.6))

for ax_idx, med_code in enumerate(medium_codes):
    ax = axes5[ax_idx]
    sub = df[df["Medium"] == med_code]

    ct_m = pd.crosstab(
        sub["obs_class"].map(CLASS_NAMES),
        sub["null_class"].map(CLASS_NAMES),
    )
    ct_m = ct_m.reindex(index=ordered, columns=ordered, fill_value=0)

    sns.heatmap(
        ct_m, annot=True, fmt="d", cmap="Greys",
        linewidths=0.5, linecolor="white",
        ax=ax, cbar=False,
        annot_kws={"size": 9},
    )
    med_label = medium_labels_short[med_code]
    med_color = get_medium_color(med_code)
    ax.set_title(med_label, fontsize=8, color=med_color, fontweight="bold")
    ax.set_xlabel("Null class" if ax_idx == 1 else "", fontsize=8)
    ax.set_ylabel("Observed class" if ax_idx == 0 else "", fontsize=8)
    ax.tick_params(axis="both", labelsize=6)
    if ax_idx > 0:
        ax.set_yticklabels([])

# Panel labels
for i, ax in enumerate(axes5):
    ax.text(-0.15, 1.05, f"({chr(ord('a') + i)})", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="right")

fig5.tight_layout()
save_fig(fig5, "fig5_per_medium_contingency")
plt.close(fig5)


# ===================================================================
# FIGURE 6: Scatter — observed (x,y) vs null (x,y) in asymmetricity space
# ===================================================================
print("Generating Figure 6: Asymmetricity space comparison...")

fig6, (ax6a, ax6b) = plt.subplots(1, 2, figsize=(5.5, 2.5), sharey=True, sharex=True)

for ax_i, (ax, col_prefix, title) in enumerate([
    (ax6a, "obs", "Observed"),
    (ax6b, "null", "Null (additive)"),
]):
    for _, row in df.iterrows():
        cx = row[f"{col_prefix}_x"]
        cy = row[f"{col_prefix}_y"]
        cls = int(row[f"{col_prefix}_class"])
        ax.scatter(cx**2, cy, color=CLASS_COLORS[cls], s=8, alpha=0.5,
                   edgecolors="none")

    # Draw classification boundaries
    ax.axvline(0.5, color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.axhline(0.5, color="k", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("$x^2$ (retention)", fontsize=8)
    ax.set_title(title, fontsize=8)
    sns.despine(ax=ax)

ax6a.set_ylabel("$y$ (asymmetry)", fontsize=8)

# Panel labels
for i, ax in enumerate([ax6a, ax6b]):
    ax.text(-0.15, 1.05, f"({chr(ord('a') + i)})", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="right")

fig6.tight_layout()
save_fig(fig6, "fig6_asymmetricity_space")
plt.close(fig6)


# ===================================================================
# COMBINED RESPONSE FIGURE: Base-medium raw-count null comparison
# ===================================================================
print("Generating combined response figure: additive null comparison...")

base_rows = []
for _, event in coal_df.iterrows():
    sid = event["SampleIDX"]
    if sid in exception_list or event["Medium"] != "M":
        continue
    sid_sub1 = event["SampleIDX_Sub1"]
    sid_sub2 = event["SampleIDX_Sub2"]
    n_C_raw = raw_lookup.get(sid)
    n_A_raw = raw_lookup.get(sid_sub1)
    n_B_raw = raw_lookup.get(sid_sub2)
    if n_C_raw is None or n_A_raw is None or n_B_raw is None:
        continue

    obs = classify_vector(n_C_raw, n_A_raw, n_B_raw)
    null = classify_vector(n_A_raw + n_B_raw, n_A_raw, n_B_raw)
    if obs is None or null is None:
        continue
    obs_cls, obs_x, obs_y, obs_u, obs_v, obs_k = obs
    null_cls, null_x, null_y, null_u, null_v, null_k = null
    base_rows.append(
        {
            "SampleIDX": sid,
            "obs_class": obs_cls,
            "null_class": null_cls,
            "obs_u": obs_u,
            "obs_v": obs_v,
            "obs_x": obs_x,
            "obs_y": obs_y,
            "null_u": null_u,
            "null_v": null_v,
            "null_x": null_x,
            "null_y": null_y,
            "delta_y": obs_y - null_y,
        }
    )

base_raw_df = pd.DataFrame(base_rows)
base_raw_df.to_csv(os.path.join(SAVE_DIR, "base_raw_count_additive_null_events.csv"), index=False)
base_delta = base_raw_df["delta_y"]
base_wstat, base_wpval = stats.wilcoxon(
    base_raw_df["obs_y"],
    base_raw_df["null_y"],
    alternative="greater",
)
representative_ids = (
    base_raw_df.sample(5, random_state=23)["SampleIDX"].sort_values().tolist()
)
representative = base_raw_df[base_raw_df["SampleIDX"].isin(representative_ids)]
print("\nBase raw-count additive null:")
print(pd.crosstab(base_raw_df["obs_class"].map(CLASS_NAMES), base_raw_df["null_class"].map(CLASS_NAMES)).to_string())
print(
    f"  n={len(base_raw_df)}, median delta y={base_delta.median():.4f}, "
    f"delta>0={(base_delta > 0).sum()}/{len(base_raw_df)}, "
    f"one-sided Wilcoxon p={base_wpval:.2e}"
)
print(f"  Representative paired events: {', '.join(representative_ids)}")

fig_combined = plt.figure(figsize=(8.4, 2.9))
gs = fig_combined.add_gridspec(1, 3, width_ratios=[1.08, 0.75, 0.95], wspace=0.42)
axc0 = fig_combined.add_subplot(gs[0, 0])
axc1 = fig_combined.add_subplot(gs[0, 1])
axc2 = fig_combined.add_subplot(gs[0, 2])

def draw_similarity_boundaries(ax):
    theta = np.linspace(0, np.pi / 2, 300)
    ax.plot(np.cos(theta), np.sin(theta), color="black", linewidth=1.5)
    manuscript_r = np.sqrt(0.5)
    ax.plot(
        manuscript_r * np.cos(theta),
        manuscript_r * np.sin(theta),
        color="black",
        linewidth=0.8,
        linestyle=(0, (4, 4)),
    )
    for ang in (np.pi / 8, 3 * np.pi / 8):
        ax.plot(
            [0, np.cos(ang)],
            [0, np.sin(ang)],
            color="black",
            linewidth=0.8,
            linestyle=(0, (4, 4)),
        )
    ax.plot([0, 0], [0, 1], color="black", linewidth=1.0)
    ax.plot([0, 1], [0, 0], color="black", linewidth=1.0)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.tick_params(length=0, pad=1)
    for spine in ax.spines.values():
        spine.set_visible(False)


# Panel A: similarity map with count-vector null and observed outcomes.
draw_similarity_boundaries(axc0)
for _, row in representative.iterrows():
    axc0.annotate(
        "",
        xy=(row["obs_u"], row["obs_v"]),
        xytext=(row["null_u"], row["null_v"]),
        arrowprops={
            "arrowstyle": "->",
            "color": "#cc3d3d",
            "linewidth": 1.2,
            "alpha": 0.95,
            "shrinkA": 2.0,
            "shrinkB": 2.0,
        },
        zorder=4,
    )
axc0.scatter(
    base_raw_df["null_u"],
    base_raw_df["null_v"],
    s=15,
    facecolors="white",
    edgecolors="#d97706",
    linewidths=0.7,
    alpha=0.85,
    label="additive null",
    zorder=3,
)
axc0.scatter(
    base_raw_df["obs_u"],
    base_raw_df["obs_v"],
    s=16,
    color="#2563a8",
    alpha=0.78,
    linewidths=0,
    label="observed",
    zorder=3,
)
axc0.set_xlabel("similarity to parent A", fontsize=8)
axc0.set_ylabel("similarity to parent B", fontsize=8)
axc0.set_title("Base raw-count similarity map", fontsize=8)
axc0.legend(
    frameon=False,
    fontsize=5.8,
    loc="lower left",
    bbox_to_anchor=(0.02, 0.10),
    handlelength=1.0,
    labelspacing=0.25,
)

# Panel B: per-event shift in direction-independent parental asymmetry.
x_null = np.zeros(len(base_raw_df))
x_obs = np.ones(len(base_raw_df))
for _, row in base_raw_df.iterrows():
    axc1.annotate(
        "",
        xy=(1, row["obs_y"]),
        xytext=(0, row["null_y"]),
        arrowprops={
            "arrowstyle": "->",
            "color": "0.72",
            "linewidth": 0.55,
            "alpha": 0.48,
            "shrinkA": 1.5,
            "shrinkB": 1.5,
        },
        zorder=1,
    )
axc1.scatter(x_null, base_raw_df["null_y"], s=12, facecolors="white", edgecolors="#d97706", linewidths=0.6, zorder=3)
axc1.scatter(x_obs, base_raw_df["obs_y"], s=13, color="#2563a8", alpha=0.75, linewidths=0, zorder=3)
axc1.axhline(0.5, color="black", linestyle=(0, (4, 4)), linewidth=0.8, alpha=0.6)
axc1.set_xlim(-0.25, 1.25)
axc1.set_ylim(-0.03, 1.03)
axc1.set_xticks([0, 1])
axc1.set_xticklabels(["Null", "Exp."])
axc1.set_ylabel("asymmetry $y = |2\\mathrm{PDI}-1|$", fontsize=8)
axc1.set_title("Paired asymmetry", fontsize=8)
axc1.text(
    0.04,
    0.96,
    f"{(base_delta > 0).sum()}/{len(base_raw_df)} higher asymmetry\nWilcoxon $p=5.9\\times10^{{-14}}$",
    transform=axc1.transAxes,
    ha="left",
    va="top",
    fontsize=6,
    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
)
sns.despine(ax=axc1)

# Panel C: class transition counts from raw-count null to experiment.
ordered_classes = ["Dominance", "Mixture", "Restructuring"]
short_classes = ["Dom", "Mix", "Rest"]
transition_counts = pd.crosstab(
    base_raw_df["obs_class"].map(CLASS_NAMES),
    base_raw_df["null_class"].map(CLASS_NAMES),
).reindex(index=ordered_classes, columns=ordered_classes, fill_value=0)
sns.heatmap(
    transition_counts,
    annot=True,
    fmt="d",
    cmap="Greys",
    cbar=False,
    linewidths=0.5,
    linecolor="white",
    ax=axc2,
    annot_kws={"size": 9},
    vmin=0,
    vmax=transition_counts.to_numpy().max(),
)
axc2.set_xticklabels(short_classes, rotation=0, fontsize=8)
axc2.set_yticklabels(short_classes, rotation=0, fontsize=8)
axc2.set_xlabel("Additive-null class", fontsize=8)
axc2.set_ylabel("Observed class", fontsize=8)
axc2.set_title("Class transitions", fontsize=8)

for label, ax in zip(["A", "B", "C"], [axc0, axc1, axc2]):
    ax.text(
        -0.16,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="bottom",
        ha="right",
    )

fig_combined.savefig(
    os.path.join(SAVE_DIR, "Fig_R3_2_additive_null_comparison.pdf"),
    bbox_inches="tight",
)
fig_combined.savefig(
    os.path.join(SAVE_DIR, "Fig_R3_2_additive_null_comparison.png"),
    bbox_inches="tight",
    dpi=300,
)
fig_combined.savefig(
    os.path.join(SAVE_DIR, "Fig_R3_2_additive_null_comparison.svg"),
    bbox_inches="tight",
)
print("  Saved Fig_R3_2_additive_null_comparison.{svg,pdf,png}")
plt.close(fig_combined)


# ===================================================================
# FIGURE 7: Per-medium delta-PDI histograms
# ===================================================================
print("Generating Figure 7: Per-medium delta-PDI histograms...")

fig7, axes7 = plt.subplots(1, 3, figsize=(7.5, 2.2), sharey=True, sharex=True)

for ax_idx, med_code in enumerate(medium_codes):
    ax = axes7[ax_idx]
    sub = df[df["Medium"] == med_code]
    delta_m = sub["obs_PDI"] - sub["null_PDI"]

    med_color = get_medium_color(med_code)
    ax.hist(delta_m, bins=20, color=med_color, edgecolor="k", linewidth=0.3, alpha=0.7)
    ax.axvline(0, color="k", linestyle="--", linewidth=0.6)
    ax.axvline(delta_m.mean(), color=med_color, linestyle="-", linewidth=1.2)

    med_label = medium_labels_short[med_code]
    ax.set_title(f"{med_label} (n={len(sub)})", fontsize=8, color=med_color,
                 fontweight="bold")
    ax.set_xlabel("$\\Delta$PDI" if ax_idx == 1 else "", fontsize=8)

    # Stats annotation
    if len(sub) > 5:
        t, p = stats.ttest_rel(sub["obs_PDI"], sub["null_PDI"])
        sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "n.s."))
        ax.text(0.95, 0.92, f"mean={delta_m.mean():.3f}\n{sig} (p={p:.1e})",
                transform=ax.transAxes, fontsize=5, ha="right", va="top")
    sns.despine(ax=ax)

axes7[0].set_ylabel("Count", fontsize=8)

# Panel labels
for i, ax in enumerate(axes7):
    ax.text(-0.15, 1.05, f"({chr(ord('a') + i)})", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="right")

fig7.tight_layout()
save_fig(fig7, "fig7_per_medium_delta_PDI")
plt.close(fig7)


# ===================================================================
# FIGURE 8: Transition diagram (alluvial-style)
# ===================================================================
print("Generating Figure 8: Transition diagram...")

fig8, ax8 = plt.subplots(figsize=(3.5, 3.0))

# Count transitions
transitions = {}
for _, row in df.iterrows():
    key = (int(row["null_class"]), int(row["obs_class"]))
    transitions[key] = transitions.get(key, 0) + 1

# Position classes vertically
y_positions = {0: 2.5, 1: 1.5, 2: 0.5}
x_null = 0.2
x_obs = 0.8

# Draw class boxes
for cls in [0, 1, 2]:
    null_n = int((df["null_class"] == cls).sum())
    obs_n = int((df["obs_class"] == cls).sum())
    yp = y_positions[cls]

    # Null side
    ax8.add_patch(plt.Rectangle((x_null - 0.08, yp - 0.15), 0.16, 0.30,
                                facecolor=CLASS_COLORS[cls], edgecolor="k",
                                linewidth=0.5, alpha=0.5))
    ax8.text(x_null, yp, f"{null_n}", ha="center", va="center", fontsize=8,
             fontweight="bold")

    # Obs side
    ax8.add_patch(plt.Rectangle((x_obs - 0.08, yp - 0.15), 0.16, 0.30,
                                facecolor=CLASS_COLORS[cls], edgecolor="k",
                                linewidth=0.5))
    ax8.text(x_obs, yp, f"{obs_n}", ha="center", va="center", fontsize=8,
             fontweight="bold")

    # Labels
    ax8.text(x_null, yp + 0.22, CLASS_NAMES[cls], ha="center", va="bottom",
             fontsize=6, color=CLASS_COLORS[cls])
    ax8.text(x_obs, yp + 0.22, CLASS_NAMES[cls], ha="center", va="bottom",
             fontsize=6, color=CLASS_COLORS[cls])

# Draw arrows for transitions
max_count = max(transitions.values()) if transitions else 1
for (src, dst), count in transitions.items():
    lw = 0.5 + 4.0 * (count / max_count)
    alpha = 0.2 + 0.6 * (count / max_count)
    y_src = y_positions[src]
    y_dst = y_positions[dst]

    ax8.annotate(
        "", xy=(x_obs - 0.09, y_dst), xytext=(x_null + 0.09, y_src),
        arrowprops=dict(
            arrowstyle="->,head_width=0.15,head_length=0.08",
            color=CLASS_COLORS[dst],
            linewidth=lw,
            alpha=alpha,
            connectionstyle="arc3,rad=0.1" if src != dst else "arc3,rad=0.0",
        ),
    )
    # Count label on arrow
    mid_x = (x_null + 0.09 + x_obs - 0.09) / 2
    mid_y = (y_src + y_dst) / 2
    offset_y = 0.08 if src != dst else -0.12
    if count >= 3:
        ax8.text(mid_x, mid_y + offset_y, str(count), ha="center", va="center",
                 fontsize=5, color=CLASS_COLORS[dst], alpha=0.9)

ax8.text(x_null, 3.1, "Null\n(additive)", ha="center", va="bottom", fontsize=8,
         fontweight="bold")
ax8.text(x_obs, 3.1, "Observed", ha="center", va="bottom", fontsize=8,
         fontweight="bold")

ax8.set_xlim(-0.05, 1.05)
ax8.set_ylim(-0.1, 3.5)
ax8.axis("off")
fig8.tight_layout()
save_fig(fig8, "fig8_transition_diagram")
plt.close(fig8)


# ===================================================================
# FIGURE 9: Mixing-ratio sweep (α * n_A + (1-α) * n_B, α ∈ [0,1])
# ===================================================================
# Addresses the recommendation to add a mixing-ratio sweep:
# The additive null (α=0.5) always classifies as Mixture by mathematical
# necessity (equal parental contributions → u=v → Mixture). Showing the
# full α sweep demonstrates which ratios produce Dominance vs Mixture.
print("Generating Figure 9: Mixing-ratio sweep...")

alphas = np.linspace(0, 1, 21)  # 0.00, 0.05, ..., 1.00
sweep_results = {med_code: {cls: [] for cls in [0, 1, 2]} for med_code in ["L", "M", "H"]}
sweep_counts = {med_code: [] for med_code in ["L", "M", "H"]}

for alpha in alphas:
    per_medium = {med_code: {cls: 0 for cls in [0, 1, 2]} for med_code in ["L", "M", "H"]}
    per_medium_n = {med_code: 0 for med_code in ["L", "M", "H"]}

    for _, row in coal_df.iterrows():
        sid = row["SampleIDX"]
        if sid in exception_list:
            continue
        sid_sub1 = row["SampleIDX_Sub1"]
        sid_sub2 = row["SampleIDX_Sub2"]
        medium = row["Medium"]

        n_C_orig = seq_lookup.get(sid)
        n_A = seq_lookup.get(sid_sub1)
        n_B = seq_lookup.get(sid_sub2)
        if n_C_orig is None or n_A is None or n_B is None:
            continue

        # Swept mixture: alpha fraction from A, (1-alpha) from B
        n_mix_alpha = alpha * np.array(n_A) + (1 - alpha) * np.array(n_B)
        s = np.sum(n_mix_alpha)
        if s < 1e-12:
            continue
        n_mix_alpha = n_mix_alpha / s

        result = classify_vector(n_mix_alpha, n_A, n_B)
        if result is None:
            continue
        cls = result[0]
        per_medium[medium][cls] += 1
        per_medium_n[medium] += 1

    for med_code in ["L", "M", "H"]:
        n_total = per_medium_n[med_code]
        for cls in [0, 1, 2]:
            frac = per_medium[med_code][cls] / n_total if n_total > 0 else 0
            sweep_results[med_code][cls].append(frac)
        sweep_counts[med_code].append(n_total)

fig9, axes9 = plt.subplots(1, 3, figsize=(7.5, 2.5), sharey=True)

for ax_idx, med_code in enumerate(["L", "M", "H"]):
    ax = axes9[ax_idx]
    med_color = get_medium_color(med_code)

    for cls, cls_name, cls_color in zip([0, 1, 2],
                                         ["Dominance", "Mixture", "Restructuring"],
                                         [CLASS_COLORS[0], CLASS_COLORS[1], CLASS_COLORS[2]]):
        ax.plot(alphas, sweep_results[med_code][cls], '-',
                color=cls_color, linewidth=1.5, label=cls_name)

    ax.axvline(0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.6,
               label='Additive null (α=0.5)')
    ax.set_xlabel('Mixing ratio α\n(α·nA + (1−α)·nB)', fontsize=7)
    ax.set_title(f'{medium_labels_short[med_code]}', fontsize=8, color=med_color,
                 fontweight='bold')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    if ax_idx == 0:
        ax.set_ylabel('Fraction of events', fontsize=8)
    ax.legend(fontsize=5, frameon=False, loc='upper center')
    sns.despine(ax=ax)

# Panel labels
for i, ax in enumerate(axes9):
    ax.text(-0.15, 1.05, f"({chr(ord('a') + i)})", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="right")

fig9.tight_layout()
save_fig(fig9, "fig9_mixing_ratio_sweep")
plt.close(fig9)


# ===================================================================
# Print key takeaway
# ===================================================================
print("\n" + "=" * 60)
print("KEY TAKEAWAY")
print("=" * 60)
dom_obs = (df["obs_class"] == 0).sum()
dom_null = (df["null_class"] == 0).sum()
dom_both = ((df["obs_class"] == 0) & (df["null_class"] == 0)).sum()
dom_obs_only = ((df["obs_class"] == 0) & (df["null_class"] != 0)).sum()
print(f"Dominance events:  Observed={dom_obs},  Null={dom_null}")
print(f"  Both obs & null = Dominance (potential artifact): {dom_both}")
print(f"  Obs=Dominance but null!=Dominance (true selection): {dom_obs_only}")
if dom_obs > 0:
    print(f"  Fraction of observed Dominance that is 'true': "
          f"{dom_obs_only}/{dom_obs} = {dom_obs_only / dom_obs:.1%}")

# Per medium breakdown
for med_code in medium_codes:
    sub = df[df["Medium"] == med_code]
    med_label = MEDIUM_MAP[med_code]
    d_obs = (sub["obs_class"] == 0).sum()
    d_null = (sub["null_class"] == 0).sum()
    d_both = ((sub["obs_class"] == 0) & (sub["null_class"] == 0)).sum()
    d_true = ((sub["obs_class"] == 0) & (sub["null_class"] != 0)).sum()
    print(f"\n  {med_label}: obs_Dom={d_obs}, null_Dom={d_null}, "
          f"artifact={d_both}, true_selection={d_true}")

print("\nDone.")
