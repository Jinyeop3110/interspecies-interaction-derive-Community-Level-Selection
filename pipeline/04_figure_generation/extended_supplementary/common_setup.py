"""Shared harness for the Extended Data and Supplementary analysis scripts.

The scripts in this directory were written as working analysis code and each
begins with ``from common_setup import *``.  This module provides what they
expect: the plotting defaults, the loaded data tables, the community-lookup
helpers, and the similarity-decomposition functions.

It is a port of ``common_setup.py`` from the working tree, changed in three
ways so that the scripts run outside the machine they were written on:

* Data is loaded through :mod:`coalescence.io`, so ``COALESCENCE_DATA`` works
  and the hardcoded ``../../Analyzed/`` paths are gone.
* Nothing is created or written at import time.
* ``exception_list`` is :data:`coalescence.io.EXCLUDED_SAMPLES` rather than a
  separate inline copy, so there is one definition of the quality-control
  exclusions in the repository.

The helper signatures are unchanged, so the scripts do not need editing.

For new work prefer the library in ``src/coalescence`` directly; this module
exists to keep the existing analysis scripts runnable.
"""

import os
import random
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from coalescence import io  # noqa: E402
from coalescence.decomposition import (  # noqa: E402
    decompose as _decompose, normalize,
)

# --------------------------------------------------------------------------
# Plotting defaults
# --------------------------------------------------------------------------

sns.set_style("ticks")

#: Figure dimensions in the scripts are given in mm.
mm = 0.1 / 2.54

#: Medium colours, in the order Nutr-, Base, Nutr+.
colors = ["#A7216A", "#802000", "#E24912"]

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

my_cmap = ListedColormap(
    ["#00429d", "#4771b2", "#73a2c6", "#a5d5d8", "#ffffe0",
     "#ffbcaf", "#f4777f", "#cf3759", "#93003a"],
    name="my_cmap",
)
my_cmap_r = my_cmap.reversed()
for _cmap in (my_cmap, my_cmap_r):
    if _cmap.name not in mpl.colormaps:
        mpl.colormaps.register(cmap=_cmap)

np.random.seed(1)

# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------

Coalescence_data = pd.concat([
    io.load_coalescence_events("synthetic"),
    io.load_coalescence_events("natural"),
])
Communities_data = pd.concat([
    io.load_communities("synthetic"),
    io.load_communities("natural"),
])
Coalescence_data_sim = io._read("processed_CoalescenceEvent_simulation_uniform.xlsx")

Coalescence_recipe = io.load_coalescence_recipe("synthetic")
Coalescence_natural_recipe = io.load_coalescence_recipe("natural")
Metadata = io.load_metadata()
Processed_sequences_synthetic = io.load_sequences("synthetic").reset_index()
Processed_sequences_natural = io.load_sequences("natural").reset_index()

#: Quality-control exclusions; see :data:`coalescence.io.EXCLUSION_REASONS`.
exception_list = sorted(io.EXCLUDED_SAMPLES)

# --------------------------------------------------------------------------
# Community lookup
# --------------------------------------------------------------------------

#: Community index ranges per richness class, as ``(parental, coalesced)``.
_RICHNESS_RANGES = {
    6: ((0, 9), (0, 14)),
    12: ((9, 18), (14, 41)),
    24: ((18, 30), (41, 47)),
}


def getIDX(data, IDX):
    return np.array(np.flatnonzero([x == IDX for x in data["SampleIDX"]]))


def getAbundance(Processed_sequences, IDX):
    match = np.flatnonzero([x == IDX for x in Processed_sequences["SampleIDX"]])
    if not any(match):
        return
    return Processed_sequences.iloc[match].values.tolist()[0][1:]


def getAbundance_fromList(Processed_sequences, IDX_list):
    rows = np.squeeze([np.where(Processed_sequences["SampleIDX"] == x) for x in IDX_list])
    return (np.array(Processed_sequences.iloc[rows].values)[:, 1:]).astype(float)


def Community_PermutateList_simulation(I, S, sample_num=500):
    """Row indices into the simulation table for one (interaction, richness) cell."""
    return (3 * (I - 1) + (S - 1)) * sample_num + np.arange(0, sample_num)


def CommunityPermutate(Timepoint, CommunityOrigin, Medium, CoalescenceType):
    idx = np.where(
        (Metadata["Timepoint"] == Timepoint)
        & (Metadata["CommunityOrigin"] == CommunityOrigin)
        & (Metadata["Medium"] == Medium)
        & (Metadata["CoalescenceType"] == CoalescenceType)
    )[0]
    return np.concatenate(Metadata["SampleIDX"][idx])


def _richness_mask(base, CoalescenceType, species_pool_num):
    """Restrict a selection to one richness class."""
    if species_pool_num == 0:
        return base
    community_index = np.array([int(x) for x in Metadata["CommunityIDX"]])
    parental, coalesced = _RICHNESS_RANGES[species_pool_num]
    low, high = parental if CoalescenceType == "S" else coalesced
    return base & (community_index > low) & (community_index <= high)


def CommunityPermutate_withSpeciesPoolsize(Timepoint, CommunityOrigin, Medium,
                                           CoalescenceType, species_pool_num):
    idx = (
        (Metadata["Timepoint"] == Timepoint)
        & (Metadata["CommunityOrigin"] == CommunityOrigin)
        & (Metadata["Medium"] == Medium)
        & (Metadata["CoalescenceType"] == CoalescenceType)
    )
    idx = _richness_mask(idx, CoalescenceType, species_pool_num)
    return list(set(Metadata["SampleIDX"][idx]) - set(exception_list))


def Community_PermutateList(Timepoint, CommunityOrigin, Medium, CoalescenceType,
                            species_pool_num=0, Replicate=-1):
    """Sample IDs matching the given experimental factors, exclusions applied."""
    idx = (
        (Metadata["Timepoint"] == Timepoint)
        & (Metadata["CommunityOrigin"] == CommunityOrigin)
        & (Metadata["Medium"] == Medium)
        & (Metadata["CoalescenceType"] == CoalescenceType)
    )
    if Replicate in (1, 2):
        idx = idx & (Metadata["Replicate"] == Replicate)
    elif Replicate != -1:
        raise ValueError("Replicate must be 1, 2, or -1 for both")

    if CommunityOrigin == "S":
        idx = _richness_mask(idx, CoalescenceType, species_pool_num)

    return list(set(Metadata["SampleIDX"][idx]) - set(exception_list))


def getCommunities(Communities_data, IDX_list):
    rows = np.squeeze([np.where(Communities_data["SampleIDX"] == x) for x in IDX_list])
    return Communities_data.iloc[rows]


def getCoalescence(Coalescence_data, IDX_list):
    rows = np.squeeze([np.where(Coalescence_data["SampleIDX"] == x) for x in IDX_list])
    return Coalescence_data.iloc[rows]


# --------------------------------------------------------------------------
# Similarity decomposition
#
# These delegate to coalescence.decomposition so there is one implementation.
# --------------------------------------------------------------------------

def metric_VectorDecomposition_onlyPositive(u, v, m):
    """Decompose ``m`` onto the parental basis; see ``coalescence.decomposition``."""
    return _decompose(u, v, m)


def calculate_assymetricity(u, v, k=None):
    """Retention magnitude and angular asymmetry from decomposition coefficients."""
    x = np.sqrt(np.array(u) ** 2 + np.array(v) ** 2)
    y = np.abs(np.abs(np.arctan2(np.array(u), np.array(v))) - np.pi / 4) / (np.pi / 4)
    return x, y


def characterize_case(x, y):
    """0 Dominance, 1 Mixture, 2 Restructuring."""
    if x ** 2 <= 0.5:
        return 2
    return 0 if y > 0.5 else 1


def uniform_distribution(u, o):
    return (2 * u + 2 * o) * np.random.random() - o


def input_distribution(k=0.8):
    return np.random.exponential(k)


# --------------------------------------------------------------------------
# Precomputed index sets, as the scripts expect them
# --------------------------------------------------------------------------

I = [0.33, 0.66, 0.99]
species_pools = [6, 12, 24]
mediums_pools = ["LN", "MN", "HN"]

Syn_Simulation_IDX = {
    f"I{i}_S{s}": Community_PermutateList_simulation(i, s)
    for i in (1, 2, 3) for s in (1, 2, 3)
}
Syn_Coal_IDX = {
    f"{label}_{size}": Community_PermutateList("F", "S", code, "C", size)
    for label, code in zip(mediums_pools, "LMH") for size in species_pools
}
Syn_Sub_IDX = {
    f"{label}_{size}": Community_PermutateList("F", "S", code, "S", size)
    for label, code in zip(mediums_pools, "LMH") for size in species_pools
}
Nat_Coal_IDX = {
    label: Community_PermutateList("F", "N", code, "C")
    for label, code in zip(mediums_pools, "LMH")
}
Nat_Sub_IDX = {
    label: Community_PermutateList("F", "N", code, "S")
    for label, code in zip(mediums_pools, "LMH")
}
