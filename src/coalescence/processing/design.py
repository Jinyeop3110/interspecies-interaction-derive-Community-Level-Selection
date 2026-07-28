"""Stage 0: the experimental design tables.

Python port of ``MetadataGenerator.m``, ``RecipeGenerator.m`` and
``CoalescenceRecipe.m``.

These scripts contain no measurement — they transcribe the plate layouts and the
parental pairings that define the experiment.  Every array below is a literal
transcription of the corresponding MATLAB array, kept in the same grouping and
order so the two can be compared line by line.

Produces ``Metadata.xlsx`` (structural columns; the OD, pH and growth-curve
columns are appended by :mod:`coalescence.processing.measurements`) and
``CoalescenceRecipe.xlsx``.

Sample naming
-------------
Each community occupies one well of one plate.  ``SampleIDX`` is
``P<plate>-<well>`` with the well zero-padded to two digits, and ``IDX`` is the
same coordinates as ``100 * plate + well``.

Group keys
----------
``F_S_L_S_1`` reads as timepoint ``F`` (final), origin ``S`` (synthetic),
medium ``L`` (Nutr-), coalescence type ``S`` (parental) or ``C`` (coalesced),
replicate ``1``.
"""

import numpy as np
import pandas as pd

__all__ = [
    "PLATE_LAYOUT",
    "COALESCENCE_PAIRS",
    "NATURAL_COALESCENCE_PAIRS",
    "build_metadata",
    "build_coalescence_recipe",
]


def _seq(start, count):
    """MATLAB's ``start + (1:count)``."""
    return [start + i for i in range(1, count + 1)]


# --------------------------------------------------------------------------
# Plate layout: group key -> (plate numbers, well numbers)
# Transcribed from MetadataGenerator.m lines 23-83.
# --------------------------------------------------------------------------

PLATE_LAYOUT = {
    # Synthetic parental communities
    "F_S_L_S_1": ([1] * 9 + [1] * 9 + [2] * 9 + [1] * 3,
                  _seq(0, 9) + _seq(24, 9) + _seq(48, 9) + _seq(48, 3)),
    "F_S_L_S_2": ([1] * 9 + [1] * 9 + [2] * 9 + [1] * 3,
                  _seq(12, 9) + _seq(36, 9) + _seq(60, 9) + _seq(51, 3)),
    "F_S_M_S_1": ([8] * 9 + [8] * 9 + [2] * 9 + [1] * 3,
                  _seq(48, 9) + _seq(72, 9) + _seq(72, 9) + _seq(54, 3)),
    "F_S_M_S_2": ([8] * 9 + [8] * 9 + [2] * 9 + [1] * 3,
                  _seq(60, 9) + _seq(84, 9) + _seq(84, 9) + _seq(57, 3)),
    "F_S_H_S_1": ([2] * 9 + [2] * 9 + [3] * 12,
                  _seq(0, 9) + _seq(24, 9) + _seq(72, 12)),
    "F_S_H_S_2": ([2] * 9 + [2] * 9 + [3] * 12,
                  _seq(12, 9) + _seq(36, 9) + _seq(84, 12)),

    # Synthetic coalesced communities.  Well 42 is absent from every plate.
    "F_S_L_C_1": ([4] * 47, list(range(1, 42)) + list(range(43, 49))),
    "F_S_L_C_2": ([4] * 47, [97 - w for w in list(range(1, 42)) + list(range(43, 49))]),
    "F_S_M_C_1": ([5] * 47, list(range(1, 42)) + list(range(43, 49))),
    "F_S_M_C_2": ([5] * 47, [97 - w for w in list(range(1, 42)) + list(range(43, 49))]),
    "F_S_H_C_1": ([6] * 47, list(range(1, 42)) + list(range(43, 49))),
    "F_S_H_C_2": ([6] * 47, [97 - w for w in list(range(1, 42)) + list(range(43, 49))]),

    # Natural parental communities, six environmental samples spaced 12 wells apart
    "F_N_L_S_1": ([3] * 6, [1, 13, 25, 37, 49, 61]),
    "F_N_L_S_2": ([3] * 6, [1 + w for w in [1, 13, 25, 37, 49, 61]]),
    "F_N_M_S_1": ([3] * 6, [2 + w for w in [1, 13, 25, 37, 49, 61]]),
    "F_N_M_S_2": ([3] * 6, [3 + w for w in [1, 13, 25, 37, 49, 61]]),
    "F_N_H_S_1": ([3] * 6, [4 + w for w in [1, 13, 25, 37, 49, 61]]),
    "F_N_H_S_2": ([3] * 6, [5 + w for w in [1, 13, 25, 37, 49, 61]]),
}

# Natural coalesced communities share one base well pattern, offset per medium
# and replicate.
_NATURAL_COALESCED_WELLS = [1, 13, 25, 37, 49, 61, 73, 85] + \
                           [1 + w for w in [1, 13, 25, 37, 49, 61, 73]]

for _key, _offset in [
    ("F_N_L_C_1", 0), ("F_N_L_C_2", 2),
    ("F_N_M_C_1", 4), ("F_N_M_C_2", 6),
    ("F_N_H_C_1", 8), ("F_N_H_C_2", 10),
]:
    PLATE_LAYOUT[_key] = ([7] * 15, [_offset + w for w in _NATURAL_COALESCED_WELLS])


# --------------------------------------------------------------------------
# Coalescence pairings.  Transcribed from RecipeGenerator.m.
# --------------------------------------------------------------------------

#: Parental community pairs by richness block, as ``(offset, pairs)``.
#: Community indices 1-9 are 6-species pools, 10-18 are 12-species and
#: 19-30 are 24-species; the offsets place each block's pairs correctly.
_PAIRS_S6 = [(1, 2), (1, 3), (1, 7), (2, 5), (2, 8), (3, 4), (3, 8),
             (4, 5), (4, 9), (5, 6), (5, 9), (6, 7), (6, 9), (7, 8)]

_PAIRS_S12 = [(1, 3), (1, 4), (1, 5), (1, 6), (1, 8), (1, 9),
              (2, 3), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9),
              (3, 5), (3, 7), (3, 8), (3, 9),
              (4, 5), (4, 6), (4, 7), (4, 8), (4, 9),
              (5, 7), (5, 9), (6, 7), (6, 8), (6, 9), (7, 8)]

_PAIRS_S24 = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12)]

#: ``(coalesced index offset, pair offset, pairs)`` per richness block.
COALESCENCE_PAIRS = [
    (0, 0, _PAIRS_S6),
    (14, 9, _PAIRS_S12),
    (41, 18, _PAIRS_S24),
]

#: Natural communities: six samples, fifteen pairings.
NATURAL_COALESCENCE_PAIRS = [
    (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1), (1, 4), (3, 6),
    (1, 3), (2, 4), (3, 5), (4, 6), (5, 1), (6, 2), (2, 5),
]


def build_metadata():
    """Build the structural columns of ``Metadata.xlsx``.

    Returns one row per community with ``SampleIDX``, ``IDX``, the five parsed
    design factors, and ``CommunityIDX``.  The OD, pH and growth-curve columns
    are appended separately; see
    :func:`coalescence.processing.measurements.attach_measurements`.
    """
    rows = []
    for key, (plates, wells) in PLATE_LAYOUT.items():
        timepoint, origin, medium, coalescence_type, replicate = key.split("_")
        for community_index, (plate, well) in enumerate(zip(plates, wells), start=1):
            rows.append(
                {
                    "SampleIDX": f"P{plate}-{well:02d}",
                    "IDX": 100 * plate + well,
                    "Timepoint": timepoint,
                    "CommunityOrigin": origin,
                    "Medium": medium,
                    "CoalescenceType": coalescence_type,
                    "Replicate": int(replicate),
                    "CommunityIDX": community_index,
                }
            )
    return pd.DataFrame(rows)


def build_coalescence_recipe():
    """Build both sheets of ``CoalescenceRecipe.xlsx``.

    Returns ``(synthetic, natural)``, each mapping a coalesced community index
    to the two parental community indices that were mixed to produce it.
    """
    coalesced, sub1, sub2 = [], [], []
    for coal_offset, pair_offset, pairs in COALESCENCE_PAIRS:
        for i, (a, b) in enumerate(pairs, start=1):
            coalesced.append(coal_offset + i)
            sub1.append(pair_offset + a)
            sub2.append(pair_offset + b)

    synthetic = pd.DataFrame(
        {
            "CommunityIDX_Coal": coalesced,
            "CommunityIDX_Sub1": sub1,
            "CommunityIDX_Sub2": sub2,
        }
    )

    natural = pd.DataFrame(
        {
            "CommunityIDX_Coal_natural": list(range(1, len(NATURAL_COALESCENCE_PAIRS) + 1)),
            "CommunityIDX_Sub1_natural": [a for a, _ in NATURAL_COALESCENCE_PAIRS],
            "CommunityIDX_Sub2_natural": [b for _, b in NATURAL_COALESCENCE_PAIRS],
        }
    )

    return synthetic, natural
