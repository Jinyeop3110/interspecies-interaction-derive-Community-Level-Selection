"""Stage 4: build the per-community tables.

Python port of ``AnalyzeCommunities_synthetic.m`` and
``AnalyzeCommunities_natural.m``.

Evaluates diversity metrics for every community — parental first, then
coalesced — at each of the five abundance thresholds.  Produces
``processed_Communities_*.xlsx``.
"""

import numpy as np
import pandas as pd

from . import summary_metrics as sm

__all__ = ["analyze_communities"]


def _abundance(sequences, columns, sample_idx):
    match = sequences.index == str(sample_idx)
    if not match.any():
        return np.zeros(len(columns))
    return np.nan_to_num(
        sequences.loc[str(sample_idx), columns].to_numpy(dtype=float).ravel()[: len(columns)]
    )


def analyze_communities(
    metadata,
    sequences,
    media=("L", "M", "H"),
    timepoint="F",
    community_origin="S",
):
    """Run stage 4.

    Communities are emitted parental-first then coalesced, each block ordered by
    medium, matching the row order of the archived tables.

    Returns a DataFrame matching ``processed_Communities_*.xlsx``.
    """
    abundance_columns = [c for c in sequences.columns if c.startswith("NormalizedAbundance")]
    rows = []

    for coalescence_type in ("S", "C"):
        for medium in media:
            selection = metadata[
                (metadata.Timepoint == timepoint)
                & (metadata.CommunityOrigin == community_origin)
                & (metadata.Medium == medium)
                & (metadata.CoalescenceType == coalescence_type)
            ]

            for _, community in selection.iterrows():
                a = _abundance(sequences, abundance_columns, community.SampleIDX)

                row = community.to_dict()
                row["Threshold_num"] = len(sm.THRESHOLD_LEVELS)

                for level, threshold in enumerate(sm.THRESHOLD_LEVELS, start=1):
                    row[f"Threshold_level_{level}"] = threshold
                    row[f"DiversityR_{level}"] = sm.diversity_richness(a, threshold)
                    row[f"DiversitySN_{level}"] = sm.diversity_shannon(a, threshold)
                    row[f"DiversitySS_{level}"] = sm.diversity_simpson(a, threshold)

                rows.append(row)

    return pd.DataFrame(rows)
