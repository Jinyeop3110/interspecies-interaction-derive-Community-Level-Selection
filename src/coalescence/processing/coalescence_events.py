"""Stage 3: build the per-coalescence-event tables.

Python port of ``AnalyzeCoalescence_synthetic.m`` and
``AnalyzeCoalescence_natural.m``.

For every coalesced community, look up its two parental communities through the
coalescence recipe, fetch all three composition vectors, and evaluate the
summary metrics at each of the five abundance thresholds.  Produces
``processed_CoalescenceEvent_*.xlsx``.
"""

import numpy as np
import pandas as pd

from . import summary_metrics as sm

__all__ = ["analyze_coalescence_events"]


def _abundance(sequences, columns, sample_idx):
    """Composition vector for a sample; zeros when the sample is missing.

    The MATLAB ``GetAbundance`` returns a zero vector rather than failing when a
    sample has no sequencing record, so downstream metrics see an empty
    community.  Preserved here.
    """
    match = sequences.index == str(sample_idx)
    if not match.any():
        return np.zeros(len(columns))
    return np.nan_to_num(sequences.loc[str(sample_idx), columns].to_numpy(dtype=float).ravel()[: len(columns)])


def _parent_community_indices(recipe, coalesced_index, columns):
    """Map a coalesced community index to its two parental community indices."""
    coal_col, sub1_col, sub2_col = columns
    row = recipe[recipe[coal_col] == coalesced_index]
    if row.empty:
        return None, None
    return int(row[sub1_col].iloc[0]), int(row[sub2_col].iloc[0])


def analyze_coalescence_events(
    metadata,
    sequences,
    recipe,
    media=("L", "M", "H"),
    timepoint="F",
    community_origin="S",
    recipe_columns=("CommunityIDX_Coal", "CommunityIDX_Sub1", "CommunityIDX_Sub2"),
):
    """Run stage 3.

    Parameters
    ----------
    metadata : DataFrame
        ``Metadata.xlsx``: one row per sample, with the experimental factors and
        the OD/pH readings.
    sequences : DataFrame
        ``processed_Sequences_*.xlsx`` from stage 1, indexed by ``SampleIDX``.
    recipe : DataFrame
        ``CoalescenceRecipe.xlsx``: which parental communities were mixed.
    media : tuple
        Medium codes to process, in output order.
    timepoint, community_origin : str
        Metadata filters; ``F`` is the final timepoint, ``S`` synthetic
        communities and ``N`` natural ones.
    recipe_columns : tuple
        Column names in ``recipe``; natural communities use the ``_natural``
        suffixed variants.

    Returns
    -------
    DataFrame
        One row per coalescence event, matching the schema of
        ``processed_CoalescenceEvent_*.xlsx``.
    """
    abundance_columns = [c for c in sequences.columns if c.startswith("NormalizedAbundance")]
    rows = []

    # Samples with no sequencing record give zero vectors, so several metrics
    # divide by zero and yield NaN. That is what the MATLAB does, and the
    # archived tables carry those NaNs, so the warnings are silenced rather than
    # the values changed.
    with np.errstate(divide="ignore", invalid="ignore"):
        rows = _collect_events(
            metadata, sequences, recipe, media, timepoint, community_origin,
            recipe_columns, abundance_columns,
        )

    return pd.DataFrame(rows)


def _collect_events(metadata, sequences, recipe, media, timepoint, community_origin,
                    recipe_columns, abundance_columns):
    """Event loop for :func:`analyze_coalescence_events`."""
    rows = []

    for medium in media:
        selection = metadata[
            (metadata.Timepoint == timepoint)
            & (metadata.CommunityOrigin == community_origin)
            & (metadata.Medium == medium)
            & (metadata.CoalescenceType == "C")
        ]

        for _, coalesced in selection.iterrows():
            sub1_idx, sub2_idx = _parent_community_indices(
                recipe, int(coalesced.CommunityIDX), recipe_columns
            )
            if sub1_idx is None:
                continue

            # The parents share every experimental factor with their offspring
            # apart from the coalescence type and the community index.
            def parent(community_index):
                match = metadata[
                    (metadata.Timepoint == coalesced.Timepoint)
                    & (metadata.CommunityOrigin == coalesced.CommunityOrigin)
                    & (metadata.Medium == coalesced.Medium)
                    & (metadata.CoalescenceType == "S")
                    & (metadata.Replicate == coalesced.Replicate)
                    & (metadata.CommunityIDX == community_index)
                ]
                return match.iloc[0] if len(match) else None

            parent1, parent2 = parent(sub1_idx), parent(sub2_idx)
            if parent1 is None or parent2 is None:
                continue

            a = _abundance(sequences, abundance_columns, coalesced.SampleIDX)
            a1 = _abundance(sequences, abundance_columns, parent1.SampleIDX)
            a2 = _abundance(sequences, abundance_columns, parent2.SampleIDX)

            row = coalesced.to_dict()
            row["SampleIDX_Sub1"] = parent1.SampleIDX
            row["SampleIDX_Sub2"] = parent2.SampleIDX
            row["Threshold_num"] = len(sm.THRESHOLD_LEVELS)

            for level, threshold in enumerate(sm.THRESHOLD_LEVELS, start=1):
                row[f"Threshold_level_{level}"] = threshold

                # Column order matters: the archived tables list every
                # SimilarityTo1 column before any SimilarityTo2 column.
                for name, similarity in sm.SIMILARITY_FUNCTIONS.items():
                    row[f"SimilarityTo1_{name}_{level}"] = sm.similarity_to_1(
                        a, a1, a2, threshold, similarity)
                for name, similarity in sm.SIMILARITY_FUNCTIONS.items():
                    row[f"SimilarityTo2_{name}_{level}"] = sm.similarity_to_2(
                        a, a1, a2, threshold, similarity)

                row[f"Additivity1_{level}"] = sm.additivity_1(a, a1, a2, threshold)
                row[f"Additivity2_{level}"] = sm.additivity_2(a, a1, a2, threshold)
                row[f"Additivity3_{level}"] = sm.additivity_3(a, a1, a2, threshold)

                for name, similarity in sm.SIMILARITY_FUNCTIONS.items():
                    row[f"Assymetricity_{name}_{level}"] = sm.asymmetry_from_similarity(
                        a, a1, a2, threshold, similarity)

                row[f"Assymetricity1_{level}"] = sm.asymmetry_1(a, a1, a2, threshold)
                row[f"Assymetricity2_{level}"] = sm.asymmetry_2(a, a1, a2, threshold)
                row[f"Assymetricity3_{level}"] = sm.asymmetry_3(a, a1, a2, threshold)

                row[f"ASVs_notinsub1_{level}"] = sm.asvs_fraction_not_in_parents_count(
                    a, a1, a2, threshold)
                row[f"ASVs_notinsub2_{level}"] = sm.asvs_fraction_not_in_parents_abundance(
                    a, a1, a2, threshold)
                row[f"ASVs_overlap1_{level}"] = sm.asvs_fraction_overlap_count(
                    a, a1, a2, threshold)
                row[f"ASVs_overlap2_{level}"] = sm.asvs_fraction_overlap_abundance(
                    a, a1, a2, threshold)

            rows.append(row)

    return rows
