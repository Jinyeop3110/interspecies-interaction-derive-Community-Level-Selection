"""Event-level pipeline: processed tables in, classified outcomes out."""

import numpy as np
import pandas as pd

from . import io
from .decomposition import (
    classify, decompose, parental_dominance_index, retention_and_asymmetry,
)

__all__ = ["outcome_table", "outcome_fractions"]


def outcome_table(source="synthetic", apply_exclusions=True):
    """Classify every coalescence event in one experimental dataset.

    Returns a frame with one row per usable event and the columns
    ``SampleIDX``, ``Medium``, ``Replicate``, ``x1``, ``x2``, ``residual``,
    ``r``, ``asymmetry``, ``PDI`` and ``outcome``.

    Events are skipped when the coalesced community or either parent is
    missing from the sequence table or carries no assigned reads.

    ``apply_exclusions`` drops events involving a sample in
    ``io.EXCLUDED_SAMPLES``, the quality-control list used throughout the
    paper.  Leave it on to reproduce published values; turn it off to see the
    unfiltered data.
    """
    events = io.load_coalescence_events(source)
    sequences = io.load_sequences(source)
    columns = io.composition_matrix(sequences)

    rows = []
    for _, event in events.iterrows():
        if apply_exclusions and not {
            str(event.SampleIDX),
            str(event.SampleIDX_Sub1),
            str(event.SampleIDX_Sub2),
        }.isdisjoint(io.EXCLUDED_SAMPLES):
            continue
        m = io.sample_vector(sequences, columns, event.SampleIDX)
        u = io.sample_vector(sequences, columns, event.SampleIDX_Sub1)
        v = io.sample_vector(sequences, columns, event.SampleIDX_Sub2)
        if m is None or u is None or v is None:
            continue

        x1, x2, residual = decompose(u, v, m)
        if not np.isfinite(x1) or not np.isfinite(x2):
            continue

        r, asymmetry = retention_and_asymmetry(x1, x2)
        rows.append(
            {
                "SampleIDX": event.SampleIDX,
                "Medium": event.Medium,
                "Replicate": event.get("Replicate", np.nan),
                "x1": x1,
                "x2": x2,
                "residual": residual,
                "r": r,
                "asymmetry": asymmetry,
                "PDI": float(parental_dominance_index(x1, x2)),
                "outcome": classify(r, asymmetry),
            }
        )

    return pd.DataFrame(rows)


def outcome_fractions(table):
    """Fraction of each outcome per medium, ordered along the nutrient gradient.

    Returns a frame indexed by medium label with columns ``Dominance``,
    ``Mixture``, ``Restructuring`` and ``n``.
    """
    rows = {}
    for medium in io.MEDIA_ORDER:
        group = table[table.Medium == medium]
        if group.empty:
            continue
        rows[io.MEDIUM_LABELS[medium]] = {
            "Dominance": (group.outcome == 0).mean(),
            "Mixture": (group.outcome == 1).mean(),
            "Restructuring": (group.outcome == 2).mean(),
            "n": len(group),
        }
    return pd.DataFrame(rows).T
