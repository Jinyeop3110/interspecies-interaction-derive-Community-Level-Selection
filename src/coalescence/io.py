"""Loading of the processed data tables.

Every figure script starts here.  The tables themselves are not in this
repository; download them from Dryad and point ``COALESCENCE_DATA`` at the
directory (see ``data/README.md``).
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "data_dir",
    "load_coalescence_events",
    "load_communities",
    "load_metadata",
    "load_coalescence_recipe",
    "load_sequences",
    "sample_vector",
    "composition_matrix",
    "MEDIUM_LABELS",
    "MEDIA_ORDER",
    "EXCLUDED_SAMPLES",
]

#: Medium codes as they appear in the processed tables, mapped to the names
#: used in the paper.
MEDIUM_LABELS = {"L": "Nutr-", "M": "Base", "H": "Nutr+"}

#: Canonical left-to-right ordering along the nutrient gradient.
MEDIA_ORDER = ["L", "M", "H"]

#: Samples excluded from the analysis on quality-control grounds, with the
#: reason each was dropped.
#:
#: Applied uniformly across every figure in the study.  Plates 4, 5 and 6 hold
#: the coalesced communities in Nutr-, Base and Nutr+ respectively; plate 7 the
#: natural coalesced communities; plate 8 Base parental communities.
#:
#: Applying this list is what makes the outcome fractions match the published
#: values exactly (n = 90/83/90 for Nutr-/Base/Nutr+).
EXCLUSION_REASONS = {
    "P4-02": "no reads in the coalescence result",
    "P4-03": "no reads in the coalescence result",
    "P4-23": "no reads in the coalescence result",
    "P4-24": "no reads in the coalescence result",
    "P7-97": "no reads in the coalescence result",
    "P8-12": "no reads in the coalescence result",
    "P8-91": "sequencing file missing",
    "P5-56": "parental community MN E7 missing",
    "P5-59": "parental community MN E7 missing",
    "P5-61": "parental community MN E7 missing",
    "P5-64": "parental community MN E7 missing",
    "P5-69": "parental community MN E7 missing",
    "P5-73": "parental community MN E7 missing",
    "P5-47": "unexpected ASVs, likely mislabelling",
    "P5-50": "unexpected ASVs, likely mislabelling",
    "P5-39": "more than 30% of abundance from ASVs in neither parent",
    "P5-54": "more than 30% of abundance from ASVs in neither parent",
    "P5-87": "more than 30% of abundance from ASVs in neither parent",
    "P6-02": "more than 30% of abundance from ASVs in neither parent",
    "P6-47": "more than 30% of abundance from ASVs in neither parent",
    "P6-57": "more than 30% of abundance from ASVs in neither parent",
    "P6-74": "more than 30% of abundance from ASVs in neither parent",
}

EXCLUDED_SAMPLES = frozenset(EXCLUSION_REASONS)

_DEFAULT_DIR = Path(__file__).resolve().parents[2] / "data"


def data_dir():
    """Directory holding the processed tables.

    Override with the ``COALESCENCE_DATA`` environment variable.
    """
    return Path(os.environ.get("COALESCENCE_DATA", _DEFAULT_DIR))


def _read(name):
    path = data_dir() / name
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found.\n"
            "Download the processed data from Dryad and either place it in "
            f"{data_dir()} or set COALESCENCE_DATA to its location. "
            "See data/README.md."
        )
    return pd.read_excel(path)


def load_coalescence_events(source="synthetic"):
    """One row per coalescence event, with per-event summary metrics.

    ``source`` is ``"synthetic"``, ``"natural"`` or ``"simulation"``.
    """
    return _read(f"processed_CoalescenceEvent_{source}.xlsx")


def load_communities(source="synthetic"):
    """One row per community (parental and coalesced), with diversity metrics."""
    return _read(f"processed_Communities_{source}.xlsx")


def load_metadata():
    """One row per community: design factors plus OD, pH and growth-curve readings."""
    return _read("Metadata.xlsx")


def load_coalescence_recipe(source="synthetic"):
    """Which parental communities were mixed to produce each coalesced community."""
    sheet = 0 if source == "synthetic" else 1
    path = data_dir() / "CoalescenceRecipe.xlsx"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; see data/README.md")
    return pd.read_excel(path, sheet_name=sheet)


def load_sequences(source="synthetic"):
    """Per-sample relative abundance across ASVs, indexed by ``SampleIDX``.

    This is the table the similarity-space decomposition operates on; the
    per-event summary tables do not carry the full composition vectors.
    """
    return _read(f"processed_Sequences_{source}.xlsx").set_index("SampleIDX")


def composition_matrix(sequences):
    """Return the abundance column names of a sequence table, in order."""
    return [c for c in sequences.columns if c.startswith("NormalizedAbundance")]


def sample_vector(sequences, columns, sample_idx):
    """Relative-abundance vector for one sample, or ``None`` if unusable.

    Returns ``None`` when the sample is absent from the table or has no
    assigned reads, which is the condition every figure script filters on.
    """
    if sample_idx not in sequences.index:
        return None
    row = sequences.loc[sample_idx, columns].to_numpy(dtype=float)
    row = np.nan_to_num(np.atleast_2d(row)[0])
    if row.sum() <= 0:
        return None
    return row
