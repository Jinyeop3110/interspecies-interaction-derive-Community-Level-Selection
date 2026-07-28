"""Check the Python pipeline against the archived MATLAB output.

The processing stage was ported from MATLAB; these tests are what make that
port trustworthy. Each one regenerates a table from the same inputs the MATLAB
consumed and compares it cell by cell against the archived result.

Requires the raw ASV tables and the design tables in addition to the processed
data, so the tests skip unless ``COALESCENCE_RAW`` points at a tree containing:

    SEQanalysis/excludeNatural/   SEQanalysis/onlyNatural/
    Postprocessed/               Analyzed/
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from coalescence.processing.coalescence_events import analyze_coalescence_events  # noqa: E402
from coalescence.processing.communities import analyze_communities  # noqa: E402
from coalescence.processing.sequences import NATURAL, SYNTHETIC, process_sequences  # noqa: E402

#: Floating-point tolerance. The ports are exact up to accumulation order.
TOLERANCE = 1e-9


def raw_root():
    root = os.environ.get("COALESCENCE_RAW")
    if not root:
        pytest.skip("set COALESCENCE_RAW to the working-tree root to run parity tests")
    return Path(root)


def assert_frames_match(ported, archived, label):
    assert ported.shape == archived.shape, (
        f"{label}: shape {ported.shape} != archived {archived.shape}"
    )
    assert list(ported.columns) == list(archived.columns), f"{label}: column mismatch"

    numeric = [c for c in archived.columns if pd.api.types.is_numeric_dtype(archived[c])]
    difference = np.abs(
        ported[numeric].to_numpy(dtype=float) - archived[numeric].to_numpy(dtype=float)
    )
    worst = np.nanmax(difference)
    assert worst < TOLERANCE, f"{label}: largest difference {worst:.3e}"


@pytest.fixture(scope="module")
def metadata():
    return pd.read_excel(raw_root() / "Postprocessed" / "Metadata.xlsx")


def load_sequences(source):
    frame = pd.read_excel(
        raw_root() / "Postprocessed" / f"processed_Sequences_{source}.xlsx"
    ).set_index("SampleIDX")
    frame.index = frame.index.astype(str)
    return frame


@pytest.mark.parametrize(
    "source,asv_dir,config",
    [
        ("synthetic", "excludeNatural", SYNTHETIC),
        ("natural", "onlyNatural", NATURAL),
    ],
)
def test_sequences_stage(source, asv_dir, config):
    root = raw_root()
    directory = root / "SEQanalysis" / asv_dir
    ported, _ = process_sequences(
        pd.read_csv(directory / "M_OTUtableGreenGenes.csv"),
        pd.read_csv(directory / "M_TAXAtableGreenGenes.csv"),
        pd.read_csv(directory / "M_UniqueSequences.csv"),
        config,
    )
    archived = pd.read_excel(root / "Postprocessed" / f"processed_Sequences_{source}.xlsx")
    assert_frames_match(ported, archived, f"sequences/{source}")


@pytest.mark.parametrize(
    "source,origin,sheet,columns",
    [
        ("synthetic", "S", 0,
         ("CommunityIDX_Coal", "CommunityIDX_Sub1", "CommunityIDX_Sub2")),
        ("natural", "N", 1,
         ("CommunityIDX_Coal_natural", "CommunityIDX_Sub1_natural",
          "CommunityIDX_Sub2_natural")),
    ],
)
def test_coalescence_events_stage(metadata, source, origin, sheet, columns):
    root = raw_root()
    recipe = pd.read_excel(root / "Postprocessed" / "CoalescenceRecipe.xlsx", sheet_name=sheet)
    ported = analyze_coalescence_events(
        metadata, load_sequences(source), recipe,
        community_origin=origin, recipe_columns=columns,
    )
    archived = pd.read_excel(root / "Analyzed" / f"processed_CoalescenceEvent_{source}.xlsx")
    assert_frames_match(ported, archived, f"coalescence_events/{source}")


@pytest.mark.parametrize("source,origin", [("synthetic", "S"), ("natural", "N")])
def test_communities_stage(metadata, source, origin):
    ported = analyze_communities(metadata, load_sequences(source), community_origin=origin)
    archived = pd.read_excel(
        raw_root() / "Analyzed" / f"processed_Communities_{source}.xlsx"
    )
    assert_frames_match(ported, archived, f"communities/{source}")
