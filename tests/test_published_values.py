"""Check that the shared library still reproduces values printed in the paper.

Run with ``pytest`` after downloading the processed data.  These are guard
rails against silent drift in the metric or the classification boundaries, not
a full reproduction of every figure.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from coalescence import io, outcomes  # noqa: E402

# Fig. 4d, main text.  Tolerance is 3 percentage points: a small number of
# events are excluded by the per-figure QC filter that is applied inside the
# figure scripts rather than in the shared library, so fractions computed from
# the unfiltered tables land close to, but not exactly on, the published
# values.  See figures/fig4/README.md.
PUBLISHED_DOMINANCE = {"Nutr-": 0.39, "Base": 0.65, "Nutr+": 0.76}
PUBLISHED_MIXTURE = {"Nutr-": 0.53, "Base": 0.04, "Nutr+": 0.06}
TOLERANCE = 0.03


@pytest.fixture(scope="module")
def fractions():
    try:
        table = outcomes.outcome_table("synthetic")
    except FileNotFoundError as exc:
        pytest.skip(str(exc))
    return outcomes.outcome_fractions(table)


@pytest.mark.parametrize("medium", list(PUBLISHED_DOMINANCE))
def test_dominance_fraction(fractions, medium):
    assert fractions.loc[medium, "Dominance"] == pytest.approx(
        PUBLISHED_DOMINANCE[medium], abs=TOLERANCE
    )


@pytest.mark.parametrize("medium", list(PUBLISHED_MIXTURE))
def test_mixture_fraction(fractions, medium):
    assert fractions.loc[medium, "Mixture"] == pytest.approx(
        PUBLISHED_MIXTURE[medium], abs=TOLERANCE
    )


def test_dominance_increases_along_nutrient_gradient(fractions):
    ordered = [fractions.loc[io.MEDIUM_LABELS[m], "Dominance"] for m in io.MEDIA_ORDER]
    assert ordered == sorted(ordered)


def test_decomposition_is_normalized():
    """x1, x2 and the residual should lie on the unit sphere by construction."""
    from coalescence.decomposition import decompose

    rng = np.random.default_rng(0)
    for _ in range(50):
        u, v, m = (rng.random(20) for _ in range(3))
        x1, x2, residual = decompose(u, v, m)
        assert x1**2 + x2**2 + residual**2 == pytest.approx(1.0)
