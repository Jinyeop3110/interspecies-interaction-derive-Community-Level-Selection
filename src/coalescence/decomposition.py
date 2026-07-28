"""Similarity-space decomposition and coalescence outcome classification.

This is the canonical implementation of the metric used throughout the paper.
It was consolidated from the per-figure notebooks, where it was previously
duplicated (``metric_VectorDecomposition_onlyPositive`` /
``calculate_assymetricity`` / ``characterize_case``).

Given a coalesced community ``m`` and its two parental communities ``u`` and
``v``, all as relative-abundance vectors over a shared ASV index, we express
``m`` in the (generally non-orthogonal) basis spanned by ``u`` and ``v``.  The
two resulting coefficients are the coordinates of the coalescence event in the
two-parent similarity map; the residual is the part of the coalesced community
that neither parent explains.

Note that the parental basis is oblique: ``u`` and ``v`` are not orthogonal
whenever the parental communities share species, so the coefficients come from
solving the 2x2 Gram system rather than from taking two independent cosine
similarities.  For disjoint parental communities the two agree.
"""

import numpy as np

__all__ = [
    "normalize",
    "decompose",
    "retention_and_asymmetry",
    "classify",
    "OUTCOME_LABELS",
    "RETENTION_BOUNDARY",
    "ASYMMETRY_BOUNDARY",
]

#: Outcome codes, in the order used by the stacked-fraction panels.
OUTCOME_LABELS = {0: "Dominance", 1: "Mixture", 2: "Restructuring"}

#: Restructuring boundary, applied to r**2 (equivalently r <= 1/sqrt(2)).
RETENTION_BOUNDARY = 0.5

#: Dominance/Mixture boundary on the normalized angular asymmetry.
ASYMMETRY_BOUNDARY = 0.5


def normalize(v):
    """Return ``v`` scaled to unit L2 norm; zero vectors are returned as-is."""
    v = np.asarray(v, dtype=float)
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def decompose(u, v, m):
    """Decompose community ``m`` onto the parental basis ``(u, v)``.

    Parameters
    ----------
    u, v : array_like
        Relative-abundance vectors of the two parental communities.
    m : array_like
        Relative-abundance vector of the coalesced community.

    Returns
    -------
    x1, x2 : float
        Non-negative coordinates in the two-parent similarity map, rescaled so
        that ``x1**2 + x2**2 + residual**2 == 1``.
    residual : float
        Norm of the component of ``m`` outside the parental span.

    Negative coefficients are clipped to zero before rescaling: a parental
    community cannot contribute a negative amount of composition, and in
    practice negative values arise only for events already deep in the
    Restructuring region.
    """
    u, v, m = normalize(u), normalize(v), normalize(m)

    gram = np.array([[u @ u, u @ v], [u @ v, v @ v]])
    coef = np.linalg.inv(gram) @ np.array([m @ u, m @ v])

    x1 = coef[0] * (coef[0] > 0)
    x2 = coef[1] * (coef[1] > 0)
    residual = np.linalg.norm(m - coef[0] * u - coef[1] * v)

    if x1 == 0 and x2 == 0:
        # Neither parent contributes positively: the event carries no parental
        # signal at all, so it sits at the origin of the map (Restructuring).
        return 0.0, 0.0, residual

    scale = np.sqrt((1 - residual ** 2) / (x1 ** 2 + x2 ** 2))
    return scale * x1, scale * x2, residual


def retention_and_asymmetry(x1, x2):
    """Convert map coordinates to retention magnitude and angular asymmetry.

    Returns
    -------
    r : float
        Retention magnitude, ``sqrt(x1**2 + x2**2)``.  How much of the
        coalesced composition the two parents jointly explain.
    a : float
        Angular asymmetry in ``[0, 1]``.  0 means the event sits on the
        diagonal (both parents contribute equally); 1 means it sits on an axis
        (one parent explains everything).

    The parental dominance index reported in the paper is a monotone
    reparameterization of the same angle:
    ``PDI = x1 / (x1 + x2)``, so ``a`` near 1 corresponds to PDI near 0 or 1.
    """
    x1, x2 = np.asarray(x1, dtype=float), np.asarray(x2, dtype=float)
    r = np.sqrt(x1 ** 2 + x2 ** 2)
    # arctan2 gives the same angle as arctan(x1/x2) on the positive quadrant
    # but stays finite when x2 == 0 (complete dominance by parent 1).
    angle = np.arctan2(x1, x2)
    a = np.abs(np.abs(angle) - np.pi / 4) / (np.pi / 4)
    return r, a


def classify(r, a):
    """Classify one coalescence outcome from retention ``r`` and asymmetry ``a``.

    Returns 0 (Dominance), 1 (Mixture) or 2 (Restructuring); see
    ``OUTCOME_LABELS``.
    """
    if r ** 2 <= RETENTION_BOUNDARY:
        return 2
    return 0 if a > ASYMMETRY_BOUNDARY else 1
