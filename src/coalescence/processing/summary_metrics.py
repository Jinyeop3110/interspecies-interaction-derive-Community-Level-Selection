"""Per-event and per-community summary metrics.

Faithful ports of the metric functions defined at the bottom of
``AnalyzeCoalescence_*.m`` and ``AnalyzeCommunities_*.m``.  They populate the
summary columns of the processed tables.

These are *not* the metric behind the published similarity maps — that is
``coalescence.decomposition``, which the figure scripts call directly on the
composition vectors.  The columns computed here are retained for cross-checking
and for the supplementary threshold-sensitivity analyses.

Every function takes relative-abundance vectors and an abundance ``threshold``.
``A`` is the coalesced community, ``A1`` and ``A2`` its two parents.
"""

import numpy as np

__all__ = [
    "THRESHOLD_LEVELS",
    "diversity_richness",
    "diversity_shannon",
    "diversity_simpson",
    "similarity_bray_curtis",
    "similarity_jensen_shannon",
    "similarity_jaccard",
    "similarity_dot_product",
    "similarity_to_1",
    "similarity_to_2",
    "additivity_1",
    "additivity_2",
    "additivity_3",
    "asymmetry_from_similarity",
    "asymmetry_1",
    "asymmetry_2",
    "asymmetry_3",
    "asvs_fraction_not_in_parents_count",
    "asvs_fraction_not_in_parents_abundance",
    "asvs_fraction_overlap_count",
    "asvs_fraction_overlap_abundance",
    "SIMILARITY_FUNCTIONS",
]

#: The five abundance thresholds the summary columns are evaluated at.
THRESHOLD_LEVELS = [0.1, 0.033, 0.01, 0.0033, 0.001]


# --------------------------------------------------------------------------
# Diversity
# --------------------------------------------------------------------------

def diversity_richness(a, threshold):
    """Number of ASVs above threshold."""
    return np.sum(a > threshold)


def diversity_shannon(a, threshold):
    """Shannon entropy. Note the 1e-6 offset inside the log, as in MATLAB."""
    return -np.sum(a * np.log(a + 1e-6))


def diversity_simpson(a, threshold):
    """Inverse Simpson index; zero for samples that lost most of their mass."""
    if np.sum(a) < 0.9:
        return 0.0
    return 1.0 / np.sum(a ** 2)


# --------------------------------------------------------------------------
# Pairwise similarity
# --------------------------------------------------------------------------

def similarity_bray_curtis(a, b, threshold):
    """Shared abundance: sum of the elementwise minimum."""
    a = a * (a > threshold)
    b = b * (b > threshold)
    return np.sum(np.minimum(a, b))


def similarity_jensen_shannon(a, b, threshold):
    """One minus the Jensen-Shannon distance."""
    a = a * (a > threshold)
    b = b * (b > threshold)
    a_, b_ = a + 1e-10, b + 1e-10
    mean = (a_ + b_) / 2
    divergence = 0.5 * np.sum(a_ * (np.log(a_) - np.log(mean))) + \
                 0.5 * np.sum(b_ * (np.log(b_) - np.log(mean)))
    return 1 - np.sqrt(divergence)


def similarity_jaccard(a, b, threshold):
    """Jaccard index over presence/absence."""
    a = a * (a > threshold)
    b = b * (b > threshold)
    union = np.sum((a > threshold) | (b > threshold))
    return np.sum((a > threshold) & (b > threshold)) / union


def similarity_dot_product(a, b, threshold):
    """The ``DP`` columns.

    Despite the name this is not a cosine similarity: each vector is divided by
    its *sum of squares* and then square-rooted before the inner product, giving
    a Bhattacharyya-style coefficient.  Reproduced as written, since these
    columns exist in the archived tables.
    """
    a = a * (a > threshold)
    b = b * (b > threshold)
    return np.sum(np.sqrt(a / np.sum(a ** 2)) * np.sqrt(b / np.sum(b ** 2)))


#: Suffixes used in the processed-table column names.
SIMILARITY_FUNCTIONS = {
    "BC": similarity_bray_curtis,
    "J": similarity_jaccard,
    "JS": similarity_jensen_shannon,
    "DP": similarity_dot_product,
}


def similarity_to_1(a, a1, a2, threshold, similarity):
    """Share of total similarity attributable to parent 1."""
    s1 = similarity(a, a1, threshold)
    s2 = similarity(a, a2, threshold)
    return s1 / (s1 + s2)


def similarity_to_2(a, a1, a2, threshold, similarity):
    """Share of total similarity attributable to parent 2."""
    s1 = similarity(a, a1, threshold)
    s2 = similarity(a, a2, threshold)
    return s2 / (s1 + s2)


# --------------------------------------------------------------------------
# Additivity: how much of the combined parental richness survives
# --------------------------------------------------------------------------

def additivity_1(a, a1, a2, threshold):
    """Coalesced richness over the union richness of the parents."""
    return np.sum(a > threshold) / (
        np.sum(a1 > threshold) + np.sum(a2 > threshold)
        - np.sum((a1 > threshold) * (a2 > threshold))
    )


def additivity_2(a, a1, a2, threshold):
    """Fraction of the parental union that survives coalescence."""
    present = (a1 > threshold) | (a2 > threshold)
    return np.sum((a > threshold) * present) / np.sum(present)


def additivity_3(a, a1, a2, threshold):
    """As ``additivity_2`` but normalized by summed rather than union richness."""
    present = (a1 > threshold) | (a2 > threshold)
    return np.sum((a > threshold) * present) / (
        np.sum(a1 > threshold) + np.sum(a2 > threshold)
    )


# --------------------------------------------------------------------------
# Asymmetry: which parent the coalesced community leans toward
# --------------------------------------------------------------------------

def asymmetry_from_similarity(a, a1, a2, threshold, similarity):
    """Rescaled parental share, 0 (balanced) to 1 (one parent only)."""
    s1 = similarity(a, a1, threshold)
    s2 = similarity(a, a2, threshold)
    return 2 * abs(s1 / (s1 + s2) - 0.5)


def asymmetry_1(a, a1, a2, threshold):
    """Sign-counted asymmetry over ASVs unique to one parent."""
    numerator = np.sum(np.sign(a * (a1 > threshold) - a * (a2 > threshold)))
    denominator = np.sum(np.sign(a * (a1 > threshold) + a * (a2 > threshold)))
    return abs(numerator / denominator)


def asymmetry_2(a, a1, a2, threshold):
    """Abundance-weighted asymmetry, shared ASVs counted once."""
    numerator = np.sum(a * (a1 > threshold) - a * (a2 > threshold))
    denominator = np.sum(
        a * (a1 > threshold) + a * (a2 > threshold)
        - a * (a1 > threshold) * (a2 > threshold)
    )
    return abs(numerator / denominator)


def asymmetry_3(a, a1, a2, threshold):
    """As ``asymmetry_2`` but shared ASVs fully discounted."""
    numerator = np.sum(a * (a1 > threshold) - a * (a2 > threshold))
    denominator = np.sum(
        a * (a1 > threshold) + a * (a2 > threshold)
        - 2 * a * (a1 > threshold) * (a2 > threshold)
    )
    return abs(numerator / denominator)


# --------------------------------------------------------------------------
# Novel and shared ASVs
# --------------------------------------------------------------------------

def asvs_fraction_not_in_parents_count(a, a1, a2, threshold):
    """Fraction of coalesced ASVs found in neither parent, by count."""
    novel = a * (1 - ((a1 > threshold) | (a2 > threshold)))
    return np.sum(novel > threshold) / (np.sum(a > threshold) + 1e-10)


def asvs_fraction_not_in_parents_abundance(a, a1, a2, threshold):
    """Total abundance of ASVs found in neither parent."""
    return np.sum(a * (1 - ((a1 > threshold) | (a2 > threshold))))


def asvs_fraction_overlap_count(a, a1, a2, threshold):
    """Fraction of coalesced ASVs shared by both parents, by count."""
    shared = a * (a1 > threshold) * (a2 > threshold)
    return np.sum(shared > threshold) / (np.sum(a > threshold) + 1e-10)


def asvs_fraction_overlap_abundance(a, a1, a2, threshold):
    """Total abundance of ASVs shared by both parents."""
    return np.sum(a * (a1 > threshold) * (a2 > threshold))
