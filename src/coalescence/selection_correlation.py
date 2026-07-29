"""Pairwise selection correlation: do species from the same parent share fates?

Implements the fate-concordance metric defined in the Supplementary Methods and
applied in Supplementary Note 7.  It is the measure behind Fig. 2d and
Extended Data Figs. 5 and 6.

The question is whether coalescence selects on communities or on species.  If
species from the same parental community tend to persist or go extinct
*together*, their fates are correlated and selection is acting above the level
of the individual species.  If fates are independent of parental origin, it is
not.

Method
------
Within one coalescence event, every species is given a parental affiliation:

* detected in one parent only        -> that parent
* detected in both parents           -> the parent where it was more abundant
                                        before coalescence (exact ties go to
                                        parent A by convention)
* detected in neither                -> unassigned, excluded from pairing

Each pair of assigned species either shares a fate (both present in the
offspring, or both absent) or does not.  The shared-fate rate is rescaled to a
correlation-style value::

    rho = 2 * shared_fate_rate - 1

so +1 means every pair shares its fate and -1 means none do.  ``rho`` is
computed separately for same-affiliation and cross-affiliation pairs, then
averaged over events.  Significance comes from shuffling affiliation labels
within each event.
"""

import numpy as np

__all__ = [
    "assign_parental_affiliation",
    "concordance",
    "event_concordance",
    "selection_correlation",
    "permutation_test",
]

#: Relative abundance above which a species counts as present.
PRESENCE_THRESHOLD = 1e-3


def assign_parental_affiliation(parent_a, parent_b, threshold=PRESENCE_THRESHOLD):
    """Label each species by which parental community it came from.

    Returns an integer array: ``0`` for parent A, ``1`` for parent B, ``-1``
    for species present in neither parent, which take no part in pairing.
    """
    parent_a = np.asarray(parent_a, dtype=float)
    parent_b = np.asarray(parent_b, dtype=float)

    in_a = parent_a > threshold
    in_b = parent_b > threshold

    affiliation = np.full(len(parent_a), -1, dtype=int)
    affiliation[in_a & ~in_b] = 0
    affiliation[in_b & ~in_a] = 1

    # Shared species go to whichever parent carried more of them. Exact ties go
    # to parent A, as stated in the Supplementary Methods.
    shared = in_a & in_b
    affiliation[shared] = np.where(parent_a[shared] >= parent_b[shared], 0, 1)

    return affiliation


def concordance(survived, pair_indices):
    """Rescale a shared-fate rate over the given pairs to ``[-1, 1]``.

    Returns ``nan`` when there are no pairs to average over.
    """
    if not len(pair_indices):
        return np.nan
    first, second = pair_indices[:, 0], pair_indices[:, 1]
    shared_fate = survived[first] == survived[second]
    return 2 * shared_fate.mean() - 1


def _pairs_within_and_across(affiliation):
    """Index pairs of assigned species, split by same vs different parent."""
    assigned = np.flatnonzero(affiliation >= 0)
    same, cross = [], []
    for i in range(len(assigned)):
        for j in range(i + 1, len(assigned)):
            a, b = assigned[i], assigned[j]
            (same if affiliation[a] == affiliation[b] else cross).append((a, b))
    return np.array(same, dtype=int).reshape(-1, 2), np.array(cross, dtype=int).reshape(-1, 2)


def event_concordance(parent_a, parent_b, coalesced, threshold=PRESENCE_THRESHOLD):
    """Concordance for one coalescence event.

    Returns ``(rho_same, rho_cross)``, either of which may be ``nan`` if that
    class of pair does not occur in this event.
    """
    affiliation = assign_parental_affiliation(parent_a, parent_b, threshold)
    survived = np.asarray(coalesced, dtype=float) > threshold

    same_pairs, cross_pairs = _pairs_within_and_across(affiliation)
    return concordance(survived, same_pairs), concordance(survived, cross_pairs)


def selection_correlation(events, threshold=PRESENCE_THRESHOLD):
    """Average concordance across events.

    ``events`` is an iterable of ``(parent_a, parent_b, coalesced)`` triples of
    relative-abundance vectors.

    Returns ``(rho_same, rho_cross, delta)``, averaged over events that yielded
    a value.  ``delta`` is the quantity the permutation test evaluates.
    """
    same, cross = [], []
    for parent_a, parent_b, coalesced in events:
        rho_same, rho_cross = event_concordance(parent_a, parent_b, coalesced, threshold)
        if np.isfinite(rho_same):
            same.append(rho_same)
        if np.isfinite(rho_cross):
            cross.append(rho_cross)

    mean_same = float(np.mean(same)) if same else np.nan
    mean_cross = float(np.mean(cross)) if cross else np.nan
    return mean_same, mean_cross, mean_same - mean_cross


def permutation_test(events, n_permutations=1000, threshold=PRESENCE_THRESHOLD, seed=0):
    """Permutation test on the same-minus-cross concordance difference.

    The null shuffles parental-affiliation labels *within* each event, which
    preserves each event's composition and survival pattern and destroys only
    the association between parental origin and fate.

    Returns ``(observed_delta, p_value, null_distribution)``.  The p-value is
    one-sided: the fraction of permutations whose delta is at least the
    observed one, with the usual ``+1`` correction so it is never zero.
    """
    rng = np.random.default_rng(seed)
    observed = selection_correlation(events, threshold)[2]

    prepared = []
    for parent_a, parent_b, coalesced in events:
        affiliation = assign_parental_affiliation(parent_a, parent_b, threshold)
        survived = np.asarray(coalesced, dtype=float) > threshold
        prepared.append((affiliation, survived))

    null = np.empty(n_permutations)
    for permutation in range(n_permutations):
        same, cross = [], []
        for affiliation, survived in prepared:
            shuffled = affiliation.copy()
            assigned = np.flatnonzero(shuffled >= 0)
            shuffled[assigned] = rng.permutation(shuffled[assigned])

            same_pairs, cross_pairs = _pairs_within_and_across(shuffled)
            rho_same = concordance(survived, same_pairs)
            rho_cross = concordance(survived, cross_pairs)
            if np.isfinite(rho_same):
                same.append(rho_same)
            if np.isfinite(rho_cross):
                cross.append(rho_cross)
        null[permutation] = (np.mean(same) if same else np.nan) - \
                            (np.mean(cross) if cross else np.nan)

    p_value = (np.sum(null >= observed) + 1) / (n_permutations + 1)
    return observed, float(p_value), null
