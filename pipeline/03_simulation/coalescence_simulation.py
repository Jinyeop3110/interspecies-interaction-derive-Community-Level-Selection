#!/usr/bin/env python3
"""Generalized Lotka-Volterra simulation of community coalescence.

This is the simulation behind Figures 2 and 3.  It replaces the demonstration
script that previously stood in for it, which did not implement the sampling
design used in the paper.

Model
-----
Species grow logistically and compete pairwise::

    dn_i/dt = g_i n_i (1 - sum_j I_ij n_j / k_i)

with growth rates and carrying capacities fixed to 1 and self-interaction
``I_ii = 1``.  Off-diagonal competition coefficients are drawn from
``U(mu - spread, mu + spread + ...)``; see :func:`uniform_interaction`.  The
parameter ``mu`` therefore sets both the mean and the spread of competitive
effects.

Design
------
A pool of ``n_communities * species_per_community`` species is partitioned into
non-overlapping parental communities of equal size.  Each parental community is
assembled to equilibrium on its own, then pairs are mixed at their equilibrium
abundances and re-equilibrated.  This mirrors the experimental protocol of
separate parental assembly followed by pairwise coalescence.

Reproducibility
---------------
Every replicate seeds NumPy's global generator with the replicate index before
drawing anything, so a sweep is fully determined by its parameters.  The order
of random draws is load-bearing: the interaction matrix consumes ``N * N``
values and the initial abundances ``N`` more, in that order.  Changing the
order changes the results even at the same seed.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from coalescence.decomposition import (  # noqa: E402
    classify, decompose, retention_and_asymmetry,
)
from coalescence.lv import run_lotka_volterra  # noqa: E402

#: Abundance below which a species is treated as extinct.
EXTINCTION_THRESHOLD = 1e-3

#: Integration window for each equilibration.
TIME_SPAN = [0, 5000]


def uniform_interaction(mu, spread=0.0):
    """Draw one competition coefficient.

    Reproduces the ``uniform_distribution(u, o)`` helper from the notebooks:
    ``(2*mu + 2*spread) * random() - spread``.  With ``spread = 0`` this is
    ``U(0, 2*mu)``, the distribution described in the Methods.
    """
    return (2 * mu + 2 * spread) * np.random.random() - spread


def initialize_species_pool(n_species, mu, spread=0.0):
    """Build the interaction matrix, growth rates and carrying capacities.

    Consumes exactly ``n_species ** 2`` random draws, matching the loop order
    of ``InitializeSpeceiesPool`` in the notebooks.  Growth rates and carrying
    capacities are constant and consume no randomness.
    """
    interactions = np.zeros((n_species, n_species))
    for i in range(n_species):
        for j in range(n_species):
            interactions[i, j] = uniform_interaction(mu, spread)
    np.fill_diagonal(interactions, 1.0)

    growth_rates = np.ones(n_species)
    carrying_capacities = np.ones(n_species)
    return interactions, growth_rates, carrying_capacities


def initialize_community_pool(n_species, n_communities, species_per_community):
    """Partition the pool into non-overlapping parental communities.

    Community ``i`` takes the consecutive block of species
    ``[i * species_per_community, (i + 1) * species_per_community)``, so no two
    parental communities share a species.
    """
    library = np.zeros((n_communities, n_species))
    for i in range(n_communities):
        block = np.arange(species_per_community * i, species_per_community * (i + 1))
        library[i, block] = 1
    return library


def simulate_one_replicate(mu, seed, n_communities=2, species_per_community=12,
                           spread=0.0, pool_size=None, return_interactions=False):
    """Assemble parental communities and coalesce every pair.

    ``pool_size`` is the number of species in the pool.  It defaults to exactly
    ``n_communities * species_per_community``, but may be larger, in which case
    the surplus species are present in the pool and in the interaction matrix
    but belong to no parental community.  The manuscript draws four 12-species
    communities from a pool of 54, so the pool is six species larger than the
    communities consume.  Pool size changes both the dynamics and the number of
    random draws, so it must match to reproduce a given sweep.

    Returns a list of ``(parent_a, parent_b, coalesced)`` equilibrium abundance
    vectors, one per unordered pair of parental communities.
    """
    n_species = pool_size or n_communities * species_per_community
    if n_species < n_communities * species_per_community:
        raise ValueError("pool_size is smaller than the communities require")
    np.random.seed(seed)

    interactions, growth_rates, carrying_capacities = initialize_species_pool(
        n_species, mu, spread
    )
    library = initialize_community_pool(n_species, n_communities, species_per_community)
    initial = np.random.rand(n_species) * 0.1

    # Assemble each parental community separately.
    parents = []
    for community in library:
        equilibrium = run_lotka_volterra(
            initial, TIME_SPAN, community, interactions, growth_rates, carrying_capacities
        )
        equilibrium[equilibrium < EXTINCTION_THRESHOLD] = 0
        parents.append(equilibrium)

    if return_interactions:
        return _pairwise_results(parents, library, interactions, growth_rates,
                                 carrying_capacities, n_communities), interactions, library

    return _pairwise_results(parents, library, interactions, growth_rates,
                             carrying_capacities, n_communities)


def _pairwise_results(parents, library, interactions, growth_rates,
                      carrying_capacities, n_communities):
    """Coalesce every pair of parental communities at their equilibria."""
    results = []
    for a in range(n_communities):
        for b in range(a + 1, n_communities):
            mixed = np.array(parents[a] + parents[b])
            survived = mixed > EXTINCTION_THRESHOLD
            coalesced = run_lotka_volterra(
                mixed, TIME_SPAN, survived, interactions, growth_rates, carrying_capacities
            )
            coalesced[coalesced < EXTINCTION_THRESHOLD] = 0
            results.append((parents[a], parents[b], coalesced))
    return results


def mean_interaction_before_after(interactions, library, parents):
    """Mean off-diagonal competition coefficient before and after assembly.

    ``before`` averages over all species pairs within a parental community as
    it was seeded; ``after`` averages over only those that survived to
    equilibrium.  Competitive exclusion removes strongly competing species, so
    the surviving set should show weaker mutual interactions — the assembly
    effect in Fig. 2c.

    Returns one ``(before, after)`` pair per parental community.
    """
    summaries = []
    for community, equilibrium in zip(library, parents):
        seeded = np.flatnonzero(community > 0)
        survived = np.flatnonzero((community > 0) & (equilibrium > EXTINCTION_THRESHOLD))
        if len(seeded) < 2 or len(survived) < 2:
            continue

        def mean_off_diagonal(members):
            block = interactions[np.ix_(members, members)]
            off = block[~np.eye(len(members), dtype=bool)]
            return float(off.mean())

        summaries.append((mean_off_diagonal(seeded), mean_off_diagonal(survived)))
    return summaries


def sweep(mu_values, n_replicates=100, n_communities=2, species_per_community=12,
          spread=0.0, pool_size=None, progress=True):
    """Run the interaction-strength sweep.

    Returns one row per coalescence event with its similarity-map coordinates,
    retention, asymmetry and outcome class.
    """
    rows = []
    for mu in mu_values:
        for replicate in range(n_replicates):
            for parent_a, parent_b, coalesced in simulate_one_replicate(
                mu, replicate, n_communities, species_per_community, spread,
                pool_size,
            ):
                x1, x2, residual = decompose(parent_a, parent_b, coalesced)
                r, asymmetry = retention_and_asymmetry(x1, x2)
                rows.append({
                    "mu": mu,
                    "replicate": replicate,
                    "x1": x1,
                    "x2": x2,
                    "residual": residual,
                    "r": r,
                    "asymmetry": asymmetry,
                    "outcome": classify(r, asymmetry),
                })
        if progress:
            print(f"  mu = {mu:.2f} done ({len(rows)} events so far)")
    return pd.DataFrame(rows)


def outcome_fractions(table):
    """Fraction of each outcome class at each value of mu."""
    counts = table.groupby("mu").outcome.value_counts(normalize=True).unstack(fill_value=0)
    counts.columns = [{0: "Dominance", 1: "Mixture", 2: "Restructuring"}[c]
                      for c in counts.columns]
    counts["n"] = table.groupby("mu").size()
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mu-min", type=float, default=0.0)
    parser.add_argument("--mu-max", type=float, default=1.2)
    parser.add_argument("--mu-step", type=float, default=0.1)
    parser.add_argument("--replicates", type=int, default=100,
                        help="replicate pools per value of mu")
    parser.add_argument("--communities", type=int, default=2,
                        help="parental communities per pool")
    parser.add_argument("--species", type=int, default=12,
                        help="species per parental community")
    parser.add_argument("--spread", type=float, default=0.0)
    parser.add_argument("--pool-size", type=int, default=None,
                        help="species in the pool; defaults to communities x species")
    parser.add_argument("--out", type=Path, default=Path("simulation_outcomes.csv"))
    args = parser.parse_args()

    mu_values = np.arange(args.mu_min, args.mu_max, args.mu_step)
    print(f"sweeping mu over {len(mu_values)} values, "
          f"{args.replicates} replicates, "
          f"{args.communities} x {args.species} species")

    table = sweep(mu_values, args.replicates, args.communities, args.species,
                  args.spread, args.pool_size)
    table.to_csv(args.out, index=False)
    print(f"\nwrote {args.out}  ({len(table)} events)")
    print(outcome_fractions(table).to_string(float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
