# Generalized Lotka-Volterra simulations

`coalescence_simulation.py` is the simulation behind Figures 2 and 3.

```bash
# interaction-strength sweep, manuscript configuration
python pipeline/03_simulation/coalescence_simulation.py \
    --mu-min 0 --mu-max 1.2 --mu-step 0.1 \
    --replicates 200 --communities 4 --species 12 --pool-size 54 \
    --out simulation_outcomes.csv
```

It replaces `gLV_simulation.py`, a demonstration script that stood in for the
real simulation in earlier versions of this repository. That script used ten
species drawn from a single shared pool and did not implement the sampling
design, the interaction-strength sweep, or the outcome classification, so it
could not reproduce any published panel.

## Model

Species grow logistically and compete pairwise:

    dn_i/dt = g_i n_i (1 - sum_j I_ij n_j / k_i)

with growth rates and carrying capacities fixed to 1 and self-interaction
`I_ii = 1`. Off-diagonal competition coefficients are drawn from `U(0, 2*mu)`,
so `mu` sets both the mean and the spread of competitive effects. Species below
`1e-3` are treated as extinct.

## Design

A pool of species is partitioned into non-overlapping parental communities of
equal size. Each is assembled to equilibrium alone; pairs are then mixed at
their equilibrium abundances and re-equilibrated. This mirrors the experimental
protocol of separate parental assembly followed by pairwise coalescence.

The manuscript configuration is four 12-species communities drawn from a pool of
54, giving six pairwise coalescences per replicate; 200 replicate pools
therefore give 1,200 coalescence events per value of `mu`. The pool is larger
than the communities consume, so six species belong to no parental community
while still appearing in the interaction matrix. Pool size affects both the
dynamics and the number of random draws, so `--pool-size` must match to
reproduce a given sweep.

## Reproducibility

Each replicate seeds NumPy's global generator with the replicate index before
drawing anything, so a sweep is fully determined by its parameters — no
simulation outputs need to be archived to reproduce the figures.

The **order** of random draws is load-bearing: the interaction matrix consumes
`N * N` values and the initial abundances `N` more, in that order. Reordering
them changes the results at the same seed.

## Verification

Reproduces the cached notebook outputs exactly. For the two-community,
12-species sweep the similarity-map coordinates agree to five decimal places
across the range of `mu`:

| mu index | replicate | this code | cached |
|---|---|---|---|
| 0 | 0 | (0.70743, 0.70678, 0.00123) | (0.70743, 0.70678, 0.00123) |
| 3 | 0 | (0.68319, 0.54653, 0.48431) | (0.68319, 0.54653, 0.48431) |
| 6 | 11 | (0.22487, 0.57178, 0.78899) | (0.22487, 0.57178, 0.78899) |
| 11 | 99 | (1.00000, 0.00000, 0.00000) | (1.00000, 0.00000, 0.00000) |

Source: `Main_Fig_2/Impt_Generate_Fig2_5.ipynb` and
`Main_Fig_2/Generate_Fig2_1..2_4.ipynb`, whose `InitializeCommunityPool`,
`uniform_distribution` and sweep loop are consolidated here. The species-pool
and dynamics helpers are `coalescence.species_pool` and `coalescence.lv`.
